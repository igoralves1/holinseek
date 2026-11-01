"""Aprendizado de máquina sobre os candidatos.

Duas coisas distintas acontecem aqui, e vale separá-las bem:

**Não supervisionado — sempre disponível.** Agrupa os candidatos no espaço de
atributos físico-químicos. Não precisa de rótulo nenhum e serve para enxergar
estrutura: candidatos que se parecem entre si caem no mesmo grupo, e um
candidato isolado de todos os outros merece um olhar.

**Supervisionado — exige rótulos seus.** Treina um classificador a partir de
um CSV de exemplos que você mesmo monta e valida. O projeto não embarca
rótulos prontos de propósito: um conjunto de treino herdado sem auditoria é a
forma mais rápida de um classificador aprender a base de dados em vez da
biologia.

A validação cruzada é **agrupada por genoma**. Sem isso, duas holinas do mesmo
fago caem uma no treino e outra no teste, e a acurácia medida vira ficção.
"""

from __future__ import annotations

import csv
import json
from dataclasses import dataclass
from pathlib import Path

from .scoring import Candidate

FEATURE_NAMES = (
    "length",
    "gravy",
    "n_tmd",
    "tm_fraction",
    "max_tm_peak",
    "tail_length",
    "tail_charge",
    "tail_charge_density",
    "net_charge_ph7",
    "isoelectric_point",
    "aromaticity",
)


class SklearnMissing(RuntimeError):
    """scikit-learn não está instalado."""


def _require_sklearn():
    try:
        import sklearn  # noqa: F401
    except ImportError as exc:  # pragma: no cover
        raise SklearnMissing(
            "scikit-learn não instalado — rode: pip install -r requirements.txt"
        ) from exc


def feature_matrix(candidates: list[Candidate]) -> tuple[list[list[float]], list[str]]:
    """Monta a matriz de atributos e a lista de identificadores, na mesma ordem."""
    rows: list[list[float]] = []
    ids: list[str] = []
    for c in candidates:
        d = c.features.as_dict()
        rows.append([float(d[name]) for name in FEATURE_NAMES])
        ids.append(c.protein_id)
    return rows, ids


@dataclass
class ClusterResult:
    labels: dict[str, int]
    n_clusters: int
    silhouette: float | None
    outliers: list[str]
    coords: dict[str, tuple[float, float]]

    def to_json(self) -> str:
        return json.dumps(
            {
                "n_clusters": self.n_clusters,
                "silhouette": self.silhouette,
                "labels": self.labels,
                "outliers": self.outliers,
                "coords": {k: list(v) for k, v in self.coords.items()},
            },
            indent=2,
            ensure_ascii=False,
        )


def cluster(candidates: list[Candidate], n_clusters: int = 0) -> ClusterResult:
    """Agrupa os candidatos por similaridade físico-química.

    Com `n_clusters=0`, escolhe o número de grupos pelo melhor coeficiente de
    silhueta entre 2 e 6. Devolve também coordenadas em duas dimensões (PCA),
    úteis para a figura do relatório.
    """
    _require_sklearn()
    from sklearn.cluster import KMeans
    from sklearn.decomposition import PCA
    from sklearn.metrics import silhouette_score
    from sklearn.preprocessing import StandardScaler

    rows, ids = feature_matrix(candidates)
    if len(rows) < 4:
        return ClusterResult({i: 0 for i in ids}, 1, None, [], {})

    X = StandardScaler().fit_transform(rows)

    best_k, best_score, best_labels = 1, -1.0, [0] * len(rows)
    candidate_ks = [n_clusters] if n_clusters >= 2 else range(2, min(7, len(rows)))
    for k in candidate_ks:
        if k >= len(rows):
            continue
        labels = KMeans(n_clusters=k, n_init=10, random_state=0).fit_predict(X)
        try:
            s = silhouette_score(X, labels)
        except ValueError:
            continue
        if s > best_score:
            best_k, best_score, best_labels = k, float(s), labels.tolist()

    coords2d = PCA(n_components=2, random_state=0).fit_transform(X)

    # Um candidato isolado é aquele cujo grupo tem um único membro.
    counts: dict[int, int] = {}
    for lab in best_labels:
        counts[lab] = counts.get(lab, 0) + 1
    outliers = [ids[i] for i, lab in enumerate(best_labels) if counts[lab] == 1]

    return ClusterResult(
        labels=dict(zip(ids, best_labels)),
        n_clusters=best_k,
        silhouette=round(best_score, 4) if best_score > -1 else None,
        outliers=outliers,
        coords={ids[i]: (round(float(coords2d[i][0]), 4), round(float(coords2d[i][1]), 4)) for i in range(len(ids))},
    )


def load_labels(path: str | Path) -> dict[str, int]:
    """Lê um CSV de rótulos com as colunas protein_id,label (1 = holina, 0 = não)."""
    labels: dict[str, int] = {}
    with open(path, newline="", encoding="utf-8") as fh:
        for row in csv.DictReader(fh):
            pid = (row.get("protein_id") or "").strip()
            raw = (row.get("label") or "").strip().lower()
            if not pid or not raw:
                continue
            labels[pid] = 1 if raw in {"1", "true", "sim", "yes", "holin", "holina"} else 0
    return labels


@dataclass
class TrainResult:
    n_samples: int
    n_positive: int
    cv_folds: int
    mean_auc: float | None
    mean_precision: float | None
    mean_recall: float | None
    feature_importance: dict[str, float]
    model_path: str
    warning: str | None = None


def train_classifier(
    candidates: list[Candidate],
    labels: dict[str, int],
    model_path: str | Path = "models/holin_classifier.joblib",
) -> TrainResult:
    """Treina um classificador, com validação cruzada agrupada por genoma."""
    _require_sklearn()
    import joblib
    import numpy as np
    from sklearn.ensemble import RandomForestClassifier
    from sklearn.model_selection import GroupKFold, cross_validate
    from sklearn.pipeline import Pipeline
    from sklearn.preprocessing import StandardScaler

    labelled = [c for c in candidates if c.protein_id in labels]
    if len(labelled) < 10:
        raise ValueError(
            f"apenas {len(labelled)} exemplos rotulados — são necessários ao menos 10 "
            "para que a validação cruzada diga alguma coisa"
        )

    rows, ids = feature_matrix(labelled)
    X = np.array(rows)
    y = np.array([labels[i] for i in ids])
    groups = np.array([c.genome_id for c in labelled])

    n_positive = int(y.sum())
    warning = None
    if n_positive < 5 or (len(y) - n_positive) < 5:
        warning = (
            f"conjunto desequilibrado ({n_positive} positivos de {len(y)}); "
            "as métricas abaixo são indicativas, não conclusivas"
        )

    model = Pipeline(
        [
            ("scale", StandardScaler()),
            (
                "forest",
                RandomForestClassifier(
                    n_estimators=400,
                    min_samples_leaf=2,
                    class_weight="balanced",
                    random_state=0,
                    n_jobs=-1,
                ),
            ),
        ]
    )

    n_groups = len(set(groups))
    folds = max(2, min(5, n_groups))
    metrics: dict[str, float | None] = {"roc_auc": None, "precision": None, "recall": None}
    if n_groups >= 2 and 0 < n_positive < len(y):
        scores = cross_validate(
            model,
            X,
            y,
            groups=groups,
            cv=GroupKFold(n_splits=folds),
            scoring=("roc_auc", "precision", "recall"),
            error_score="raise",
        )
        for key in metrics:
            metrics[key] = round(float(np.nanmean(scores[f"test_{key}"])), 4)
    else:
        warning = (warning or "") + " | genomas insuficientes para validação agrupada"

    model.fit(X, y)
    importance = dict(
        sorted(
            zip(FEATURE_NAMES, (round(float(v), 4) for v in model.named_steps["forest"].feature_importances_)),
            key=lambda kv: kv[1],
            reverse=True,
        )
    )

    model_path = Path(model_path)
    model_path.parent.mkdir(parents=True, exist_ok=True)
    joblib.dump({"model": model, "features": list(FEATURE_NAMES)}, model_path)

    return TrainResult(
        n_samples=len(y),
        n_positive=n_positive,
        cv_folds=folds,
        mean_auc=metrics["roc_auc"],
        mean_precision=metrics["precision"],
        mean_recall=metrics["recall"],
        feature_importance=importance,
        model_path=str(model_path),
        warning=warning.strip(" |") if warning else None,
    )


def predict_proba(candidates: list[Candidate], model_path: str | Path) -> dict[str, float]:
    """Aplica um classificador treinado e devolve a probabilidade por candidato."""
    _require_sklearn()
    import joblib
    import numpy as np

    bundle = joblib.load(model_path)
    model = bundle["model"]
    rows, ids = feature_matrix(candidates)
    probs = model.predict_proba(np.array(rows))[:, 1]
    return {pid: round(float(p), 4) for pid, p in zip(ids, probs)}
