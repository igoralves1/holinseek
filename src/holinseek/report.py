"""Saídas do pipeline: tabela, JSON e relatório legível."""

from __future__ import annotations

import csv
import json
from datetime import datetime, timezone
from pathlib import Path

from .interpret import Interpretation
from .scoring import Candidate


def write_table(candidates: list[Candidate], path: str | Path) -> Path:
    """Grava a tabela de candidatos em TSV, uma linha por proteína."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    if not candidates:
        path.write_text("sem candidatos\n", encoding="utf-8")
        return path

    rows = [c.as_row() for c in candidates]
    with open(path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    return path


def write_json(candidates: list[Candidate], path: str | Path, meta: dict) -> Path:
    """Grava o resultado completo, com os vereditos critério a critério."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "meta": meta,
        "candidates": [
            {
                "protein_id": c.protein_id,
                "genome_id": c.genome_id,
                "score": c.score,
                "tier": c.tier,
                "features": c.features.as_dict(),
                "checks": [chk.as_dict() for chk in c.checks],
                "cassette": {
                    "endolysin_nearby": c.cassette.endolysin_nearby,
                    "endolysin_id": c.cassette.endolysin_id,
                    "endolysin_distance_bp": c.cassette.endolysin_distance,
                    "endolysin_reason": c.cassette.endolysin_reason,
                    "spanin_nearby": c.cassette.spanin_nearby,
                    "overlapping_orf": c.cassette.overlapping_orf,
                    "neighbours_examined": c.cassette.neighbour_count,
                },
                "sequence": c.sequence,
            }
            for c in candidates
        ],
    }
    path.write_text(json.dumps(payload, indent=2, ensure_ascii=False), encoding="utf-8")
    return path


def write_report(
    candidates: list[Candidate],
    interpretations: list[Interpretation],
    path: str | Path,
    meta: dict,
    narrative: str = "",
    narrative_status: str = "",
) -> Path:
    """Escreve o relatório em Markdown, pronto para leitura humana."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)

    tiers = {"alta": 0, "média": 0, "baixa": 0, "descartado": 0}
    for c in candidates:
        tiers[c.tier] = tiers.get(c.tier, 0) + 1

    lines = [
        "# HolinSeek — relatório de triagem",
        "",
        f"Gerado em {datetime.now(timezone.utc).isoformat(timespec='seconds')} (UTC)",
        "",
        "## Como este resultado foi produzido",
        "",
    ]
    for key, value in meta.items():
        lines.append(f"- **{key}**: {value}")

    lines += [
        "",
        "## Resumo",
        "",
        f"- candidatos acima do corte: **{len(candidates)}**",
        f"- prioridade alta: {tiers['alta']} · média: {tiers['média']} · baixa: {tiers['baixa']}",
        "",
        "> Tudo abaixo é **predição computacional**. Nenhum candidato é uma holina "
        "até que a bancada mostre que é.",
        "",
    ]

    if candidates:
        lines += [
            "## Tabela",
            "",
            "| # | proteína | genoma | nota | faixa | aa | TMD | GRAVY | cauda | endolisina |",
            "|---|---|---|---|---|---|---|---|---|---|",
        ]
        for i, c in enumerate(candidates, 1):
            f = c.features
            lines.append(
                f"| {i} | `{c.protein_id}` | {c.genome_id} | {c.score:.2f} | {c.tier} | "
                f"{f.length} | {f.n_tmd} | {f.gravy:+.2f} | {f.tail_charge:+d} | "
                f"{'sim' if c.cassette.endolysin_nearby else 'não'} |"
            )
        lines.append("")

    if narrative:
        lines += ["## Leitura do conjunto", "", narrative, "", f"*{narrative_status}*", ""]
    elif narrative_status:
        lines += [f"*{narrative_status}*", ""]

    if interpretations:
        lines += ["## Candidato a candidato", ""]
        for interp in interpretations:
            lines.append(interp.to_text())
            lines.append("")

    path.write_text("\n".join(lines), encoding="utf-8")
    return path
