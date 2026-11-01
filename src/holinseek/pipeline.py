"""Orquestração: do genoma bruto à lista de candidatos priorizada."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from .biophysics import compute_features
from .cassette import analyse_neighbourhood
from .config import Settings
from .interpret import Interpretation, interpret_candidate, narrative_summary
from .orfs import predict_genes
from .report import write_json, write_report, write_table
from .scoring import Candidate, evaluate, rank
from .seqio import Record, iter_fasta_dir, read_fasta


@dataclass
class PipelineResult:
    genomes: int
    orfs: int
    candidates: list[Candidate]
    interpretations: list[Interpretation]
    gene_method: str
    outputs: dict[str, Path]
    narrative_status: str


def _load_genomes(path: Path) -> list[Record]:
    if path.is_dir():
        records: list[Record] = []
        for _file, recs in iter_fasta_dir(path):
            records.extend(recs)
        return records
    return read_fasta(path)


def run(
    input_path: str | Path,
    settings: Settings,
    min_score: float = 0.5,
    window_bp: int = 6000,
    prefer_prodigal: bool = True,
    max_candidates: int = 200,
    on_progress=None,
) -> PipelineResult:
    """Roda o pipeline inteiro sobre um FASTA ou um diretório de FASTAs."""
    def say(msg: str) -> None:
        if on_progress:
            on_progress(msg)

    settings.ensure_dirs()
    input_path = Path(input_path)

    say("lendo genomas")
    genomes = _load_genomes(input_path)
    if not genomes:
        raise ValueError(f"nenhuma sequência encontrada em {input_path}")

    say(f"predizendo genes em {len(genomes)} genoma(s)")
    orfs, gene_method = predict_genes(
        genomes,
        min_aa=max(20, settings.criteria.min_length - 20),
        max_aa=600,
        prefer_prodigal=prefer_prodigal,
    )

    # Anotações vindas do cabeçalho do FASTA, quando houver.
    annotations = {g.id: g.description for g in genomes}

    say(f"avaliando {len(orfs)} ORFs")
    crit = settings.criteria
    shortlist = [
        o for o in orfs if crit.min_length <= o.length <= crit.max_length
    ]

    candidates: list[Candidate] = []
    for orf in shortlist:
        features = compute_features(
            orf.protein, window=crit.tm_window, threshold=crit.tm_threshold
        )
        # Só vale examinar a vizinhança de quem já é plausível pela física.
        if features.n_tmd == 0 and features.gravy < crit.min_gravy:
            continue
        cassette = analyse_neighbourhood(orf, orfs, window_bp=window_bp, annotations=annotations)
        candidates.append(
            evaluate(
                protein_id=orf.protein_id,
                genome_id=orf.genome_id,
                sequence=orf.protein,
                features=features,
                cassette=cassette,
                criteria=crit,
            )
        )

    say("ordenando candidatos")
    ranked = rank(candidates, min_score=min_score)[:max_candidates]

    tm_method = "pyrodigal + hidropatia de Kyte-Doolittle" if gene_method == "pyrodigal" else "hidropatia de Kyte-Doolittle"
    interpretations = [interpret_candidate(c, tm_method=tm_method) for c in ranked]

    say("montando o relatório")
    narrative, narrative_status = narrative_summary(ranked, interpretations, settings.llm)

    meta = {
        "entrada": str(input_path),
        "genomas": len(genomes),
        "ORFs preditas": len(orfs),
        "método de predição de genes": gene_method,
        "método de topologia": tm_method,
        "nota mínima": min_score,
        "janela de sintenia": f"{window_bp} pb",
        "critérios": (
            f"{crit.min_length}–{crit.max_length} aa · {crit.min_tmd}–{crit.max_tmd} TMD · "
            f"GRAVY > {crit.min_gravy} · cauda ≥ +{crit.min_tail_charge}"
        ),
        "camada de prosa": narrative_status,
    }

    results_dir = settings.results_dir
    outputs = {
        "tabela": write_table(ranked, results_dir / "candidatos.tsv"),
        "json": write_json(ranked, results_dir / "candidatos.json", meta),
        "relatorio": write_report(
            ranked, interpretations, results_dir / "relatorio.md", meta, narrative, narrative_status
        ),
    }

    return PipelineResult(
        genomes=len(genomes),
        orfs=len(orfs),
        candidates=ranked,
        interpretations=interpretations,
        gene_method=gene_method,
        outputs=outputs,
        narrative_status=narrative_status,
    )
