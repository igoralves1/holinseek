"""Linha de comando do HolinSeek."""

from __future__ import annotations

from pathlib import Path

import typer
from rich.console import Console
from rich.panel import Panel
from rich.table import Table

from .config import load_settings
from .llm import describe
from .ml import SklearnMissing, cluster, load_labels, predict_proba, train_classifier
from .ncbi import PRESET_QUERIES, NCBIError, download
from .pipeline import run as run_pipeline

app = typer.Typer(
    add_completion=False,
    help="Descoberta computacional de holinas em genomas de bacteriófago.",
)
console = Console()


@app.command("check")
def check_cmd() -> None:
    """Mostra a configuração ativa e o que está instalado."""
    settings = load_settings()

    table = Table(title="Configuração", show_header=False, box=None)
    table.add_row("e-mail NCBI", settings.ncbi_email or "[red]não definido[/red]")
    table.add_row("chave NCBI", "definida" if settings.ncbi_api_key else "ausente (3 req/s)")
    table.add_row("diretório de dados", str(settings.data_dir))
    table.add_row("diretório de resultados", str(settings.results_dir))
    c = settings.criteria
    table.add_row("comprimento", f"{c.min_length}–{c.max_length} aa")
    table.add_row("hélices TM", f"{c.min_tmd}–{c.max_tmd}")
    table.add_row("GRAVY mínimo", str(c.min_gravy))
    table.add_row("carga da cauda", f"≥ +{c.min_tail_charge}")
    table.add_row("camada de prosa", describe(settings.llm))
    console.print(table)

    deps = Table(title="Bibliotecas", show_header=False, box=None)
    for name, label in (
        ("numpy", "numérico"),
        ("sklearn", "aprendizado de máquina"),
        ("requests", "downloads"),
        ("Bio", "biopython"),
        ("pyrodigal", "predição de genes (opcional)"),
        ("pyhmmer", "perfis de domínio (opcional)"),
        ("torch", "embeddings (opcional)"),
        ("anthropic", "provedor de prosa (opcional)"),
    ):
        try:
            __import__(name)
            deps.add_row(f"[green]ok[/green]  {name}", label)
        except ImportError:
            deps.add_row(f"[yellow]--[/yellow]  {name}", label)
    console.print(deps)


@app.command("fetch")
def fetch_cmd(
    query: str = typer.Option(
        "fagos-refseq",
        "--query",
        "-q",
        help="Consulta ao NCBI, ou o nome de uma consulta pronta.",
    ),
    limit: int = typer.Option(20, "--limit", "-n", help="Quantos genomas baixar."),
    name: str = typer.Option("genomas", "--name", help="Nome base do arquivo de saída."),
    listar: bool = typer.Option(False, "--list", help="Lista as consultas prontas e sai."),
) -> None:
    """Baixa genomas de bacteriófago do NCBI."""
    if listar:
        for key, value in PRESET_QUERIES.items():
            console.print(f"[bold]{key}[/bold]\n  {value}\n")
        raise typer.Exit()

    settings = load_settings()
    settings.ensure_dirs()
    term = PRESET_QUERIES.get(query, query)

    try:
        with console.status(f"consultando o NCBI: {term[:60]}…"):
            result = download(
                query=term,
                out_dir=settings.raw_dir,
                email=settings.ncbi_email,
                api_key=settings.ncbi_api_key,
                limit=limit,
                name=name,
            )
    except NCBIError as exc:
        console.print(f"[red]{exc}[/red]")
        raise typer.Exit(code=1)

    console.print(
        Panel(
            f"{len(result.accessions)} genoma(s)\n"
            f"FASTA:     {result.fasta_path}\n"
            f"manifesto: {result.manifest_path}",
            title="download concluído",
        )
    )


@app.command("run")
def run_cmd(
    entrada: Path = typer.Argument(..., help="Arquivo FASTA ou diretório de FASTAs."),
    min_score: float = typer.Option(0.5, "--min-score", help="Nota mínima para entrar na lista."),
    window: int = typer.Option(6000, "--window", help="Janela de sintenia, em pares de base."),
    top: int = typer.Option(20, "--top", help="Quantos candidatos mostrar no terminal."),
    sem_prodigal: bool = typer.Option(False, "--sem-prodigal", help="Força o localizador de ORFs interno."),
) -> None:
    """Roda o pipeline completo sobre genomas locais."""
    settings = load_settings()

    with console.status("rodando…") as status:
        result = run_pipeline(
            entrada,
            settings,
            min_score=min_score,
            window_bp=window,
            prefer_prodigal=not sem_prodigal,
            on_progress=lambda m: status.update(m),
        )

    table = Table(title=f"{len(result.candidates)} candidato(s) acima do corte")
    for col in ("#", "proteína", "genoma", "nota", "faixa", "aa", "TMD", "GRAVY", "cauda", "endolisina"):
        table.add_column(col)
    for i, c in enumerate(result.candidates[:top], 1):
        f = c.features
        table.add_row(
            str(i),
            c.protein_id,
            c.genome_id[:22],
            f"{c.score:.2f}",
            c.tier,
            str(f.length),
            str(f.n_tmd),
            f"{f.gravy:+.2f}",
            f"{f.tail_charge:+d}",
            "sim" if c.cassette.endolysin_nearby else "não",
        )
    console.print(table)

    console.print(
        Panel(
            f"genomas: {result.genomes}   ORFs: {result.orfs}   método: {result.gene_method}\n"
            f"tabela:    {result.outputs['tabela']}\n"
            f"json:      {result.outputs['json']}\n"
            f"relatório: {result.outputs['relatorio']}\n\n"
            f"[dim]{result.narrative_status}[/dim]",
            title="saídas",
        )
    )
    console.print(
        "[yellow]Lembrete:[/yellow] tudo acima é predição. Nenhum candidato é uma "
        "holina até a confirmação experimental."
    )


@app.command("cluster")
def cluster_cmd(
    entrada: Path = typer.Argument(..., help="Arquivo FASTA ou diretório de FASTAs."),
    grupos: int = typer.Option(0, "--grupos", "-k", help="Número de grupos (0 = escolher automaticamente)."),
) -> None:
    """Agrupa os candidatos por similaridade físico-química."""
    settings = load_settings()
    result = run_pipeline(entrada, settings, on_progress=lambda m: None)

    try:
        cl = cluster(result.candidates, n_clusters=grupos)
    except SklearnMissing as exc:
        console.print(f"[red]{exc}[/red]")
        raise typer.Exit(code=1)

    out = settings.results_dir / "grupos.json"
    out.write_text(cl.to_json(), encoding="utf-8")

    console.print(
        Panel(
            f"grupos: {cl.n_clusters}\n"
            f"silhueta: {cl.silhouette if cl.silhouette is not None else 'n/d'}\n"
            f"isolados: {', '.join(cl.outliers) if cl.outliers else 'nenhum'}\n"
            f"arquivo: {out}",
            title="agrupamento",
        )
    )


@app.command("train")
def train_cmd(
    entrada: Path = typer.Argument(..., help="Arquivo FASTA ou diretório de FASTAs."),
    rotulos: Path = typer.Option(..., "--rotulos", "-r", help="CSV com protein_id,label."),
    saida: Path = typer.Option(Path("models/holin_classifier.joblib"), "--saida", "-o"),
) -> None:
    """Treina um classificador a partir dos seus próprios rótulos."""
    settings = load_settings()
    result = run_pipeline(entrada, settings, min_score=0.0, on_progress=lambda m: None)

    try:
        res = train_classifier(result.candidates, load_labels(rotulos), model_path=saida)
    except (SklearnMissing, ValueError) as exc:
        console.print(f"[red]{exc}[/red]")
        raise typer.Exit(code=1)

    table = Table(title="treino", show_header=False, box=None)
    table.add_row("exemplos", f"{res.n_samples} ({res.n_positive} positivos)")
    table.add_row("validação", f"GroupKFold por genoma, {res.cv_folds} dobras")
    table.add_row("AUC médio", str(res.mean_auc))
    table.add_row("precisão média", str(res.mean_precision))
    table.add_row("revocação média", str(res.mean_recall))
    table.add_row("modelo", res.model_path)
    console.print(table)

    imp = Table(title="importância dos atributos")
    imp.add_column("atributo")
    imp.add_column("peso")
    for k, v in list(res.feature_importance.items())[:8]:
        imp.add_row(k, f"{v:.4f}")
    console.print(imp)

    if res.warning:
        console.print(f"[yellow]{res.warning}[/yellow]")


@app.command("predict")
def predict_cmd(
    entrada: Path = typer.Argument(..., help="Arquivo FASTA ou diretório de FASTAs."),
    modelo: Path = typer.Option(Path("models/holin_classifier.joblib"), "--modelo", "-m"),
    top: int = typer.Option(20, "--top"),
) -> None:
    """Aplica um classificador já treinado."""
    settings = load_settings()
    result = run_pipeline(entrada, settings, min_score=0.0, on_progress=lambda m: None)

    try:
        probs = predict_proba(result.candidates, modelo)
    except (SklearnMissing, FileNotFoundError) as exc:
        console.print(f"[red]{exc}[/red]")
        raise typer.Exit(code=1)

    table = Table(title="probabilidade pelo classificador")
    table.add_column("proteína")
    table.add_column("probabilidade")
    table.add_column("nota por regras")
    for pid, p in sorted(probs.items(), key=lambda kv: kv[1], reverse=True)[:top]:
        rule = next((c.score for c in result.candidates if c.protein_id == pid), 0.0)
        table.add_row(pid, f"{p:.3f}", f"{rule:.2f}")
    console.print(table)


@app.command("demo")
def demo_cmd() -> None:
    """Roda o pipeline sobre o genoma sintético de exemplo, sem rede."""
    exemplo = Path(__file__).resolve().parents[2] / "data" / "exemplo" / "fago_sintetico.fasta"
    if not exemplo.exists():
        console.print(f"[red]arquivo de exemplo não encontrado em {exemplo}[/red]")
        raise typer.Exit(code=1)
    console.print(
        Panel(
            "Este é um genoma SINTÉTICO, construído para exercitar o pipeline.\n"
            "Não é um fago real e não serve como resultado biológico.",
            title="demonstração",
        )
    )
    run_cmd(entrada=exemplo, min_score=0.4, window=6000, top=10, sem_prodigal=False)


if __name__ == "__main__":
    app()
