"""Download de genomas de bacteriófago das bases públicas do NCBI.

Usa as E-utilities (esearch + efetch), a interface programática oficial do
NCBI. Duas regras da casa, respeitadas aqui:

  * toda requisição identifica um e-mail de contato (`NCBI_EMAIL` no .env);
  * sem chave de API, no máximo 3 requisições por segundo; com chave, 10.

O módulo não embarca nenhum conjunto de genomas: tudo é baixado na hora, da
fonte primária, e fica em `data/raw/` com o registro da consulta que o gerou.
"""

from __future__ import annotations

import json
import time
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from urllib.parse import urlencode

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

# Consultas prontas, úteis como ponto de partida.
PRESET_QUERIES: dict[str, str] = {
    "fagos-refseq": (
        '"Caudoviricetes"[Organism] AND srcdb_refseq[PROP] AND '
        '("complete genome"[Title] OR "complete sequence"[Title])'
    ),
    "fagos-ecoli": (
        'bacteriophage[Title] AND "Escherichia coli"[Title] AND '
        'srcdb_refseq[PROP] AND "complete genome"[Title]'
    ),
    "fagos-staphylococcus": (
        'bacteriophage[Title] AND Staphylococcus[Title] AND '
        'srcdb_refseq[PROP] AND "complete genome"[Title]'
    ),
    "fagos-salmonella": (
        'bacteriophage[Title] AND Salmonella[Title] AND '
        'srcdb_refseq[PROP] AND "complete genome"[Title]'
    ),
}


class NCBIError(RuntimeError):
    pass


@dataclass
class DownloadResult:
    query: str
    accessions: list[str]
    fasta_path: Path
    manifest_path: Path
    retrieved_at: str


def _request(endpoint: str, params: dict, email: str, api_key: str, timeout: int = 60) -> str:
    try:
        import requests
    except ImportError as exc:  # pragma: no cover
        raise NCBIError("pacote 'requests' não instalado") from exc

    if not email:
        raise NCBIError(
            "NCBI_EMAIL não definido. O NCBI exige um e-mail de contato em toda "
            "requisição programática — preencha no arquivo .env."
        )

    params = {**params, "email": email, "tool": "holinseek"}
    if api_key:
        params["api_key"] = api_key

    url = f"{EUTILS}/{endpoint}?{urlencode(params)}"
    resp = requests.get(url, timeout=timeout)
    if resp.status_code >= 400:
        raise NCBIError(f"NCBI respondeu {resp.status_code}: {resp.text[:200]}")
    return resp.text


def search(
    query: str,
    email: str,
    api_key: str = "",
    retmax: int = 20,
    db: str = "nuccore",
) -> list[str]:
    """Busca acessos que satisfazem a consulta. Devolve a lista de identificadores."""
    raw = _request(
        "esearch.fcgi",
        {"db": db, "term": query, "retmax": retmax, "retmode": "json"},
        email,
        api_key,
    )
    try:
        data = json.loads(raw)
        return list(data["esearchresult"]["idlist"])
    except (ValueError, KeyError) as exc:
        raise NCBIError(f"resposta inesperada do esearch: {raw[:200]}") from exc


def fetch_fasta(
    accessions: list[str],
    email: str,
    api_key: str = "",
    db: str = "nuccore",
    batch_size: int = 20,
) -> str:
    """Baixa as sequências em FASTA, em lotes, respeitando o limite de taxa."""
    if not accessions:
        return ""

    delay = 0.12 if api_key else 0.4
    chunks: list[str] = []

    for i in range(0, len(accessions), batch_size):
        batch = accessions[i : i + batch_size]
        chunks.append(
            _request(
                "efetch.fcgi",
                {"db": db, "id": ",".join(batch), "rettype": "fasta", "retmode": "text"},
                email,
                api_key,
            )
        )
        time.sleep(delay)

    return "\n".join(c.strip() for c in chunks if c.strip()) + "\n"


def download(
    query: str,
    out_dir: str | Path,
    email: str,
    api_key: str = "",
    limit: int = 20,
    name: str = "genomas",
) -> DownloadResult:
    """Busca, baixa e grava genomas junto com o registro da consulta."""
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    accessions = search(query, email, api_key, retmax=limit)
    if not accessions:
        raise NCBIError(f"nenhum resultado para a consulta: {query}")

    fasta = fetch_fasta(accessions, email, api_key)
    fasta_path = out_dir / f"{name}.fasta"
    fasta_path.write_text(fasta, encoding="utf-8")

    retrieved_at = datetime.now(timezone.utc).isoformat(timespec="seconds")
    manifest = {
        "query": query,
        "database": "nuccore",
        "accessions": accessions,
        "count": len(accessions),
        "retrieved_at": retrieved_at,
        "source": "NCBI E-utilities",
        "note": (
            "Registro de proveniência: guarde junto com os resultados para que a "
            "análise possa ser refeita exatamente sobre o mesmo conjunto."
        ),
    }
    manifest_path = out_dir / f"{name}.manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2, ensure_ascii=False), encoding="utf-8")

    return DownloadResult(query, accessions, fasta_path, manifest_path, retrieved_at)
