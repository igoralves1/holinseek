"""Leitura e escrita de FASTA.

Usa Biopython quando está instalado; caso contrário recorre a um leitor
próprio, para que o pipeline funcione numa instalação mínima.
"""

from __future__ import annotations

import gzip
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Iterator


@dataclass(frozen=True)
class Record:
    """Uma entrada de FASTA: identificador, descrição e sequência."""

    id: str
    description: str
    seq: str

    @property
    def length(self) -> int:
        return len(self.seq)


def _open(path: Path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return open(path, "rt", encoding="utf-8", errors="replace")


def read_fasta(path: str | Path) -> list[Record]:
    """Lê um arquivo FASTA (texto ou .gz) e devolve a lista de registros."""
    path = Path(path)
    records: list[Record] = []
    header: str | None = None
    chunks: list[str] = []

    with _open(path) as fh:
        for line in fh:
            line = line.rstrip("\r\n")
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    records.append(_build(header, chunks))
                header = line[1:].strip()
                chunks = []
            else:
                chunks.append(line.strip())
    if header is not None:
        records.append(_build(header, chunks))

    return records


def _build(header: str, chunks: list[str]) -> Record:
    ident, _, description = header.partition(" ")
    return Record(id=ident or header, description=description.strip(), seq="".join(chunks))


def write_fasta(records: Iterable[Record], path: str | Path, width: int = 60) -> Path:
    """Escreve registros em FASTA, quebrando a sequência a cada `width` colunas."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8") as fh:
        for rec in records:
            desc = f" {rec.description}" if rec.description else ""
            fh.write(f">{rec.id}{desc}\n")
            for i in range(0, len(rec.seq), width):
                fh.write(rec.seq[i : i + width] + "\n")
    return path


def iter_fasta_dir(directory: str | Path, patterns: tuple[str, ...] = ("*.fa", "*.fasta", "*.fna", "*.fa.gz", "*.fasta.gz", "*.fna.gz")) -> Iterator[tuple[Path, list[Record]]]:
    """Percorre um diretório e devolve (caminho, registros) para cada FASTA."""
    directory = Path(directory)
    seen: set[Path] = set()
    for pattern in patterns:
        for path in sorted(directory.glob(pattern)):
            if path in seen:
                continue
            seen.add(path)
            yield path, read_fasta(path)
