"""Predição de genes em genomas de bacteriófago.

Usa Prodigal (via pyrodigal, em modo metagenômico/viral) quando disponível.
Sem ele, recorre a um localizador de ORFs próprio, de seis quadros, com o
código genético 11 (bacteriano e de arqueias) — o mesmo usado por fagos.

O localizador interno é mais permissivo que o Prodigal: encontra mais ORFs
espúrias. Para a triagem isso é aceitável, porque os filtros físico-químicos
a jusante descartam a maioria. Para um resultado publicável, instale pyrodigal.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable

from .seqio import Record

# Código genético 11. Difere do padrão apenas nos códons de início alternativos.
CODON_TABLE: dict[str, str] = {}
_BASES = "TCAG"
_AAS = (
    "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
)
for _i, _b1 in enumerate(_BASES):
    for _j, _b2 in enumerate(_BASES):
        for _k, _b3 in enumerate(_BASES):
            CODON_TABLE[_b1 + _b2 + _b3] = _AAS[_i * 16 + _j * 4 + _k]

START_CODONS = ("ATG", "GTG", "TTG")
STOP_CODONS = ("TAA", "TAG", "TGA")
COMPLEMENT = str.maketrans("ACGTRYKMBDHVNacgtrykmbdhvn", "TGCAYRMKVHDBNtgcayrmkvhdbn")


@dataclass(frozen=True)
class ORF:
    """Uma fase de leitura aberta e sua tradução."""

    genome_id: str
    index: int
    start: int          # 1-based, na fita direta do genoma
    end: int            # 1-based, inclusivo
    strand: int         # +1 ou -1
    protein: str

    @property
    def protein_id(self) -> str:
        return f"{self.genome_id}_orf{self.index:04d}"

    @property
    def length(self) -> int:
        return len(self.protein)

    @property
    def midpoint(self) -> float:
        return (self.start + self.end) / 2


def reverse_complement(dna: str) -> str:
    return dna.translate(COMPLEMENT)[::-1]


def translate(dna: str) -> str:
    """Traduz DNA em proteína pelo código 11, parando no primeiro códon de parada."""
    dna = dna.upper().replace("U", "T")
    out: list[str] = []
    for i in range(0, len(dna) - 2, 3):
        aa = CODON_TABLE.get(dna[i : i + 3], "X")
        if aa == "*":
            break
        out.append(aa)
    return "".join(out)


def translate_orf(dna: str) -> str:
    """Traduz uma ORF completa, tratando o códon de início.

    No código 11, GTG e TTG também iniciam tradução. Na posição inicial eles
    carregam formil-metionina, não valina ou leucina — por isso o primeiro
    resíduo vira sempre M. É o que o Prodigal faz, e sem isso a verificação de
    metionina inicial rejeitaria proteínas perfeitamente válidas: vários
    cassetes de lise de fago começam justamente em GTG ou TTG.
    """
    protein = translate(dna)
    if protein and dna[:3].upper() in START_CODONS:
        return "M" + protein[1:]
    return protein


def find_orfs(
    genome: Record,
    min_aa: int = 30,
    max_aa: int = 400,
) -> list[ORF]:
    """Encontra ORFs nos seis quadros de leitura.

    Para cada códon de parada, caminha para trás até o códon de início mais a
    montante dentro do mesmo quadro — a convenção que o Prodigal também adota
    para a ORF mais longa.
    """
    seq = genome.seq.upper().replace("U", "T")
    results: list[ORF] = []
    index = 0

    for strand in (1, -1):
        work = seq if strand == 1 else reverse_complement(seq)
        n = len(work)

        for frame in range(3):
            codons = [work[i : i + 3] for i in range(frame, n - 2, 3)]
            start_pos: int | None = None

            for ci, codon in enumerate(codons):
                if start_pos is None and codon in START_CODONS:
                    start_pos = ci
                elif codon in STOP_CODONS:
                    if start_pos is not None:
                        protein = translate_orf("".join(codons[start_pos:ci]))
                        if min_aa <= len(protein) <= max_aa:
                            a = frame + start_pos * 3
                            b = frame + ci * 3 + 2
                            if strand == 1:
                                gstart, gend = a + 1, b + 1
                            else:
                                gstart, gend = n - b, n - a
                            index += 1
                            results.append(
                                ORF(
                                    genome_id=genome.id,
                                    index=index,
                                    start=gstart,
                                    end=gend,
                                    strand=strand,
                                    protein=protein,
                                )
                            )
                    start_pos = None

    results.sort(key=lambda o: (o.start, o.end))
    return [
        ORF(o.genome_id, i + 1, o.start, o.end, o.strand, o.protein)
        for i, o in enumerate(results)
    ]


def find_orfs_pyrodigal(genome: Record, min_aa: int = 30, max_aa: int = 400) -> list[ORF] | None:
    """Predição com Prodigal, se pyrodigal estiver instalado. Senão, None."""
    try:
        import pyrodigal
    except ImportError:
        return None

    finder = pyrodigal.GeneFinder(meta=True)
    out: list[ORF] = []
    for i, gene in enumerate(finder.find_genes(genome.seq.encode()), start=1):
        protein = gene.translate().rstrip("*")
        if min_aa <= len(protein) <= max_aa:
            out.append(
                ORF(
                    genome_id=genome.id,
                    index=i,
                    start=gene.begin,
                    end=gene.end,
                    strand=1 if gene.strand > 0 else -1,
                    protein=protein,
                )
            )
    return out


def predict_genes(
    genomes: Iterable[Record],
    min_aa: int = 30,
    max_aa: int = 400,
    prefer_prodigal: bool = True,
) -> tuple[list[ORF], str]:
    """Prediz genes em vários genomas. Devolve (ORFs, nome do método usado)."""
    orfs: list[ORF] = []
    method = "orf-interno"

    for genome in genomes:
        found = None
        if prefer_prodigal:
            found = find_orfs_pyrodigal(genome, min_aa, max_aa)
            if found is not None:
                method = "pyrodigal"
        if found is None:
            found = find_orfs(genome, min_aa, max_aa)
        orfs.extend(found)

    return orfs, method
