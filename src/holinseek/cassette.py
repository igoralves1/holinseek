"""Contexto genômico: o cassete de lise.

O filtro isolado mais forte não é nenhuma propriedade da proteína sozinha — é
a vizinhança. Em bacteriófagos, holina, endolisina e spanina ficam quase
sempre lado a lado, num módulo contíguo do genoma. Uma proteína pequena e
hidrofóbica ao lado de uma hidrolase de parede celular é um candidato muito
mais forte do que a mesma proteína isolada em outra região.

Este módulo procura endolisinas na vizinhança de cada candidato. A detecção é
feita por duas vias:

  1. assinaturas de sequência de hidrolases de peptidoglicano (motivos curtos
     e conservados de amidases, muramidases e transglicosilases líticas);
  2. anotação textual, quando o FASTA de entrada já traz descrições.

Para uma identificação mais confiável, rode uma busca HMMER contra Pfam/PHROG
e passe o resultado em `annotations` — ver `holinseek.domains`.
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from typing import Iterable, Sequence

from .orfs import ORF

# Palavras que, numa anotação, indicam uma enzima lítica de parede celular.
ENDOLYSIN_WORDS = (
    "endolysin", "endolisina", "lysozyme", "lisozima", "lysin",
    "amidase", "muramidase", "peptidoglycan hydrolase", "hidrolase",
    "transglycosylase", "transglicosilase", "glucosaminidase",
    "n-acetylmuramoyl", "cell wall hydrolase", "autolysin",
)

SPANIN_WORDS = ("spanin", "spanina", "rz", "rz1")
HOLIN_WORDS = ("holin", "holina", "antiholin", "antiholina")

# Motivos curtos de famílias de hidrolase de peptidoglicano. São indícios,
# não identificações: um acerto aumenta a confiança, a ausência não exclui.
ENDOLYSIN_MOTIFS: tuple[tuple[str, str], ...] = (
    (r"C[DE].{1,3}H.{10,30}[DE].{10,40}H", "sítio de amidase/peptidase dependente de zinco"),
    (r"E.{10,25}[DN].{3,12}T", "arranjo catalítico tipo muramidase"),
    (r"G.{2}[ST].{2}G.{3}[DE]", "núcleo tipo transglicosilase lítica"),
)


@dataclass(frozen=True)
class CassetteEvidence:
    """O que a vizinhança genômica diz sobre um candidato."""

    neighbour_count: int
    endolysin_nearby: bool
    endolysin_id: str | None
    endolysin_distance: int | None
    endolysin_reason: str | None
    spanin_nearby: bool
    overlapping_orf: bool

    @property
    def score(self) -> float:
        """Contribuição da sintenia para a pontuação final, de 0 a 1."""
        value = 0.0
        if self.endolysin_nearby:
            value += 0.7
            if self.endolysin_distance is not None and self.endolysin_distance <= 500:
                value += 0.1
        if self.spanin_nearby:
            value += 0.1
        if self.overlapping_orf:
            # ORFs embutidas em fase deslocada são a assinatura clássica do
            # par holina/antiholina.
            value += 0.1
        return min(value, 1.0)


def looks_like_endolysin(protein: str, description: str = "") -> tuple[bool, str | None]:
    """Decide se uma proteína parece uma hidrolase de parede celular."""
    text = description.lower()
    for word in ENDOLYSIN_WORDS:
        if word in text:
            return True, f"anotação contém «{word}»"

    # Endolisinas são tipicamente maiores que holinas.
    if not 100 <= len(protein) <= 600:
        return False, None

    for pattern, label in ENDOLYSIN_MOTIFS:
        if re.search(pattern, protein):
            return True, label
    return False, None


def _mentions(description: str, words: Sequence[str]) -> bool:
    text = description.lower()
    return any(w in text for w in words)


def analyse_neighbourhood(
    candidate: ORF,
    all_orfs: Iterable[ORF],
    window_bp: int = 6000,
    annotations: dict[str, str] | None = None,
) -> CassetteEvidence:
    """Examina as ORFs vizinhas de um candidato dentro de uma janela em pares de base."""
    annotations = annotations or {}
    neighbours = [
        o
        for o in all_orfs
        if o.genome_id == candidate.genome_id
        and o.protein_id != candidate.protein_id
        and abs(o.midpoint - candidate.midpoint) <= window_bp
    ]

    best_id: str | None = None
    best_distance: int | None = None
    best_reason: str | None = None
    spanin = False
    overlapping = False

    for n in neighbours:
        desc = annotations.get(n.protein_id, "")

        # Sobreposição física com outra ORF, em qualquer fita.
        if n.start <= candidate.end and candidate.start <= n.end:
            overlapping = True

        if _mentions(desc, SPANIN_WORDS):
            spanin = True

        is_endo, reason = looks_like_endolysin(n.protein, desc)
        if is_endo:
            distance = int(abs(n.midpoint - candidate.midpoint))
            if best_distance is None or distance < best_distance:
                best_id, best_distance, best_reason = n.protein_id, distance, reason

    return CassetteEvidence(
        neighbour_count=len(neighbours),
        endolysin_nearby=best_id is not None,
        endolysin_id=best_id,
        endolysin_distance=best_distance,
        endolysin_reason=best_reason,
        spanin_nearby=spanin,
        overlapping_orf=overlapping,
    )
