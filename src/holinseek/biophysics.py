"""Atributos físico-químicos de uma proteína.

Este é o coração analítico do projeto. Holinas praticamente não conservam
sequência entre famílias, então a busca não pode ser por homologia: tem de ser
pelas propriedades que a função impõe.

Uma holina precisa atravessar o interior apolar da membrana interna — cerca de
4 nm de hidrocarboneto — e manter uma cauda citoplasmática carregada. Isso
deixa três marcas mensuráveis na sequência de aminoácidos:

  1. hidrofobicidade média alta (GRAVY);
  2. um a quatro trechos contínuos de alta hidropatia (hélices transmembrana);
  3. carga líquida positiva na cauda C-terminal (regra do interior positivo).

Todas as funções abaixo operam sobre a sequência de aminoácidos em letra única
e não dependem de nenhuma biblioteca externa.
"""

from __future__ import annotations

from dataclasses import dataclass, asdict
from typing import Sequence

# Escala de hidropatia de Kyte & Doolittle (J. Mol. Biol., 1982).
# Positivo = hidrofóbico; negativo = hidrofílico.
KYTE_DOOLITTLE: dict[str, float] = {
    "A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5,
    "Q": -3.5, "E": -3.5, "G": -0.4, "H": -3.2, "I": 4.5,
    "L": 3.8, "K": -3.9, "M": 1.9, "F": 2.8, "P": -1.6,
    "S": -0.8, "T": -0.7, "W": -0.9, "Y": -1.3, "V": 4.2,
}

# pKa dos grupos ionizáveis (valores do EMBOSS), usados para a carga líquida.
PKA_ACID = {"C-term": 3.1, "D": 3.9, "E": 4.1, "C": 8.5, "Y": 10.1}
PKA_BASE = {"N-term": 8.0, "H": 6.5, "K": 10.8, "R": 12.5}

POSITIVE = ("K", "R")
NEGATIVE = ("D", "E")

AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY"


def clean(seq: str) -> str:
    """Remove espaços, parada de tradução e resíduos fora do alfabeto padrão."""
    return "".join(c for c in seq.strip().upper() if c in AMINO_ACIDS)


def hydropathy_profile(seq: str, window: int = 19) -> list[float]:
    """Perfil de Kyte-Doolittle: média móvel de hidropatia.

    A janela de 19 resíduos é a convencional para hélices transmembrana —
    corresponde ao comprimento de hélice necessário para cruzar a bicamada.
    Devolve uma lista com um valor por janela (vazia se a sequência for curta).
    """
    seq = clean(seq)
    if window <= 0 or len(seq) < window:
        return []

    values = [KYTE_DOOLITTLE.get(aa, 0.0) for aa in seq]
    total = sum(values[:window])
    profile = [total / window]
    for i in range(window, len(values)):
        total += values[i] - values[i - window]
        profile.append(total / window)
    return profile


def gravy(seq: str) -> float:
    """GRAVY: hidropatia média de toda a proteína (grand average of hydropathy)."""
    seq = clean(seq)
    if not seq:
        return 0.0
    return sum(KYTE_DOOLITTLE.get(aa, 0.0) for aa in seq) / len(seq)


@dataclass(frozen=True)
class Segment:
    """Um trecho transmembrana predito, em coordenadas 1-based inclusivas."""

    start: int
    end: int
    peak: float

    @property
    def length(self) -> int:
        return self.end - self.start + 1


def predict_tm_segments(
    seq: str,
    window: int = 19,
    threshold: float = 1.6,
    min_length: int = 17,
) -> list[Segment]:
    """Prediz hélices transmembrana por hidropatia de janela deslizante.

    É o método clássico de Kyte-Doolittle: toda janela cuja hidropatia média
    ultrapassa o limiar é candidata a hélice; janelas vizinhas são fundidas num
    único segmento.

    É uma heurística, e assumidamente mais simples que um preditor treinado
    como DeepTMHMM ou Phobius. A vantagem é rodar offline, em qualquer máquina,
    sem download de modelo — o suficiente para a triagem de primeira passada.
    Para o resultado final, confirme com um preditor dedicado.
    """
    profile = hydropathy_profile(seq, window)
    if not profile:
        return []

    segments: list[Segment] = []
    start: int | None = None
    peak = float("-inf")

    for i, value in enumerate(profile):
        if value >= threshold:
            if start is None:
                start = i
                peak = value
            else:
                peak = max(peak, value)
        elif start is not None:
            segments.append(Segment(start + 1, i + window - 1, peak))
            start, peak = None, float("-inf")

    if start is not None:
        segments.append(Segment(start + 1, len(profile) + window - 2, peak))

    return [s for s in segments if s.length >= min_length]


def net_charge(seq: str, ph: float = 7.0) -> float:
    """Carga líquida da proteína no pH indicado (equação de Henderson-Hasselbalch)."""
    seq = clean(seq)
    if not seq:
        return 0.0

    charge = 0.0
    # Grupos básicos: protonados em pH baixo.
    charge += 1.0 / (1.0 + 10 ** (ph - PKA_BASE["N-term"]))
    for aa in ("H", "K", "R"):
        charge += seq.count(aa) / (1.0 + 10 ** (ph - PKA_BASE[aa]))
    # Grupos ácidos: desprotonados em pH alto.
    charge -= 1.0 / (1.0 + 10 ** (PKA_ACID["C-term"] - ph))
    for aa in ("D", "E", "C", "Y"):
        charge -= seq.count(aa) / (1.0 + 10 ** (PKA_ACID[aa] - ph))
    return charge


def simple_charge(seq: str) -> int:
    """Carga inteira contando apenas Lys/Arg menos Asp/Glu.

    É a forma usada na regra do interior positivo: q = (K + R) − (D + E).
    Mais grosseira que `net_charge`, porém é exatamente o critério descrito na
    literatura de topologia de membrana, então vale manter as duas.
    """
    seq = clean(seq)
    return sum(seq.count(a) for a in POSITIVE) - sum(seq.count(a) for a in NEGATIVE)


def c_terminal_tail(seq: str, segments: Sequence[Segment]) -> str:
    """Trecho após a última hélice transmembrana — a cauda citoplasmática."""
    seq = clean(seq)
    if not segments:
        return ""
    return seq[segments[-1].end :]


def isoelectric_point(seq: str, tolerance: float = 0.01) -> float:
    """Ponto isoelétrico, por busca binária no pH que zera a carga líquida."""
    low, high = 0.0, 14.0
    while high - low > tolerance:
        mid = (low + high) / 2
        if net_charge(seq, mid) > 0:
            low = mid
        else:
            high = mid
    return round((low + high) / 2, 2)


def aromaticity(seq: str) -> float:
    """Fração de resíduos aromáticos (Phe, Trp, Tyr)."""
    seq = clean(seq)
    if not seq:
        return 0.0
    return sum(seq.count(a) for a in "FWY") / len(seq)


@dataclass(frozen=True)
class Features:
    """Todos os atributos calculados para uma proteína."""

    length: int
    gravy: float
    n_tmd: int
    tm_fraction: float
    max_tm_peak: float
    tail_length: int
    tail_charge: int
    tail_charge_density: float
    net_charge_ph7: float
    isoelectric_point: float
    aromaticity: float
    has_start_met: bool

    def as_dict(self) -> dict:
        return asdict(self)


def compute_features(
    seq: str,
    window: int = 19,
    threshold: float = 1.6,
) -> Features:
    """Calcula o conjunto completo de atributos de uma sequência."""
    aa = clean(seq)
    segments = predict_tm_segments(aa, window=window, threshold=threshold)
    tail = c_terminal_tail(aa, segments)
    tm_residues = sum(s.length for s in segments)

    return Features(
        length=len(aa),
        gravy=round(gravy(aa), 4),
        n_tmd=len(segments),
        tm_fraction=round(tm_residues / len(aa), 4) if aa else 0.0,
        max_tm_peak=round(max((s.peak for s in segments), default=0.0), 4),
        tail_length=len(tail),
        tail_charge=simple_charge(tail),
        tail_charge_density=round(simple_charge(tail) / len(tail), 4) if tail else 0.0,
        net_charge_ph7=round(net_charge(aa), 3),
        isoelectric_point=isoelectric_point(aa),
        aromaticity=round(aromaticity(aa), 4),
        has_start_met=aa.startswith("M"),
    )
