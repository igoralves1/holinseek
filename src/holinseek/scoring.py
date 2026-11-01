"""Os sete filtros, aplicados e pontuados.

Cada critério do estudo vira aqui uma verificação isolada, com o valor medido,
o valor esperado e um veredito. Guardar critério a critério — em vez de só a
nota final — é o que torna o resultado auditável: para qualquer candidato dá
para dizer exatamente por que ele subiu ou desceu na lista.
"""

from __future__ import annotations

from dataclasses import dataclass, field, asdict

from .biophysics import Features
from .cassette import CassetteEvidence
from .config import Criteria


@dataclass(frozen=True)
class Check:
    """O resultado de um único critério."""

    name: str
    passed: bool
    observed: float | int | str
    expected: str
    weight: float
    rationale: str

    def as_dict(self) -> dict:
        return asdict(self)


@dataclass
class Candidate:
    """Uma proteína avaliada contra todos os critérios."""

    protein_id: str
    genome_id: str
    sequence: str
    features: Features
    cassette: CassetteEvidence
    checks: list[Check] = field(default_factory=list)

    @property
    def passed_count(self) -> int:
        return sum(1 for c in self.checks if c.passed)

    @property
    def total_checks(self) -> int:
        return len(self.checks)

    @property
    def score(self) -> float:
        """Nota de 0 a 1: soma ponderada dos critérios atendidos."""
        total_weight = sum(c.weight for c in self.checks) or 1.0
        earned = sum(c.weight for c in self.checks if c.passed)
        return round(earned / total_weight, 4)

    @property
    def tier(self) -> str:
        """Faixa de prioridade para a bancada."""
        s = self.score
        if s >= 0.85 and self.cassette.endolysin_nearby:
            return "alta"
        if s >= 0.7:
            return "média"
        if s >= 0.5:
            return "baixa"
        return "descartado"

    def as_row(self) -> dict:
        """Uma linha achatada, pronta para ir a uma tabela."""
        row = {
            "protein_id": self.protein_id,
            "genome_id": self.genome_id,
            "score": self.score,
            "tier": self.tier,
            "checks_passed": f"{self.passed_count}/{self.total_checks}",
        }
        row.update(self.features.as_dict())
        row.update(
            {
                "endolysin_nearby": self.cassette.endolysin_nearby,
                "endolysin_id": self.cassette.endolysin_id or "",
                "endolysin_distance_bp": self.cassette.endolysin_distance or "",
                "spanin_nearby": self.cassette.spanin_nearby,
                "overlapping_orf": self.cassette.overlapping_orf,
                "neighbours_examined": self.cassette.neighbour_count,
            }
        )
        row["sequence"] = self.sequence
        return row


def evaluate(
    protein_id: str,
    genome_id: str,
    sequence: str,
    features: Features,
    cassette: CassetteEvidence,
    criteria: Criteria,
) -> Candidate:
    """Aplica os sete critérios e devolve o candidato pontuado."""
    checks: list[Check] = [
        Check(
            name="comprimento",
            passed=criteria.min_length <= features.length <= criteria.max_length,
            observed=features.length,
            expected=f"{criteria.min_length}–{criteria.max_length} aminoácidos",
            weight=1.0,
            rationale=(
                "Holinas estão entre as menores proteínas funcionais de um genoma "
                "viral. Um polipeptídeo muito maior raramente tem a arquitetura "
                "compacta de poro que o mecanismo exige."
            ),
        ),
        Check(
            name="domínios_transmembrana",
            passed=criteria.min_tmd <= features.n_tmd <= criteria.max_tmd,
            observed=features.n_tmd,
            expected=f"{criteria.min_tmd}–{criteria.max_tmd} hélices",
            weight=2.0,
            rationale=(
                "O número de hélices separa as classes I, II e III de holina. "
                "Nenhuma hélice significa que a proteína não se insere na "
                "membrana; muitas indicam um transportador, não uma holina."
            ),
        ),
        Check(
            name="hidrofobicidade",
            passed=features.gravy > criteria.min_gravy,
            observed=features.gravy,
            expected=f"GRAVY > {criteria.min_gravy}",
            weight=1.5,
            rationale=(
                "A proteína precisa particionar no interior apolar da bicamada. "
                "GRAVY é a média de Kyte-Doolittle sobre toda a cadeia e mede "
                "exatamente essa afinidade."
            ),
        ),
        Check(
            name="carga_da_cauda",
            passed=features.tail_charge >= criteria.min_tail_charge,
            observed=features.tail_charge,
            expected=f"q ≥ +{criteria.min_tail_charge}",
            weight=1.5,
            rationale=(
                "Regra do interior positivo: a cauda que fica no citoplasma "
                "concentra Lys e Arg. É um dos critérios mais discriminantes "
                "entre holinas e outras proteínas pequenas de membrana."
            ),
        ),
        Check(
            name="início_canônico",
            passed=features.has_start_met,
            observed="M" if features.has_start_met else "ausente",
            expected="metionina inicial",
            weight=0.5,
            rationale=(
                "Checagem de sanidade da predição de genes. A ausência costuma "
                "indicar uma ORF truncada, não uma proteína real."
            ),
        ),
        Check(
            name="endolisina_vizinha",
            passed=cassette.endolysin_nearby,
            observed=(
                f"{cassette.endolysin_id} a {cassette.endolysin_distance} pb"
                if cassette.endolysin_nearby
                else "nenhuma"
            ),
            expected="hidrolase de parede celular no mesmo módulo",
            weight=2.5,
            rationale=(
                "O critério de maior peso. Holina e endolisina formam um módulo "
                "funcional e viajam juntas no genoma: uma sem a outra raramente "
                "faz sentido biológico."
            ),
        ),
        Check(
            name="fase_sobreposta",
            passed=cassette.overlapping_orf,
            observed="sim" if cassette.overlapping_orf else "não",
            expected="ORF embutida em fase deslocada (indício de antiholina)",
            weight=0.5,
            rationale=(
                "Muitos cassetes codificam a antiholina numa fase deslocada "
                "dentro do mesmo trecho. Quando aparece, reforça bastante a "
                "hipótese; quando falta, não a enfraquece."
            ),
        ),
    ]

    return Candidate(
        protein_id=protein_id,
        genome_id=genome_id,
        sequence=sequence,
        features=features,
        cassette=cassette,
        checks=checks,
    )


def rank(candidates: list[Candidate], min_score: float = 0.5) -> list[Candidate]:
    """Ordena por nota decrescente e descarta o que está abaixo do corte."""
    kept = [c for c in candidates if c.score >= min_score]
    kept.sort(key=lambda c: (c.score, c.cassette.score, c.features.gravy), reverse=True)
    return kept
