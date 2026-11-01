"""Testes dos sete critérios e da pontuação."""

from holinseek.biophysics import compute_features
from holinseek.cassette import CassetteEvidence
from holinseek.config import Criteria
from holinseek.scoring import evaluate, rank

HOLINA = (
    "M" + "SNNEQ" + "LIAVGLIAGLFIAVGMLVAL" + "GSNTSG"
    + "VGAILFALIVGAVLVGFIAL" + "GKKRKRRKQAE"
)
SOLUVEL = "M" + "DEKRNQST" * 10

SEM_VIZINHANCA = CassetteEvidence(
    neighbour_count=3,
    endolysin_nearby=False,
    endolysin_id=None,
    endolysin_distance=None,
    endolysin_reason=None,
    spanin_nearby=False,
    overlapping_orf=False,
)
COM_VIZINHANCA = CassetteEvidence(
    neighbour_count=3,
    endolysin_nearby=True,
    endolysin_id="g_orf0002",
    endolysin_distance=320,
    endolysin_reason="sítio de amidase",
    spanin_nearby=True,
    overlapping_orf=True,
)


def _avaliar(seq, cassete):
    return evaluate("p1", "g1", seq, compute_features(seq), cassete, Criteria())


def test_sete_criterios_sao_avaliados():
    c = _avaliar(HOLINA, COM_VIZINHANCA)
    assert c.total_checks == 7
    nomes = {chk.name for chk in c.checks}
    assert "comprimento" in nomes
    assert "domínios_transmembrana" in nomes
    assert "endolisina_vizinha" in nomes


def test_holina_com_vizinhanca_tem_prioridade_alta():
    c = _avaliar(HOLINA, COM_VIZINHANCA)
    assert c.score >= 0.85
    assert c.tier == "alta"


def test_vizinhanca_muda_a_faixa():
    com = _avaliar(HOLINA, COM_VIZINHANCA)
    sem = _avaliar(HOLINA, SEM_VIZINHANCA)
    assert com.score > sem.score
    assert com.tier == "alta"
    assert sem.tier != "alta", "sem endolisina vizinha não deveria ser prioridade alta"


def test_proteina_soluvel_e_descartada():
    c = _avaliar(SOLUVEL, SEM_VIZINHANCA)
    assert c.features.n_tmd == 0
    assert c.score < 0.5
    assert c.tier == "descartado"


def test_cada_criterio_explica_se():
    c = _avaliar(HOLINA, COM_VIZINHANCA)
    for chk in c.checks:
        assert chk.rationale.strip(), f"critério {chk.name} sem justificativa"
        assert chk.expected.strip()


def test_rank_ordena_e_corta():
    bons = _avaliar(HOLINA, COM_VIZINHANCA)
    ruins = _avaliar(SOLUVEL, SEM_VIZINHANCA)
    ordenado = rank([ruins, bons], min_score=0.5)
    assert ordenado == [bons]


def test_peso_da_sintenia():
    assert COM_VIZINHANCA.score > SEM_VIZINHANCA.score
    assert 0.0 <= COM_VIZINHANCA.score <= 1.0
