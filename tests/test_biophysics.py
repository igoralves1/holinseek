"""Testes dos atributos físico-químicos."""

from holinseek.biophysics import (
    KYTE_DOOLITTLE,
    clean,
    compute_features,
    gravy,
    hydropathy_profile,
    isoelectric_point,
    net_charge,
    predict_tm_segments,
    simple_charge,
)


def test_escala_kyte_doolittle_completa():
    assert len(KYTE_DOOLITTLE) == 20
    assert KYTE_DOOLITTLE["I"] == 4.5   # o mais hidrofóbico
    assert KYTE_DOOLITTLE["R"] == -4.5  # o mais hidrofílico


def test_clean_remove_residuos_invalidos():
    assert clean(" mkv*x z\n") == "MKVX".replace("X", "")  # X não é padrão


def test_gravy_de_sequencia_homogenea():
    # Uma cadeia só de isoleucina tem GRAVY igual ao valor da isoleucina.
    assert gravy("I" * 20) == 4.5
    assert gravy("R" * 20) == -4.5
    assert gravy("") == 0.0


def test_gravy_hidrofobico_vs_hidrofilico():
    assert gravy("LIVFAMLIVFAM") > 0
    assert gravy("DEKRNQDEKRNQ") < 0


def test_perfil_tem_tamanho_esperado():
    seq = "A" * 30
    perfil = hydropathy_profile(seq, window=19)
    assert len(perfil) == 30 - 19 + 1
    assert hydropathy_profile("AAA", window=19) == []


def test_helice_transmembrana_detectada():
    # Um trecho hidrofóbico longo cercado de resíduos carregados.
    seq = "KKKKK" + "LIAVGLIAGLFIAVGMLVAL" * 2 + "KKKKK"
    segmentos = predict_tm_segments(seq, window=19, threshold=1.6)
    assert len(segmentos) >= 1
    assert all(s.length >= 17 for s in segmentos)
    assert all(s.peak >= 1.6 for s in segmentos)


def test_sem_helice_em_proteina_soluvel():
    seq = "MKDEKRNQSTDEKRNQSTDEKRNQSTDEKRNQSTDEKRNQST"
    assert predict_tm_segments(seq) == []


def test_carga_simples_conta_apenas_kr_de():
    assert simple_charge("KKKRRR") == 6
    assert simple_charge("DDDEEE") == -6
    assert simple_charge("KKKDDD") == 0
    assert simple_charge("AAAGGG") == 0


def test_carga_liquida_segue_o_ph():
    seq = "KKKRRRDDDEEE"
    assert net_charge(seq, ph=2.0) > net_charge(seq, ph=12.0)


def test_ponto_isoeletrico_em_faixa_plausivel():
    assert 9.0 < isoelectric_point("KKKKRRRR") <= 14.0
    assert 0.0 <= isoelectric_point("DDDDEEEE") < 5.0


def test_compute_features_de_perfil_de_holina():
    seq = "M" + "SNNEQ" + "LIAVGLIAGLFIAVGMLVAL" + "GSNTSG" + "VGAILFALIVGAVLVGFIAL" + "GKKRKRRKQAE"
    f = compute_features(seq)
    assert 50 <= f.length <= 150
    assert f.n_tmd >= 1
    assert f.gravy > 0.4
    assert f.tail_charge >= 3
    assert f.has_start_met is True
