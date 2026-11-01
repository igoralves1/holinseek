"""Testes de ponta a ponta sobre o genoma sintético de exemplo."""

from pathlib import Path

import pytest

from holinseek.cassette import looks_like_endolysin
from holinseek.config import Criteria, Settings
from holinseek.orfs import find_orfs, reverse_complement, translate
from holinseek.pipeline import run
from holinseek.seqio import read_fasta

EXEMPLO = Path(__file__).resolve().parents[1] / "data" / "exemplo" / "fago_sintetico.fasta"


def test_exemplo_existe():
    assert EXEMPLO.exists(), "o genoma sintético de exemplo deveria estar versionado"


def test_traducao_para_no_codon_de_parada():
    assert translate("ATGAAATAA") == "MK"
    assert translate("ATGATGATG") == "MMM"


def test_complemento_reverso():
    assert reverse_complement("ATGC") == "GCAT"
    assert reverse_complement(reverse_complement("ATGCGTA")) == "ATGCGTA"


def test_orfs_encontradas_no_exemplo():
    genomas = read_fasta(EXEMPLO)
    assert len(genomas) == 1
    orfs = find_orfs(genomas[0], min_aa=30, max_aa=400)
    assert len(orfs) >= 2
    assert all(o.protein.startswith("M") for o in orfs)


def test_detector_de_endolisina():
    # Contém o arranjo CxxH ... D ... H que o detector procura.
    positivo = "M" + "A" * 28 + "CDAAH" + "G" * 18 + "D" + "S" * 24 + "H" + "A" * 60
    achou, motivo = looks_like_endolysin(positivo)
    assert achou and motivo

    # Curta demais para ser uma hidrolase.
    assert looks_like_endolysin("MKVLIAG")[0] is False


def test_pipeline_completo(tmp_path):
    settings = Settings(
        data_dir=tmp_path,
        results_dir=tmp_path / "results",
        criteria=Criteria(),
    )
    resultado = run(EXEMPLO, settings, min_score=0.4)

    assert resultado.genomes == 1
    assert resultado.orfs >= 2
    assert len(resultado.candidates) >= 1

    melhor = resultado.candidates[0]
    assert melhor.score >= 0.8
    assert melhor.tier in {"alta", "média"}
    assert melhor.features.n_tmd >= 1
    assert melhor.cassette.endolysin_nearby

    for nome in ("tabela", "json", "relatorio"):
        assert resultado.outputs[nome].exists()

    texto = resultado.outputs["relatorio"].read_text(encoding="utf-8")
    assert "predição" in texto.lower()


def test_pipeline_recusa_entrada_vazia(tmp_path):
    vazio = tmp_path / "vazio.fasta"
    vazio.write_text("", encoding="utf-8")
    settings = Settings(data_dir=tmp_path, results_dir=tmp_path / "r")
    with pytest.raises(ValueError):
        run(vazio, settings)
