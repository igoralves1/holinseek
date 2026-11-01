"""HolinSeek — descoberta computacional de holinas em genomas de bacteriófago.

Holinas são proteínas pequenas e hidrofóbicas que formam poros na membrana
interna da bactéria e determinam o momento da lise mediada por fago. Entre
famílias elas quase não conservam sequência, de modo que uma busca por
homologia encontra pouco. Este projeto busca, em vez disso, as propriedades
que a função impõe: hidrofobicidade, topologia de membrana, carga da cauda
citoplasmática e vizinhança genômica com uma endolisina.

Uso rápido:

    from holinseek.config import load_settings
    from holinseek.pipeline import run

    settings = load_settings()
    resultado = run("data/raw/genomas.fasta", settings)
    print(len(resultado.candidates), "candidatos")

Ou pela linha de comando:

    holinseek run data/raw/genomas.fasta
"""

__version__ = "0.1.0"

__all__ = ["__version__"]
