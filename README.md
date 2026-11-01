# HolinSeek

**Descoberta computacional de holinas em genomas de bacteriófago.**

Site do projeto: **https://igoralves1.github.io/holinseek/**

Holinas são proteínas pequenas e hidrofóbicas que formam poros na membrana
interna da bactéria e determinam o *momento* da lise mediada por fago. Elas
interessam porque abrem um mecanismo antimicrobiano diferente dos alvos
clássicos — e porque são difíceis de achar.

O problema central: **entre famílias, holinas quase não conservam sequência.**
Uma busca por homologia encontra pouco. Este projeto busca, em vez disso, as
propriedades que a função impõe.

```
              ┌─ hidrofobicidade (GRAVY)
              │
sequência ────┼─ topologia de membrana (hélices transmembrana)
de aminoácidos│
              ├─ carga da cauda citoplasmática (regra do interior positivo)
              │
              └─ vizinhança genômica (o cassete de lise)
```

---

## Instalação

Requer Python 3.10 ou superior.

```bash
git clone https://github.com/igoralves1/holinseek.git
cd holinseek

python -m venv .venv
source .venv/bin/activate          # no Windows: .venv\Scripts\activate

pip install -r requirements.txt
pip install -e .

cp .env.example .env               # e preencha NCBI_EMAIL
```

Confira o que ficou instalado e configurado:

```bash
holinseek check
```

### Dependências opcionais

Nada abaixo é obrigatório — cada item apenas troca uma implementação interna
por uma ferramenta mais precisa:

```bash
pip install -r requirements-optional.txt
```

| Pacote | O que melhora |
|---|---|
| `pyrodigal` | Predição de genes com Prodigal em vez do localizador de ORFs interno |
| `pyhmmer` | Busca HMMER contra Pfam/PHROG sem instalar o binário |
| `torch`, `transformers` | Embeddings de modelos de linguagem de proteína |
| `umap-learn` | Projeção dos candidatos para as figuras |

---

## Uso

### Demonstração, sem rede

```bash
holinseek demo
```

Roda o pipeline inteiro sobre um genoma **sintético** versionado em
`data/exemplo/`. Serve para verificar a instalação. Não é um organismo real e
não produz resultado biológico.

### Pipeline completo

```bash
holinseek fetch --query fagos-refseq --limit 50     # baixa genomas do NCBI
holinseek run data/raw/genomas.fasta                # triagem
```

Saídas em `data/results/`:

| Arquivo | Conteúdo |
|---|---|
| `candidatos.tsv` | Uma linha por candidato, com todos os atributos |
| `candidatos.json` | O mesmo, mais o veredito critério a critério |
| `relatorio.md` | Relatório legível, com as evidências e as ressalvas |

### Agrupamento e classificação

```bash
holinseek cluster data/raw/genomas.fasta
holinseek train   data/raw/genomas.fasta --rotulos meus_rotulos.csv
holinseek predict data/raw/genomas.fasta --modelo models/holin_classifier.joblib
```

---

## Como funciona

### Os sete critérios

| Critério | Faixa esperada | Peso | Por quê |
|---|---|---|---|
| Comprimento | 50–150 aa | 1,0 | Holinas estão entre as menores proteínas funcionais de um genoma viral |
| Hélices transmembrana | 1–4 | 2,0 | Define as classes I, II e III; zero significa que não se insere na membrana |
| Hidrofobicidade | GRAVY > 0,4 | 1,5 | A proteína precisa particionar no interior apolar da bicamada |
| Carga da cauda | q ≥ +3 | 1,5 | Regra do interior positivo: a cauda citoplasmática concentra Lys e Arg |
| Metionina inicial | presente | 0,5 | Checagem de sanidade da predição de genes |
| Endolisina vizinha | mesmo módulo | **2,5** | O critério de maior peso — holina e endolisina viajam juntas |
| Fase sobreposta | ORF embutida | 0,5 | Assinatura clássica do par holina/antiholina |

Todos os limiares são configuráveis pelo `.env`, sem tocar no código.

### A física por trás

A membrana interna de *E. coli* é uma bicamada de fosfatidiletanolamina
(C₃₉H₇₆NO₈P), fosfatidilglicerol (C₄₀H₇₇O₁₀P) e cardiolipina (C₈₁H₁₅₆O₁₇P₂),
com cerca de 4 nm de interior apolar. Atravessá-la exige uma hélice de
aproximadamente 19 resíduos hidrofóbicos — e é exatamente o comprimento da
janela usada no perfil de Kyte-Doolittle.

A holina não cliva nada. O gatilho é elétrico: ao atingir a concentração
crítica na membrana, o poro se abre e colapsa o potencial

```
Δp = Δψ − (2,303·RT/F)·ΔpH
```

Só então a endolisina alcança o peptidoglicano e hidrolisa a ligação β-1,4
entre GlcNAc (C₈H₁₅NO₆) e MurNAc (C₁₁H₁₉NO₈).

### Limitações, declaradas

- O preditor de hélices transmembrana interno é a heurística de
  Kyte-Doolittle, **não** um preditor treinado. É suficiente para a primeira
  peneira; confirme com DeepTMHMM, Phobius ou TMbed.
- A detecção de endolisina por motivo é um indício, não uma identificação.
  Uma busca HMMER contra Pfam/PHROG é bem mais sensível.
- O pipeline não distingue holina de antiholina nem determina a classe
  topológica.
- **Nenhum candidato é uma holina até que a bancada mostre que é.**

---

## A camada de interpretação

O relatório de evidências é montado **localmente**, a partir dos números que o
pipeline calculou. Para cada candidato ele diz o que foi medido, qual critério
passou, onde está fraco e qual seria o próximo passo experimental. Isso não
depende de serviço externo nenhum e é reprodutível.

Opcionalmente, um modelo de linguagem pode redigir a leitura do conjunto em
prosa. Configure em `.env`:

```ini
LLM_PROVIDER=none              # none | anthropic | openai_compatible
```

O provedor `openai_compatible` aceita qualquer endpoint `/chat/completions`,
o que inclui modelos abertos servidos localmente por Ollama, vLLM ou
llama.cpp. Em nenhum caso o modelo decide alguma coisa: ele recebe os números
já calculados e escreve sobre eles.

---

## Aprendizado de máquina

**Não supervisionado** (não precisa de rótulo): agrupa os candidatos no espaço
de atributos e aponta os isolados.

**Supervisionado** (precisa dos seus rótulos): treina um classificador a partir
de um CSV `protein_id,label`. O projeto **não embarca rótulos prontos** de
propósito — um conjunto de treino herdado sem auditoria é o caminho mais curto
para um classificador que aprendeu a base de dados em vez da biologia.

A validação cruzada é **agrupada por genoma** (`GroupKFold`). Sem isso, duas
proteínas do mesmo fago caem uma no treino e outra no teste, e a acurácia
medida vira ficção. É a crítica metodológica que mais aparece em banca.

---

## Desenvolvimento

```bash
pytest                    # 25 testes
make demo                 # pipeline sobre o exemplo sintético
make site                 # serve o site local em http://localhost:8000
```

Estrutura:

```
src/holinseek/
  config.py       configuração e os sete critérios, lidos do .env
  seqio.py        leitura de FASTA
  ncbi.py         download de genomas das E-utilities
  orfs.py         predição de genes (Prodigal, ou localizador interno)
  biophysics.py   GRAVY, hélices transmembrana, carga, pI
  cassette.py     vizinhança genômica e detecção de endolisina
  scoring.py      os sete critérios aplicados e pontuados
  ml.py           agrupamento e classificador
  llm.py          provedores de linguagem (opcional)
  interpret.py    o motor de interpretação
  report.py       tabela, JSON e relatório
  pipeline.py     orquestração
  cli.py          linha de comando
docs/             site publicado em GitHub Pages
```

---

## Biossegurança

Trabalho com bacteriófagos e proteínas líticas segue as normas de
biossegurança da instituição e, no caso terapêutico, a regulação sanitária
aplicável. A triagem computacional **antecede** essa avaliação; não a
substitui nem a dispensa.

---

## Autoria

- **Deraldo Duran** — Universidade Federal do Ceará (UFC)
- **Igor Alves** — Especialista em Bioinformática, UQAM (Canadá)

## Licença

MIT. Ver [LICENSE](LICENSE).
