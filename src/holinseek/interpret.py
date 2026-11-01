"""Interpretação dos achados.

Duas camadas, nesta ordem:

1. **Determinística.** Para cada candidato, monta um texto de evidências a
   partir dos números medidos — por que ele passou, onde ficou fraco, o que
   precisa ser confirmado na bancada. Roda offline, é reprodutível e é o que
   entra no relatório oficial.

2. **Em prosa, opcional.** Se um provedor de linguagem estiver configurado no
   `.env`, ele recebe as evidências já calculadas e redige a leitura do
   conjunto. Nunca recebe a tarefa de decidir: recebe os números prontos e
   escreve sobre eles.

A separação é deliberada. Um modelo de linguagem é bom em redigir e péssimo
como fonte de verdade sobre uma medida que o próprio pipeline já computou.
"""

from __future__ import annotations

from dataclasses import dataclass

from .config import LLMSettings
from .llm import LLMUnavailable, complete, describe
from .scoring import Candidate

SYSTEM_PROMPT = (
    "Você é um bioinformata escrevendo a seção de resultados de uma triagem "
    "computacional de holinas em genomas de bacteriófago. Escreva em português "
    "do Brasil, de forma sóbria e técnica. "
    "Use SOMENTE os números fornecidos: não invente medidas, acessos, nomes de "
    "fago ou referências. Deixe claro que são predições que exigem confirmação "
    "experimental. Se a evidência for fraca, diga isso."
)


@dataclass
class Interpretation:
    """O que a camada de interpretação produziu para um candidato."""

    protein_id: str
    verdict: str
    evidence: list[str]
    caveats: list[str]
    next_steps: list[str]

    def to_text(self) -> str:
        parts = [f"## {self.protein_id}", "", self.verdict, "", "**Evidências**"]
        parts += [f"- {e}" for e in self.evidence]
        if self.caveats:
            parts += ["", "**Ressalvas**"] + [f"- {c}" for c in self.caveats]
        if self.next_steps:
            parts += ["", "**Próximo passo experimental**"] + [f"- {s}" for s in self.next_steps]
        return "\n".join(parts)


def _verdict_sentence(c: Candidate) -> str:
    f, cas = c.features, c.cassette
    if c.tier == "alta":
        return (
            f"Candidato forte. Proteína de {f.length} aminoácidos com "
            f"{f.n_tmd} hélice(s) transmembrana predita(s), hidrofobicidade "
            f"global de {f.gravy:+.2f} e cauda C-terminal de carga {f.tail_charge:+d}, "
            f"situada no mesmo módulo genômico de uma provável hidrolase de "
            f"parede celular. O conjunto reproduz a assinatura canônica de holina."
        )
    if c.tier == "média":
        return (
            f"Candidato plausível. Os atributos físico-químicos são compatíveis "
            f"({f.n_tmd} hélice(s), GRAVY {f.gravy:+.2f}, cauda {f.tail_charge:+d}), "
            f"mas {'falta a vizinhança lítica' if not cas.endolysin_nearby else 'parte dos critérios não foi atendida'}. "
            f"Vale manter na lista, abaixo dos candidatos fortes."
        )
    return (
        f"Candidato fraco. Atende {c.passed_count} de {c.total_checks} critérios. "
        f"Só faz sentido investigar se houver motivo externo — por exemplo, um "
        f"fago de interesse particular."
    )


def _evidence_lines(c: Candidate) -> list[str]:
    lines: list[str] = []
    for check in c.checks:
        mark = "atende" if check.passed else "não atende"
        lines.append(
            f"**{check.name.replace('_', ' ')}** — {mark}: "
            f"observado {check.observed}, esperado {check.expected}. {check.rationale}"
        )
    if c.cassette.endolysin_nearby and c.cassette.endolysin_reason:
        lines.append(
            f"**indício da endolisina** — {c.cassette.endolysin_id} foi classificada "
            f"como hidrolase por {c.cassette.endolysin_reason}."
        )
    return lines


def _caveat_lines(c: Candidate, tm_method: str) -> list[str]:
    out = [
        f"As hélices transmembrana vêm de {tm_method}. Confirme com um preditor "
        f"dedicado (DeepTMHMM, Phobius, TMbed) antes de qualquer conclusão.",
        "A predição não distingue holina de antiholina, nem determina a classe "
        "topológica — isso exige o mapeamento experimental da topologia.",
    ]
    if not c.cassette.endolysin_nearby:
        out.append(
            "Nenhuma endolisina foi detectada na vizinhança. Pode ser ausência real, "
            "mas também pode ser limite da detecção por motivo: uma busca HMMER "
            "contra Pfam/PHROG é mais sensível."
        )
    if c.features.n_tmd == 0:
        out.append(
            "Sem hélice predita, a hipótese de proteína de membrana fica sem base. "
            "Verifique se a ORF não está truncada."
        )
    if c.features.length < 60:
        out.append(
            "Proteínas muito curtas são as mais sujeitas a erro de predição de gene; "
            "confirme o códon de início."
        )
    return out


def _next_steps(c: Candidate) -> list[str]:
    if c.tier == "descartado":
        return []
    steps = [
        "Rodar DeepTMHMM ou Phobius para confirmar a topologia predita.",
        "Predizer a estrutura (ColabFold) e buscar parentes estruturais com Foldseek — "
        "holinas divergem em sequência, mas menos em forma.",
    ]
    if c.cassette.endolysin_nearby:
        steps.append(
            "Anotar o módulo inteiro em volta do candidato e verificar se o arranjo "
            "holina–endolisina–spanina está completo."
        )
    if c.tier == "alta":
        steps.append(
            "Sintetizar o gene e acompanhar o efeito sobre o cultivo por densidade "
            "óptica, sob promotor induzível, seguindo as normas de biossegurança da "
            "instituição."
        )
    return steps


def interpret_candidate(c: Candidate, tm_method: str = "hidropatia de Kyte-Doolittle") -> Interpretation:
    """Monta a interpretação determinística de um candidato."""
    return Interpretation(
        protein_id=c.protein_id,
        verdict=_verdict_sentence(c),
        evidence=_evidence_lines(c),
        caveats=_caveat_lines(c, tm_method),
        next_steps=_next_steps(c),
    )


def _summary_prompt(candidates: list[Candidate], interpretations: list[Interpretation]) -> str:
    counts = {"alta": 0, "média": 0, "baixa": 0, "descartado": 0}
    for c in candidates:
        counts[c.tier] = counts.get(c.tier, 0) + 1

    lines = [
        "Triagem computacional de holinas. Números medidos pelo pipeline:",
        "",
        f"- candidatos avaliados: {len(candidates)}",
        f"- prioridade alta: {counts['alta']}; média: {counts['média']}; baixa: {counts['baixa']}",
        "",
        "Os melhores candidatos, com seus atributos:",
        "",
    ]
    for c in candidates[:10]:
        f = c.features
        lines.append(
            f"- {c.protein_id} (genoma {c.genome_id}): nota {c.score:.2f}, faixa {c.tier}, "
            f"{f.length} aa, {f.n_tmd} TMD, GRAVY {f.gravy:+.2f}, cauda {f.tail_charge:+d}, "
            f"pI {f.isoelectric_point}, "
            f"endolisina vizinha: {'sim' if c.cassette.endolysin_nearby else 'não'}"
        )

    lines += [
        "",
        "Escreva de 3 a 5 parágrafos: o que esse conjunto sugere, quais candidatos "
        "merecem a bancada primeiro e por quê, e quais são os limites desta análise. "
        "Não invente nenhum dado além dos listados acima.",
    ]
    return "\n".join(lines)


def narrative_summary(
    candidates: list[Candidate],
    interpretations: list[Interpretation],
    settings: LLMSettings,
) -> tuple[str, str]:
    """Pede ao provedor configurado um resumo em prosa.

    Devolve (texto, situação). Se a camada estiver desligada ou falhar, o texto
    vem vazio e a situação explica o motivo — nunca levanta exceção, porque o
    relatório principal não depende disto.
    """
    if not settings.usable:
        return "", f"camada de prosa {describe(settings)}"

    try:
        text = complete(settings, SYSTEM_PROMPT, _summary_prompt(candidates, interpretations))
        return text, f"redigido por {describe(settings)}"
    except LLMUnavailable as exc:
        return "", f"camada de prosa indisponível ({exc}); relatório determinístico mantido"
