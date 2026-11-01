"""Camada opcional de redação em linguagem natural.

O pipeline nunca depende de um modelo de linguagem para decidir nada: a
pontuação, os vereditos e as evidências são todos calculados localmente, de
forma determinística. Esta camada apenas redige, em prosa, a leitura do que
já foi medido — e pode ser desligada sem perda de resultado.

Dois provedores vêm implementados:

  anthropic          API da Anthropic
  openai_compatible  qualquer endpoint /chat/completions — outro serviço
                     comercial ou um modelo aberto servido por você

Acrescentar um terceiro é escrever uma função `_call_<nome>` e registrá-la
em `_PROVIDERS`.
"""

from __future__ import annotations

from .config import LLMSettings


class LLMUnavailable(RuntimeError):
    """O provedor está desligado, mal configurado ou sem a biblioteca instalada."""


def _call_anthropic(settings: LLMSettings, system: str, prompt: str) -> str:
    try:
        import anthropic
    except ImportError as exc:  # pragma: no cover
        raise LLMUnavailable(
            "pacote 'anthropic' não instalado — rode: pip install anthropic"
        ) from exc

    kwargs = {}
    if settings.anthropic_api_key:
        kwargs["api_key"] = settings.anthropic_api_key
    client = anthropic.Anthropic(timeout=float(settings.timeout), **kwargs)

    try:
        response = client.messages.create(
            model=settings.anthropic_model,
            max_tokens=settings.max_tokens,
            system=system,
            output_config={"effort": settings.anthropic_effort},
            messages=[{"role": "user", "content": prompt}],
        )
    except anthropic.AuthenticationError as exc:
        raise LLMUnavailable("chave de API inválida ou ausente") from exc
    except anthropic.RateLimitError as exc:
        raise LLMUnavailable("limite de requisições atingido; tente mais tarde") from exc
    except anthropic.APIConnectionError as exc:
        raise LLMUnavailable(f"falha de conexão: {exc}") from exc
    except anthropic.APIStatusError as exc:
        raise LLMUnavailable(f"erro da API ({exc.status_code}): {exc.message}") from exc

    if response.stop_reason == "refusal":
        raise LLMUnavailable("a requisição foi recusada pelo provedor")

    return "\n".join(b.text for b in response.content if b.type == "text").strip()


def _call_openai_compatible(settings: LLMSettings, system: str, prompt: str) -> str:
    try:
        import requests
    except ImportError as exc:  # pragma: no cover
        raise LLMUnavailable("pacote 'requests' não instalado") from exc

    headers = {"Content-Type": "application/json"}
    if settings.api_key:
        headers["Authorization"] = f"Bearer {settings.api_key}"

    payload = {
        "model": settings.model,
        "max_tokens": settings.max_tokens,
        "messages": [
            {"role": "system", "content": system},
            {"role": "user", "content": prompt},
        ],
    }

    try:
        resp = requests.post(
            f"{settings.api_base}/chat/completions",
            json=payload,
            headers=headers,
            timeout=settings.timeout,
        )
    except requests.RequestException as exc:
        raise LLMUnavailable(f"falha de conexão com {settings.api_base}: {exc}") from exc

    if resp.status_code >= 400:
        raise LLMUnavailable(f"erro {resp.status_code} do endpoint: {resp.text[:200]}")

    try:
        return resp.json()["choices"][0]["message"]["content"].strip()
    except (KeyError, IndexError, ValueError) as exc:
        raise LLMUnavailable(f"resposta em formato inesperado: {resp.text[:200]}") from exc


_PROVIDERS = {
    "anthropic": _call_anthropic,
    "openai_compatible": _call_openai_compatible,
}


def complete(settings: LLMSettings, system: str, prompt: str) -> str:
    """Envia um pedido ao provedor configurado e devolve o texto da resposta."""
    if not settings.usable:
        raise LLMUnavailable(
            "nenhum provedor configurado — defina LLM_PROVIDER no .env "
            "(none | anthropic | openai_compatible)"
        )
    handler = _PROVIDERS.get(settings.provider)
    if handler is None:
        raise LLMUnavailable(f"provedor desconhecido: {settings.provider}")
    return handler(settings, system, prompt)


def describe(settings: LLMSettings) -> str:
    """Uma linha dizendo o que está configurado, para o cabeçalho do relatório."""
    if settings.provider == "anthropic":
        return f"anthropic / {settings.anthropic_model} (esforço {settings.anthropic_effort})"
    if settings.provider == "openai_compatible":
        return f"endpoint compatível / {settings.model} em {settings.api_base}"
    return "desligada — somente o relatório determinístico"
