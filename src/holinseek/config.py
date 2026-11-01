"""Configuração do pipeline, lida do ambiente e do arquivo .env.

Um único objeto `Settings` carrega tudo que o pipeline precisa saber. Os
valores vêm, nesta ordem de precedência: argumento explícito no código,
variável de ambiente, arquivo .env, padrão embutido.
"""

from __future__ import annotations

import os
from dataclasses import dataclass, field
from pathlib import Path

try:  # python-dotenv é opcional: sem ele, só o ambiente é lido.
    from dotenv import load_dotenv
except ImportError:  # pragma: no cover
    def load_dotenv(*_args, **_kwargs) -> bool:
        return False


def _env_str(name: str, default: str = "") -> str:
    return os.environ.get(name, default).strip()


def _env_int(name: str, default: int) -> int:
    raw = os.environ.get(name, "").strip()
    try:
        return int(raw) if raw else default
    except ValueError:
        return default


def _env_float(name: str, default: float) -> float:
    raw = os.environ.get(name, "").strip().replace(",", ".")
    try:
        return float(raw) if raw else default
    except ValueError:
        return default


def _env_bool(name: str, default: bool = False) -> bool:
    raw = os.environ.get(name, "").strip().lower()
    if not raw:
        return default
    return raw in {"1", "true", "yes", "y", "sim", "on"}


@dataclass(frozen=True)
class Criteria:
    """Os sete filtros de triagem descritos no estudo.

    Cada um traduz uma restrição física ou genômica em um número que pode ser
    calculado a partir da sequência de aminoácidos.
    """

    min_length: int = 50
    max_length: int = 150
    min_tmd: int = 1
    max_tmd: int = 4
    min_gravy: float = 0.4
    min_tail_charge: int = 3
    tm_window: int = 19
    tm_threshold: float = 1.6

    @classmethod
    def from_env(cls) -> "Criteria":
        return cls(
            min_length=_env_int("HOLINSEEK_MIN_LENGTH", 50),
            max_length=_env_int("HOLINSEEK_MAX_LENGTH", 150),
            min_tmd=_env_int("HOLINSEEK_MIN_TMD", 1),
            max_tmd=_env_int("HOLINSEEK_MAX_TMD", 4),
            min_gravy=_env_float("HOLINSEEK_MIN_GRAVY", 0.4),
            min_tail_charge=_env_int("HOLINSEEK_MIN_TAIL_CHARGE", 3),
            tm_window=_env_int("HOLINSEEK_TM_WINDOW", 19),
            tm_threshold=_env_float("HOLINSEEK_TM_THRESHOLD", 1.6),
        )


@dataclass(frozen=True)
class LLMSettings:
    """Camada opcional de redação em linguagem natural.

    O pipeline não depende disto: o relatório de evidências é sempre montado
    localmente a partir dos números calculados. Isto apenas acrescenta a
    leitura em prosa de cada candidato.

    Dois provedores são suportados — a API da Anthropic e qualquer endpoint
    compatível com /chat/completions, o que inclui modelos abertos servidos
    localmente.
    """

    provider: str = "none"          # none | anthropic | openai_compatible
    anthropic_api_key: str = ""
    anthropic_model: str = "claude-opus-5-5"
    anthropic_effort: str = "medium"
    api_base: str = ""
    api_key: str = ""
    model: str = ""
    timeout: int = 120
    max_tokens: int = 4000

    @property
    def usable(self) -> bool:
        if self.provider == "anthropic":
            # A chave também pode vir do ambiente lido pelo próprio SDK.
            return bool(self.anthropic_model)
        if self.provider == "openai_compatible":
            return bool(self.api_base and self.model)
        return False

    @classmethod
    def from_env(cls) -> "LLMSettings":
        provider = _env_str("LLM_PROVIDER", "none").lower() or "none"
        if provider not in {"none", "anthropic", "openai_compatible"}:
            provider = "none"
        return cls(
            provider=provider,
            anthropic_api_key=_env_str("ANTHROPIC_API_KEY"),
            anthropic_model=_env_str("ANTHROPIC_MODEL", "claude-opus-5-5"),
            anthropic_effort=_env_str("ANTHROPIC_EFFORT", "medium"),
            api_base=_env_str("LLM_API_BASE").rstrip("/"),
            api_key=_env_str("LLM_API_KEY"),
            model=_env_str("LLM_MODEL"),
            timeout=_env_int("LLM_TIMEOUT", 120),
            max_tokens=_env_int("LLM_MAX_TOKENS", 4000),
        )


@dataclass(frozen=True)
class Settings:
    ncbi_email: str = ""
    ncbi_api_key: str = ""
    data_dir: Path = Path("data")
    results_dir: Path = Path("data/results")
    criteria: Criteria = field(default_factory=Criteria)
    llm: LLMSettings = field(default_factory=LLMSettings)

    @property
    def raw_dir(self) -> Path:
        return self.data_dir / "raw"

    @property
    def interim_dir(self) -> Path:
        return self.data_dir / "interim"

    def ensure_dirs(self) -> None:
        for d in (self.raw_dir, self.interim_dir, self.results_dir):
            d.mkdir(parents=True, exist_ok=True)


def load_settings(env_file: str | os.PathLike[str] | None = ".env") -> Settings:
    """Lê o .env (se existir) e devolve as configurações do pipeline."""
    if env_file and Path(env_file).exists():
        load_dotenv(env_file, override=False)

    return Settings(
        ncbi_email=_env_str("NCBI_EMAIL"),
        ncbi_api_key=_env_str("NCBI_API_KEY"),
        data_dir=Path(_env_str("HOLINSEEK_DATA_DIR", "data")),
        results_dir=Path(_env_str("HOLINSEEK_RESULTS_DIR", "data/results")),
        criteria=Criteria.from_env(),
        llm=LLMSettings.from_env(),
    )
