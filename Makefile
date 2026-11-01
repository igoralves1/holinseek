.PHONY: help setup test demo lint site clean

PY ?= python3
VENV = .venv
BIN = $(VENV)/bin

help:
	@echo "setup   cria o ambiente virtual e instala tudo"
	@echo "test    roda a suíte de testes"
	@echo "demo    roda o pipeline sobre o genoma sintético de exemplo"
	@echo "site    serve o site em http://localhost:8000"
	@echo "clean   remove caches e resultados intermediários"

setup:
	$(PY) -m venv $(VENV)
	$(BIN)/pip install --upgrade pip
	$(BIN)/pip install -r requirements.txt
	$(BIN)/pip install -e .
	@test -f .env || cp .env.example .env
	@echo "pronto. ative com: source $(BIN)/activate"

test:
	$(BIN)/python -m pytest -q

demo:
	$(BIN)/holinseek demo

site:
	@echo "http://localhost:8000"
	cd docs && $(PY) -m http.server 8000

clean:
	rm -rf .pytest_cache **/__pycache__ src/**/__pycache__
	rm -rf data/interim/* data/results/*
