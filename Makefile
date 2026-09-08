# Two environments on purpose. PYTHON is the science environment and its versions are
# part of every result -- backend parity went red once on an openmm point release, so
# nothing that is merely a developer convenience may be installed into it. RUFF lives in
# .venv-tools, which exists only so formatting can never move openmm, torch or mace.
PYTHON ?= /home/ruigengji/miniforge3/envs/openmm_dev/bin/python
RUFF   ?= .venv-tools/bin/ruff

.PHONY: help tools fmt fmt-check lint test digest check hash

help:
	@echo "make tools      create .venv-tools and install ruff there"
	@echo "make fmt        format src/ and tests/"
	@echo "make lint       ruff check"
	@echo "make test       pytest (needs the science environment)"
	@echo "make digest     the golden-sample gate; run before and after touching src/prrs"
	@echo "make check      lint + test + digest, in that order"
	@echo "make hash       the implementation_sha256 a run manifest would record"

tools:
	$(PYTHON) -m venv .venv-tools
	.venv-tools/bin/pip install -q ruff

fmt:
	$(RUFF) format src/prrs tests
	$(RUFF) check --fix src/prrs tests

fmt-check:
	$(RUFF) format --check src/prrs tests

lint:
	$(RUFF) check src/prrs tests

test:
	$(PYTHON) -m pytest -q

# Not part of `test` because it is slower and answers a different question: pytest asks
# whether the code is correct, this asks whether it still produces the same numbers as
# the runs already on disk. A change that is meant to move them passes --allow.
digest:
	$(PYTHON) docs/experiments/digest_regression.py

check: lint test digest

hash:
	@$(PYTHON) scripts/implementation_hash.py
