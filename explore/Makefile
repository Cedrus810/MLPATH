# Same two-environment discipline as the parent repository and pars/: the
# science interpreter only ever runs tests and the module itself, formatting
# lives in the parent's .venv-tools so nothing can move torch/mace/openmm.
# RUFF points at the parent checkout by default; adjust it after extracting
# this folder.
PYTHON ?= /home/ruigengji/miniforge3/envs/openmm_dev/bin/python
RUFF   ?= ../.venv-tools/bin/ruff

.PHONY: help fmt fmt-check lint test check

help:
	@echo "make fmt        format src/explore and tests/"
	@echo "make fmt-check  verify formatting"
	@echo "make lint       ruff check"
	@echo "make test       pytest (needs the science environment)"
	@echo "make check      fmt-check + lint + test, in that order"

fmt:
	$(RUFF) format src/explore tests
	$(RUFF) check --fix src/explore tests

fmt-check:
	$(RUFF) format --check src/explore tests

lint:
	$(RUFF) check src/explore tests

test:
	$(PYTHON) -m pytest -q

check: fmt-check lint test
