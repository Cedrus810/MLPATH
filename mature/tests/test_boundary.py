"""The module contract: nothing imports ``mature`` back (static check).

``mature`` is allowed -- and required -- to import ``prrs`` one-way: the
construction plan mandates reusing ``runner.shoot`` and the read-only
``Registry.identify`` instead of copying the dynamics implementation.
The boundary that must hold is the reverse direction: neither the PRRS
tree (``src/prrs``) nor the PARS tree (``pars/src/pars``) may depend on
the maturation layer, statically or via dynamic-import strings. PARS's
own no-``prrs`` rule is re-asserted here from the maturation side so a
future edit on either side cannot break it unnoticed.

The sibling trees are scanned only when present: after standalone
extraction these checks skip, they do not silently pass on nothing.
"""

import ast
from pathlib import Path

import pytest

TESTS = Path(__file__).resolve().parent
PARENT = TESTS.parent.parent  # repository root; absent after extraction
PRRS_SRC = PARENT / "src" / "prrs"
PARS_SRC = PARENT / "pars" / "src" / "pars"


def _imported_names(tree: ast.AST):
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            for alias in node.names:
                yield alias.name
        elif isinstance(node, ast.ImportFrom):
            if node.module and node.level == 0:
                yield node.module


def _offenders(root: Path, banned: str):
    found = []
    for path in sorted(root.rglob("*.py")):
        source = path.read_text(encoding="utf-8")
        tree = ast.parse(source, filename=str(path))
        for name in _imported_names(tree):
            if name == banned or name.startswith(banned + "."):
                found.append(f"{path}: import {name}")
        if f'import_module("{banned}' in source:
            found.append(f'{path}: import_module("{banned}...")')
        if f"from {banned}" in source:
            found.append(f"{path}: 'from {banned}'")
    return found


@pytest.mark.skipif(not PRRS_SRC.is_dir(), reason="standalone extraction: no sibling prrs tree")
def test_prrs_never_imports_mature():
    assert _offenders(PRRS_SRC, "mature") == []


@pytest.mark.skipif(not PARS_SRC.is_dir(), reason="standalone extraction: no sibling pars tree")
def test_pars_never_imports_mature():
    assert _offenders(PARS_SRC, "mature") == []


@pytest.mark.skipif(not PARS_SRC.is_dir(), reason="standalone extraction: no sibling pars tree")
def test_pars_prrs_boundary_still_holds():
    # mirror of pars/tests/test_boundary.py, enforced from this side too
    assert _offenders(PARS_SRC, "prrs") == []
