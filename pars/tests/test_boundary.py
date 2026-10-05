"""The independence contract: src/pars must never import prrs (static check).

PARS is extractable only while this holds. The AST scan covers plain imports
and from-imports; the string convention the CLI uses (user-named
``module:attribute`` calculator factories) is resolved at runtime from user
input and therefore cannot create a static dependency.
"""

import ast
from pathlib import Path

SRC = Path(__file__).resolve().parents[1] / "src" / "pars"


def _imported_names(tree: ast.AST):
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            for alias in node.names:
                yield alias.name
        elif isinstance(node, ast.ImportFrom):
            if node.module and node.level == 0:
                yield node.module


def test_no_prrs_imports_anywhere_in_src():
    offenders = []
    for path in sorted(SRC.rglob("*.py")):
        tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
        for name in _imported_names(tree):
            if name == "prrs" or name.startswith("prrs."):
                offenders.append(f"{path.name}: {name}")
    assert offenders == [], f"pars must not import prrs; found: {offenders}"


def test_no_dynamic_import_of_prrs_strings():
    for path in sorted(SRC.rglob("*.py")):
        source = path.read_text(encoding="utf-8")
        assert 'import_module("prrs' not in source
        assert "from prrs" not in source
