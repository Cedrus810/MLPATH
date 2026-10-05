"""Electronic-state handoff: total charge and spin multiplicity.

Charge-aware potentials (MACE-POLAR, the b04 calculator) read the electronic
state from ``atoms.info["charge"]`` / ``["spin"]`` and silently fall back to a
neutral singlet when the keys are missing. The state is therefore a PES input
on the same footing as the model weights, and every structure PARS builds for a
calculator or writes to disk carries it explicitly. PRRS uses the same keys and
the same defaults (charge 0, multiplicity 1; see prrs/search.py).
"""

from __future__ import annotations

from numbers import Integral

KEYS = ("charge", "spin")
DEFAULTS = {"charge": 0, "spin": 1}


def _integer(value, what: str) -> int:
    if isinstance(value, bool) or not isinstance(value, Integral):
        raise ValueError(f"{what} must be an integer, got {value!r}")
    return int(value)


def checked_state(state, where: str) -> dict:
    """A validated ``{"charge": int, "spin": int >= 1}`` copy of ``state``."""
    if not isinstance(state, dict) or set(state) != set(KEYS):
        raise ValueError(
            f"{where}: electronic state must have exactly {list(KEYS)}, got {state!r}"
        )
    checked = {key: _integer(state[key], f"{where}: electronic state {key}") for key in KEYS}
    if checked["spin"] < 1:
        raise ValueError(f"{where}: electronic state spin multiplicity must be >= 1")
    return checked


def info_state(info: dict, where: str) -> dict:
    """Charge/spin recorded on one structure; a missing key stays None."""
    return {
        key: None if info.get(key) is None else _integer(info[key], f"{where}: {key}")
        for key in KEYS
    }


def stamp(info: dict, state: dict) -> None:
    info["charge"], info["spin"] = state["charge"], state["spin"]
