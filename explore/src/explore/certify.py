"""Common C2 certification backend (EVALUATION_PROTOCOL_DESIGN_2026-10-13 §1.1 C2 row,
§3.2 steps 2-5): saddle refinement, an independent full Hessian, descent on both sides
to confirmed endpoints, identity and event closure, duplicate check. Each step is one
ledger segment; the forwards of all of them are certification cost."""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from .poses import _jsonable


@dataclass
class Certification:
    status: str = "unresolved:saddle_not_found"
    saddle_positions: np.ndarray | None = None
    saddle_energy_eV: float | None = None
    unstable_mode: list | None = None
    hessian: dict | None = None
    sides: list = field(default_factory=list)
    event: dict | None = None
    duplicate_of: str | None = None
    diagnostics: dict = field(default_factory=dict)
    spent: int = 0
    chemical_pair: tuple | None = None


def _kabsch_rmsd(a, b):
    a = a - a.mean(axis=0)
    b = b - b.mean(axis=0)
    u, _, vt = np.linalg.svd(a.T @ b)
    d = np.sign(np.linalg.det(u @ vt))
    rot = u @ np.diag([1.0, 1.0, d]) @ vt
    return float(np.sqrt(np.mean(np.sum((a @ rot - b) ** 2, axis=1))))


def _pair(cert):
    return tuple(sorted(s["chemical_key"] for s in cert.sides))


def pair_symmetries(saddle, endpoints, config):
    """Atom permutations under which two saddles of one chemical pair are compared: the
    union of the graph automorphisms of the saddle and of its endpoint minima. Each one
    only permutes like atoms, so each is an exact symmetry of the PES. The saddle graph
    alone loses the minimum's symmetry when a moving atom bridges two groups (b05 E2:
    |Aut| 144 at the saddle, 1296 at the minima; equivalent saddles 1.673 A apart under
    the former, 0.0021 A under the latter -- DONE_LOG 2026-10-28). Returns
    (permutations, info), info being the first budget overflow met (else {})."""
    from prrs.chemistry import BUDGET_REASONS, automorphisms

    seen, perms, info = set(), [], {}
    for atoms in (saddle, *endpoints):
        group, ainfo = automorphisms(
            atoms, config.bond_scale, config.active_atoms, config.automorphism_limit
        )
        if not info and ainfo.get("reason") in BUDGET_REASONS:
            info = ainfo
        for p in group:
            key = tuple(int(i) for i in p)
            if key not in seen:
                seen.add(key)
                perms.append(p)
    return perms, info


def certify(
    guess,
    factory,
    config,
    ledger,
    task,
    seed,
    ceiling_eV,
    reference_eV,
    direction=None,
    known=(),
    attempt=0,
) -> Certification:
    """Certify the first-order saddle near `guess` as a C2 connection (or say why not).

    `known` is a sequence of `(id, Certification)`; a C2 among them with the same
    unordered chemical-key pair, saddle energy within `basin_energy_eV` and saddle
    geometry within `basin_rmsd_A` under the automorphisms of the saddle and endpoint
    graphs (`pair_symmetries`) makes this one a `duplicate_c2`. A descent that relayed through further saddles is `unresolved:relayed`
    (weaker evidence, protocol §3.2 step 4), never C2.

    v1 limitations: the negative-mode direction is not part of the duplicate rule
    (protocol §3.2 step 5 lists it; omitting it only makes merging easier). The endpoint
    minimum checks are charged to the "descent" stage because prrs `descend_saddle` runs
    them internally. With hessian_source="analytic", prrs runner.py ~684 swallows
    BudgetExhausted, so certify may report gate_rejected:analytic_hessian_unavailable
    instead of budget_exhausted."""
    from prrs.chemistry import (
        BUDGET_REASONS,
        IdentityUnavailable,
        chemical_key,
        symmetric_rmsd,
    )
    from prrs.reliability import GateRejected
    from prrs.runner import confirm_minimum, descend_saddle, follow_min_mode, side_seed
    from prrs.state import classify, same_basin

    from .ledger import BudgetExhausted

    guess = guess.copy()
    guess.set_constraint()
    cert = Certification()
    start = ledger.spent
    if config.minimum_check != "hessian":
        cert.status = "unresolved:no_hessian"
        return cert

    def seg(stage):
        return ledger.segment(task=task, stage=stage, attempt=attempt)

    try:
        # 1. saddle refinement
        with seg("saddle_refine") as s:
            saddle, info = follow_min_mode(
                guess,
                factory,
                config,
                direction=direction,
                seed=seed,
                ceiling=ceiling_eV,
                reference=reference_eV,
            )
            reason = str((info or {}).get("reason") or "")
            s.outcome = (
                "converged" if saddle is not None else f"not_found:{reason or 'unspecified'}"
            )
        cert.diagnostics["saddle_refine"] = _jsonable(info)
        if saddle is None:
            cert.status = (
                f"gate_rejected:{reason[5:]}"
                if reason.startswith("gate_")
                else f"unresolved:saddle_{reason or 'not_found'}"
            )
            return cert
        shift = saddle.positions - guess.positions
        if not isinstance(cert.diagnostics["saddle_refine"], dict):
            cert.diagnostics["saddle_refine"] = {"info": cert.diagnostics["saddle_refine"]}
        cert.diagnostics["saddle_refine"]["guess_to_saddle"] = {
            "max_atom_displacement_A": float(np.linalg.norm(shift, axis=1).max()),
            "rmsd_kabsch_A": _kabsch_rmsd(saddle.positions, guess.positions),
        }
        cert.saddle_positions = saddle.positions.copy()
        cert.saddle_energy_eV = float(saddle.get_potential_energy())

        # 2. independent full Hessian
        with seg("hessian") as s:
            _, checked = confirm_minimum(saddle, factory, config, seed)
            s.outcome = (
                f"order_{checked['saddle_order']}"
                if checked.get("saddle_order") is not None
                else str(checked.get("reason"))
            )
        cert.hessian = _jsonable(checked)
        reason = str(checked.get("reason", ""))
        order = checked.get("saddle_order")
        if reason.startswith("gate_"):
            cert.status = f"gate_rejected:{reason[5:]}"
            return cert
        if order is None:
            cert.status = f"unresolved:hessian_{reason}"
            return cert
        if order != 1:
            cert.status = f"rejected:order_{order}"
            return cert
        mode = checked.get("unstable_mode")
        if mode is None:
            cert.status = "unresolved:saddle_not_found"
            return cert
        cert.unstable_mode = _jsonable(mode)

        # 3. descent on both sides
        ends = {}
        with seg("descent") as s:
            for sign in (1, -1):
                endpoint, dchecked, relays = descend_saddle(
                    saddle,
                    np.asarray(mode),
                    factory,
                    config,
                    sign,
                    seed=side_seed(seed, sign),
                )
                ends[sign] = endpoint
                cert.sides.append(
                    {
                        "sign": sign,
                        "status": "minimum" if endpoint is not None else "failed",
                        "positions": None if endpoint is None else endpoint.positions.copy(),
                        "energy_eV": None
                        if endpoint is None
                        else float(endpoint.get_potential_energy()),
                        "relays": _jsonable(relays),
                        "chemical_key": None,
                        "diagnostics": _jsonable(dchecked),
                    }
                )
            s.outcome = "done"
        if any(e is None for e in ends.values()):
            cert.status = "unresolved:descent_failed"
            return cert
        if any(side["relays"] for side in cert.sides):
            cert.status = "unresolved:relayed"
            return cert

        # 4. identity and event
        minus, plus = ends[-1], ends[1]
        if same_basin(minus, plus, config):
            cert.status = "unresolved:both_sides_same_basin"
            return cert
        label, fingerprint = classify(minus, plus, config)
        cert.event = {"label": label, "fingerprint": _jsonable(fingerprint)}
        try:
            for side in cert.sides:
                end = ends[side["sign"]]
                side["chemical_key"] = chemical_key(
                    end,
                    config.bond_scale,
                    config.active_atoms,
                    charge_sensitive=config.charge_sensitive,
                    charge=config.total_charge,
                    multiplicity=config.multiplicity,
                    parity_tolerance=config.parity_tolerance,
                )["key"]
        except IdentityUnavailable as exc:
            cert.status = "unresolved:identity"
            cert.diagnostics["identity"] = str(exc)
            return cert
        cert.chemical_pair = _pair(cert)

        # 5. duplicates
        group = None
        for known_id, other in known:
            if other.status != "c2" or (
                other.chemical_pair is None
                or tuple(other.chemical_pair) != tuple(cert.chemical_pair)
            ):
                continue
            if abs(other.saddle_energy_eV - cert.saddle_energy_eV) > config.basin_energy_eV:
                continue
            if group is None:
                group, ainfo = pair_symmetries(saddle, (minus, plus), config)
                if ainfo.get("reason") in BUDGET_REASONS:
                    cert.status = "unresolved:identity"
                    cert.diagnostics["automorphisms"] = _jsonable(ainfo)
                    return cert
            other_atoms = saddle.copy()
            other_atoms.positions = other.saddle_positions
            if (
                symmetric_rmsd(saddle, other_atoms, group, config.active_atoms)
                <= config.basin_rmsd_A
            ):
                cert.status = "duplicate_c2"
                cert.duplicate_of = known_id
                return cert
        cert.status = "c2"
        return cert
    except BudgetExhausted:
        cert.status = "budget_exhausted"
        return cert
    except GateRejected as exc:
        cert.status = f"gate_rejected:{exc.code}"
        return cert
    finally:
        cert.spent = ledger.spent - start
