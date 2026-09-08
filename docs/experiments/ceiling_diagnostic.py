"""What is the 1.5 eV climb ceiling actually blocking?

`ceiling_blocked` was the most common way a min-mode climb ended, and after two acceptance
checker defects were fixed it became the whole of the only remaining real failure in the
P2 T=300 batch: seed 71's transfer channel has a directed observation but no transition
state, because two of its three climb seeds hit the ceiling and the third lost its target
mode.

This is a MEASUREMENT, not a tuning run. It answers one question -- where would those two
climbs have gone -- and it is not allowed to answer "what is the best ceiling".

Held fixed: the seed frame, the target direction, the model, and every other config field.
The only change is the ceiling passed to `follow_min_mode`, raised to a diagnostic 3 eV.
The reliability gates inside `GuardedCalculator` (force limit, energy rise per atom,
uncertainty) stay exactly as they were, because they are not the ceiling. The domain gate
acts at admission, and nothing here is admitted, so it is instead EVALUATED and reported on
whatever the climb reaches.

Nothing is registered: no node, no channel, no criterion, no acceptance verdict. The output
is a report.

Three outcomes are worth distinguishing, and they lead to different conclusions:

  a  a correct proton-transfer saddle at a lower energy than the ceiling would have allowed
     -> the climb was taking a detour; the problem is path or step length, not the ceiling
  b  an O-H homolysis or some other out-of-domain structure at 2-3 eV
     -> the 1.5 eV ceiling is protecting correctly and must not be raised
  c  the target mode is lost, or the energy just keeps rising
     -> a mode-following robustness problem, which raising the ceiling would not fix
"""
import json
import sys
from dataclasses import replace
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))
from prrs.calculators import load_factory                            # noqa: E402
from prrs.chemistry import chemical_key                              # noqa: E402
from prrs.config import SearchConfig                                 # noqa: E402
from prrs.runner import PathRecorder, follow_min_mode                # noqa: E402
from prrs.search import _dividing_frame                              # noqa: E402
from prrs.state import encode                                        # noqa: E402
from prrs import internal                                            # noqa: E402

RUN = Path("runs/p2_T300_seed71")
BLOCKED = ("t000001", "t000026")
DIAGNOSTIC_CEILING_eV = 3.0
OUT = Path("runs/ceiling_diagnostic")


def domain_verdict(atoms, cfg):
    """What the domain gate WOULD say. Evaluated, never applied -- nothing is admitted."""
    components = chemical_key(atoms, cfg.bond_scale, cfg.active_atoms,
                              bond_order_limit=cfg.automorphism_limit,
                              charge_sensitive=cfg.charge_sensitive,
                              charge=cfg.total_charge, multiplicity=cfg.multiplicity)
    verdict = components.get("stereo_unresolved")
    refused = verdict in ("no_assignment_satisfies_valences", "degree_exceeds_valence")
    return {"chemical_key": components["key"][:12], "fragments": components["fragments"],
            "graph_hash": components["graph_hash"], "stereo_unresolved": verdict,
            "would_be_refused_as_out_of_domain": bool(refused)}


def describe(atoms, cfg, source_energy):
    graph = encode(atoms, cfg.bond_scale, active=cfg.active_atoms)
    tracked = {spec["name"]: internal.named_values(atoms.positions,
                                                   [spec])[spec["name"]]
               for spec in cfg.tracked_coordinates}
    return {"energy_eV": float(atoms.get_potential_energy()),
            "rise_above_source_eV": float(atoms.get_potential_energy()) - source_energy,
            "tracked": {k: round(v, 4) for k, v in tracked.items()},
            "edges": [list(e) for e in sorted(graph.edges)],
            "fragments_count": graph.components,
            "domain": domain_verdict(atoms, cfg)}


def main():
    manifest = json.loads((RUN / "config.json").read_text())
    cfg = SearchConfig.from_dict(manifest["config"])
    factory = load_factory(manifest["calculator"]["factory"],
                           manifest["calculator"]["kwargs"])
    network = json.loads((RUN / "network.json").read_text())
    source_energy = network["chemical_nodes"][0]["microstates"][0]["energy_eV"]
    OUT.mkdir(parents=True, exist_ok=True)

    print(f"source E = {source_energy:.6f} eV")
    print(f"production ceiling = {cfg.saddle_search_max_rise_eV} eV above it")
    print(f"diagnostic ceiling = {DIAGNOSTIC_CEILING_eV} eV above it")
    print("every other field identical; nothing is registered\n")

    report = {"source_energy_eV": source_energy,
              "production_ceiling_eV": cfg.saddle_search_max_rise_eV,
              "diagnostic_ceiling_eV": DIAGNOSTIC_CEILING_eV,
              "implementation_sha256": manifest["implementation_sha256"],
              "meaning": ("a measurement of what the ceiling blocks; it does not propose a "
                          "ceiling and nothing here is registered"),
              "seeds": []}

    for attempt in BLOCKED:
        record = json.loads((RUN / "trials" / attempt / "result.json").read_text())
        frame, origin, tangent, source_name = _dividing_frame(
            RUN / "trials" / attempt, probe=record.get("probe"))
        entry = {"seed_attempt": attempt, "seed_frame": origin,
                 "direction_source": source_name,
                 "probe": {k: record["probe"][k]
                           for k in ("family", "indices", "amplitude", "unit")}}
        if frame is None or tangent is None:
            entry["reason"] = "seed frame or tangent unavailable"
            report["seeds"].append(entry)
            print(f"{attempt}: {entry['reason']}")
            continue
        frame.calc = factory()
        entry["seed"] = describe(frame, cfg, source_energy)
        print(f"=== {attempt}  {record['probe']['family']}{record['probe']['indices']} "
              f"amp={record['probe']['amplitude']}  frame={origin}")
        print(f"    seed sits {entry['seed']['rise_above_source_eV']:+.4f} eV above the "
              f"source, so the production ceiling left "
              f"{cfg.saddle_search_max_rise_eV - entry['seed']['rise_above_source_eV']:+.4f} eV")

        for label, ceiling in (("control_1.5", cfg.saddle_search_max_rise_eV),
                               ("diagnostic_3.0", DIAGNOSTIC_CEILING_eV)):
            recorder = PathRecorder(OUT / f"{attempt}_{label}", cfg, "climb")
            structure, diagnostics = follow_min_mode(
                frame, factory, cfg, direction=tangent,
                ceiling=source_energy + ceiling, reference=source_energy,
                recorder=recorder)
            arm = {"ceiling_eV": ceiling,
                   "reason": diagnostics.get("reason"),
                   "steps": diagnostics.get("steps"),
                   "target_overlap": diagnostics.get("target_overlap"),
                   "selected_mode_index": diagnostics.get("selected_mode_index"),
                   "initial_rise_eV": diagnostics.get("initial_rise_eV"),
                   "max_rise_eV": diagnostics.get("max_rise_eV"),
                   "backtracks": diagnostics.get("backtracks"),
                   "ceiling_hit_step": diagnostics.get("ceiling_hit_step"),
                   "saddle_order": diagnostics.get("saddle_order"),
                   "imaginary_wavenumbers_icm":
                       diagnostics.get("imaginary_wavenumbers_icm"),
                   "final_mode_overlap": diagnostics.get("final_mode_overlap")}
            if structure is not None:
                structure.calc = factory()
                arm["reached"] = describe(structure, cfg, source_energy)
            entry[label] = arm
            print(f"    {label:<15} reason={arm['reason']} steps={arm['steps']} "
                  f"max_rise={arm['max_rise_eV']} backtracks={arm['backtracks']} "
                  f"hit_at={arm['ceiling_hit_step']}")
            if structure is not None:
                r = arm["reached"]
                print(f"    {'':15} REACHED order={arm['saddle_order']} "
                      f"icm={[round(v, 1) for v in arm['imaginary_wavenumbers_icm'] or []]} "
                      f"rise={r['rise_above_source_eV']:+.4f} eV")
                print(f"    {'':15} {r['tracked']}  frags={r['fragments_count']} "
                      f"{r['domain']['fragments']} "
                      f"out_of_domain={r['domain']['would_be_refused_as_out_of_domain']}")
        report["seeds"].append(entry)
        print()

    (OUT / "report.json").write_text(json.dumps(report, indent=1), encoding="utf-8")
    print(f"wrote {OUT / 'report.json'}")


if __name__ == "__main__":
    main()
