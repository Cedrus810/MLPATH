"""Score PARS selection rules against pars_ground_truth.py records (exploration).

For each run: the refined barrier of every candidate, the distinct saddles among
them, and for each selection rule (the default, and the default with one cheap
objective left out) the representative count k, the barrier regret (best
representative barrier minus the best barrier of all candidates) and the share of
distinct saddles the representatives reach. Each rule is set against random
subsets of the same k, because a rule that does no better than chance has not
compressed anything. Also: Spearman rho between every cheap objective and the
refined barrier, i.e. whether the pre-refinement objectives carry the signal the
frontier assumes.

Saddles: top image of each converged band; two are the same saddle when their
Kabsch-aligned RMSD is below --saddle-rmsd AND their energies are within
--saddle-energy (single linkage). Atom permutations are not considered.

Post hoc by construction (see pars_ground_truth.py); nothing here is a criterion.

    python docs/experiments/pars_ground_truth_score.py --records runs/pars_ground_truth
"""

import argparse
import json
import sys
from pathlib import Path

import numpy as np
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "pars" / "src"))

import pars.frontier.archive as archive_module  # noqa: E402
from pars.config import ParsConfig  # noqa: E402
from pars.frontier.archive import ParetoArchive  # noqa: E402
from pars.frontier.objectives import evaluate_candidate, kabsch_align  # noqa: E402
from pars.paths.candidate import build_candidates  # noqa: E402
from pars.paths.ingest import ingest_run  # noqa: E402

CONVERGED = ("neb_converged", "ci_neb_converged")
OBJECTIVES = ("energy_excursion", "trigger_severity", "path_curvature", "recrossing")


def saddle_clusters(records: dict, rmsd_limit: float, energy_limit: float) -> dict:
    tops = {}
    for cid, record in records.items():
        if record["status"] in CONVERGED:
            energies = np.asarray(record["energies_eV"])
            top = int(np.argmax(energies))
            tops[cid] = (np.asarray(record["chain_A"][top]), float(energies[top]))
    ids = sorted(tops)
    label = {cid: index for index, cid in enumerate(ids)}
    for i, a in enumerate(ids):
        for b in ids[i + 1 :]:
            (pa, ea), (pb, eb) = tops[a], tops[b]
            rmsd = float(np.sqrt(np.mean(np.sum((kabsch_align(pb, pa) - pa) ** 2, axis=1))))
            if rmsd < rmsd_limit and abs(ea - eb) < energy_limit:
                old, new = label[b], label[a]
                for cid in ids:
                    if label[cid] == old:
                        label[cid] = new
    return label


def select(run_dir: Path, dropped: str | None) -> list[str]:
    config = ParsConfig()

    def evaluate(candidate, cfg):
        raw = evaluate_candidate(candidate, cfg)
        if dropped:
            raw[dropped] = None
        return raw

    archive_module.evaluate_candidate = evaluate
    try:
        candidates = build_candidates(ingest_run(run_dir, config), config)
        return ParetoArchive(candidates[0].channel_id, candidates, config).representatives
    finally:
        archive_module.evaluate_candidate = evaluate_candidate


def score(subset, barriers, labels, tolerance) -> tuple[float, float]:
    reached = [barriers[c] for c in subset if c in barriers]
    regret = (min(reached) - min(barriers.values())) if reached else float("inf")
    if regret < tolerance:  # differences below the refinement noise are not regret
        regret = 0.0
    coverage = len({labels[c] for c in subset if c in labels}) / len(set(labels.values()))
    return regret, coverage


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--records", required=True)
    parser.add_argument("--saddle-rmsd", type=float, default=0.1)
    parser.add_argument("--saddle-energy", type=float, default=0.05)
    parser.add_argument("--random-draws", type=int, default=20000)
    parser.add_argument(
        "--barrier-tolerance",
        type=float,
        default=0.005,
        help="eV; barriers closer than this count as equal (CI-NEB at fmax 0.05 "
        "reproduces one saddle to ~0.1 meV, so this is generous)",
    )
    args = parser.parse_args()
    rng = np.random.default_rng(0)
    for run_records in sorted(Path(args.records).iterdir()):
        records = {
            p.stem: json.loads(p.read_text()) for p in sorted(run_records.glob("*.json"))
        }
        if not records:  # e.g. logs/
            continue
        run_dir = ROOT / "runs" / run_records.name
        statuses = {}
        for record in records.values():
            statuses[record["status"]] = statuses.get(record["status"], 0) + 1
        barriers = {
            cid: r["barrier_eV"] for cid, r in records.items() if r["status"] in CONVERGED
        }
        labels = saddle_clusters(records, args.saddle_rmsd, args.saddle_energy)
        print(f"## {run_records.name}: {len(records)} candidates refined, statuses {statuses}")
        if not barriers:
            print("   no converged band; nothing to score")
            continue
        values = sorted(barriers.values())
        print(
            f"   barriers eV: min {values[0]:.3f} median {np.median(values):.3f} "
            f"max {values[-1]:.3f}; distinct saddles {len(set(labels.values()))}"
        )
        config = ParsConfig()
        candidates = build_candidates(ingest_run(run_dir, config), config)
        raw = {c.candidate_id: evaluate_candidate(c, config) for c in candidates}
        ids = sorted(barriers)
        if values[-1] - values[0] < args.barrier_tolerance:
            print(
                f"   barrier spread {1e3 * (values[-1] - values[0]):.2f} meV is below the "
                "tolerance: every candidate reaches the same barrier, rho is noise"
            )
        for name in OBJECTIVES if values[-1] - values[0] >= args.barrier_tolerance else ():
            xs = [raw[c][name] for c in ids]
            if None not in xs and len(set(xs)) > 1:
                rho = spearmanr(xs, [barriers[c] for c in ids]).statistic
                print(f"   rho({name}, barrier) = {rho:+.2f}")
        pool = sorted(records)
        for dropped in (None, *OBJECTIVES):
            reps = select(run_dir, dropped)
            regret, coverage = score(reps, barriers, labels, args.barrier_tolerance)
            k = len(reps)
            draws = [
                score(
                    list(rng.choice(pool, size=k, replace=False)),
                    barriers,
                    labels,
                    args.barrier_tolerance,
                )
                for _ in range(args.random_draws)
            ]
            random_regret = np.mean([d[0] for d in draws if np.isfinite(d[0])])
            p_best = np.mean([d[0] == 0.0 for d in draws])
            random_coverage = np.mean([d[1] for d in draws])
            print(
                f"   rule {'default' if dropped is None else 'without ' + dropped:28s} "
                f"k={k:2d}/{len(pool)}  regret {regret:.3f} eV (random {random_regret:.3f}, "
                f"P(best) random {p_best:.2f})  saddles {coverage:.2f} (random "
                f"{random_coverage:.2f})"
            )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
