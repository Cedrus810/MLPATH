"""Build draft tables from saved HF-6 summaries; no model or molecular jobs.

Run from any directory with Python 3.10+. Figures additionally need matplotlib.
Historical summaries are read only. The saved mirror audit is an explicit overlay.
"""

import csv
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
PAPER = ROOT / "paper"
RES = ROOT / "benchmarks/explore_d2/results"
SOURCES = {
    "local": RES / "2026-11-10/hf6_local_summary_v3.json",
    "sA": RES / "2026-11-11/hf6_e2e_summary_sA_v3.json",
    "sB": RES / "2026-11-11/hf6_e2e_summary_sB_v2.json",
    "mirror": RES / "2026-11-12/hf6_mirror_match_check.json",
    "identity": RES / "2026-11-12/dimer_identity/summary.json",
}
SYSTEMS = ("b04", "b05", "n1")


def write_csv(name, rows):
    with (PAPER / "data" / name).open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def build():
    data = {k: json.loads(p.read_text()) for k, p in SOURCES.items()}
    local, candidates, e2e, hits = [], [], [], []
    applied = []
    for system in SYSTEMS:
        spec = data["local"][system]
        overlay = {
            (r["item"], r["arm"]): r["mirror_match"] for r in data["mirror"]["local"][system]
        }
        for arm in ("hf", "dimer"):
            rows = []
            for item in spec["rows"]:
                cell = item[arm]
                key = (item["item"], arm)
                match = cell.get("match")
                if key in overlay:
                    if match or cell.get("order") != 1:
                        raise ValueError("Mirror overlay no longer fits source summary")
                    match = overlay[key]
                    applied.append(
                        {"system": system, "item": key[0], "arm": arm, "mirror_match": match}
                    )
                row = {
                    "system": system,
                    "arm": arm,
                    "item": item["item"],
                    "status": cell["status"],
                    "forwards": cell["forwards"],
                    "candidate": cell["candidate"],
                    "order": cell.get("order"),
                    "order_reason": cell.get("order_reason"),
                    "original_match": cell.get("match"),
                    "draft_match": match,
                    "strict_hit": cell.get("order") == 1 and bool(match),
                    "mirror_overlay": key in overlay,
                }
                rows.append(row)
            candidates.extend(rows)
            a = spec["arms"][arm]
            forwards = sum(r["forwards"] for r in rows)
            if forwards != a["forwards_total"]:
                raise ValueError("Per-item costs disagree with saved summary")
            local.append(
                {
                    "system": system,
                    "arm": arm,
                    "items": len(rows),
                    "candidates": sum(r["candidate"] for r in rows),
                    "order1": sum(r["order"] == 1 for r in rows),
                    "geometry_matches": sum(bool(r["draft_match"]) for r in rows),
                    "strict_reference_hits": sum(r["strict_hit"] for r in rows),
                    "distinct_strict_references": len(
                        {r["draft_match"] for r in rows if r["strict_hit"]}
                    ),
                    "forwards": forwards,
                    "source": str(SOURCES["local"].relative_to(ROOT)),
                }
            )
    for seed in ("sA", "sB"):
        for system in SYSTEMS:
            spec = data[seed][system]
            for arm, a in spec["arms"].items():
                if data["mirror"]["e2e"][seed][system]:
                    raise ValueError("New end-to-end mirror changes require review")
                found = [cost for cost in a["first_hit"].values() if cost is not None]
                if any(cost > spec["budget"] for cost in found):
                    raise ValueError("Reference hit falls outside recorded budget")
                e2e.append(
                    {
                        "seed": seed,
                        "system": system,
                        "arm": arm,
                        "budget": spec["budget"],
                        "reference_count": len(spec["references"]),
                        "strict_references_found": len(found),
                        "first_strict_hit": min(found) if found else None,
                        "first_hit_budget_fraction": min(found) / spec["budget"]
                        if found
                        else None,
                        "order1": a["order_1"],
                        "source": str(SOURCES[seed].relative_to(ROOT)),
                    }
                )
                for reference, cost in a["first_hit"].items():
                    hits.append(
                        {
                            "seed": seed,
                            "system": system,
                            "arm": arm,
                            "reference": reference,
                            "cost": cost,
                            "budget": spec["budget"],
                            "hit": cost is not None,
                        }
                    )
    (PAPER / "data").mkdir(parents=True, exist_ok=True)
    write_csv("local_results.csv", local)
    write_csv("local_items.csv", candidates)
    write_csv("end_to_end_results.csv", e2e)
    write_csv("reference_hits.csv", hits)
    evidence = {
        "scope": "development cohorts; saved validation; no new scientific calculations",
        "sources": {k: str(v.relative_to(ROOT)) for k, v in SOURCES.items()},
        "mirror_overlay_applied": applied,
        "local": local,
        "end_to_end": e2e,
        "identity_categories": {
            s: {k: v for k, v in data["identity"][s].items() if k != "rows"} for s in SYSTEMS
        },
    }
    (PAPER / "data/draft_evidence.json").write_text(json.dumps(evidence, indent=2) + "\n")
    return local, e2e, hits


def figures(local, e2e, hits):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out = PAPER / "figures"
    out.mkdir(exist_ok=True)
    colors = {"hf": "#277DA1", "dimer": "#D47736"}
    plt.rcParams.update({"font.size": 10, "pdf.fonttype": 42, "ps.fonttype": 42})
    fig, axes = plt.subplots(1, 2, figsize=(10, 4), layout="constrained")
    for index, arm in enumerate(("hf", "dimer")):
        cells = [r for r in local if r["arm"] == arm]
        x = [i + (-0.19 if index == 0 else 0.19) for i in range(3)]
        bars = axes[0].bar(
            x,
            [r["forwards"] for r in cells],
            width=0.36,
            color=colors[arm],
            label="HF" if arm == "hf" else "ASE Dimer",
        )
        axes[0].bar_label(
            bars, labels=[f"{r['forwards']:,}" for r in cells], fontsize=8, padding=3
        )
        bars = axes[1].bar(
            x, [r["strict_reference_hits"] for r in cells], width=0.36, color=colors[arm]
        )
        axes[1].bar_label(bars, padding=3)
    for ax in axes:
        ax.set_xticks(range(3), ["b04 (23)", "b05 (140)", "N1 (146)"])
        ax.spines[["top", "right"]].set_visible(False)
    axes[0].set_yscale("log")
    axes[0].set_ylim(2000, 1000000)
    axes[0].set_ylabel("Total search forward evaluations")
    axes[0].legend(frameon=False)
    axes[0].set_title("a  Cost on frozen bracket sets", loc="left")
    axes[1].set_ylabel("Index-one, reference-matching attempts")
    axes[1].set_ylim(0, 20)
    axes[1].set_yticks([0, 5, 10, 15, 20])
    axes[1].set_title("b  Repeated hits, not unique channels", loc="left")
    for extension in ("png", "pdf"):
        fig.savefig(out / f"Fig1_local.{extension}", dpi=220)
    plt.close(fig)
    fig, axes = plt.subplots(1, 3, figsize=(11, 3.6), layout="constrained")
    for ax, system in zip(axes, SYSTEMS):
        for seed, style in (("sA", "-"), ("sB", "--")):
            for arm, color in (("d2", colors["hf"]), ("db", colors["dimer"])):
                row = next(
                    r
                    for r in e2e
                    if r["seed"] == seed and r["system"] == system and r["arm"] == arm
                )
                costs = sorted(
                    r["cost"]
                    for r in hits
                    if r["seed"] == seed
                    and r["system"] == system
                    and r["arm"] == arm
                    and r["hit"]
                )
                x = [0] + [c / row["budget"] for c in costs] + [1]
                y = [0] + [(i + 1) / row["reference_count"] for i in range(len(costs))]
                y.append(y[-1])
                ax.step(
                    x,
                    y,
                    where="post",
                    color=color,
                    linestyle=style,
                    label=f"{'D2' if arm == 'd2' else 'DB'} {seed}",
                    alpha=0.85,
                )
        ax.set_title(system.upper() if system == "n1" else system)
        ax.set_xlim(0, 1.03)
        ax.set_ylim(-0.06, 1.08)
        ax.set_xlabel("Fraction of recorded budget B*")
        ax.spines[["top", "right"]].set_visible(False)
    axes[0].set_ylabel("Strict reference-set recall")
    axes[0].legend(fontsize=8, frameon=False)
    fig.suptitle(
        "Development runs: two seeds; reference-set recall excludes novel channels", fontsize=11
    )
    for extension in ("png", "pdf"):
        fig.savefig(out / f"Fig2_recall.{extension}", dpi=220)
    plt.close(fig)


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--figures", action="store_true")
    args = parser.parse_args()
    local, e2e, hits = build()
    if args.figures:
        figures(local, e2e, hits)
    print(
        f"Built {len(local)} local rows, {len(e2e)} end-to-end rows and {len(hits)} reference records."
    )
    print("Applied saved mirror audit to N1 local Dimer; historical summaries unchanged.")
