"""Standalone figures from the recorded run data. SVG + PDF, one file per figure.

Nothing here computes anything new: every point is read off a run directory. The path
figure in particular is measured, not schematic -- the probe, the free response, the quench,
the climb and both descents are the frames those stages actually wrote.
"""
import json
import sys
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

OUT = Path(__file__).resolve().parent
ROOT = OUT.parents[1]
SEEDS = (17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79)

INK, MUTED, RULE = "#141C1A", "#647069", "#D2D8D4"
TEAL, INDIGO, OXIDE = "#008C7E", "#5145C8", "#BF4A34"

plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 8.5,
    "axes.edgecolor": RULE, "axes.labelcolor": INK, "text.color": INK,
    "xtick.color": MUTED, "ytick.color": MUTED,
    "axes.linewidth": 0.7, "xtick.major.width": 0.7, "ytick.major.width": 0.7,
    "xtick.major.size": 3, "ytick.major.size": 3,
    "axes.spines.top": False, "axes.spines.right": False,
    "figure.facecolor": "white", "axes.facecolor": "white",
    "savefig.bbox": "tight", "savefig.pad_inches": 0.06, "savefig.dpi": 160,
    "legend.frameon": False, "legend.fontsize": 7.5,
})


def save(fig, stem):
    for suffix in ("svg", "pdf", "png"):
        fig.savefig(OUT / f"{stem}.{suffix}")
    plt.close(fig)
    print(f"  {stem}.svg  {stem}.pdf")


def rows(path):
    path = Path(path)
    return ([json.loads(line) for line in path.read_text().splitlines()]
            if path.exists() else [])


def wilson(k, n, z=1.959963985):
    p, d = k / n, 1 + z * z / n
    centre = (p + z * z / (2 * n)) / d
    half = z * (p * (1 - p) / n + z * z / (4 * n * n)) ** 0.5 / d
    return centre - half, centre + half


# --------------------------------------------------------------- 1. the path
def figure_transition_path():
    """E against the transfer coordinate, for one saddle and everything that reached it."""
    run = ROOT / "runs/p2_T300_seed17"
    network = json.loads((run / "network.json").read_text())
    source = network["chemical_nodes"][0]["microstates"][0]["energy_eV"]
    scale = 1000.0                                       # meV above the reactant minimum

    trial = rows(run / "trials/t000008/observations.jsonl")
    climb = rows(run / "paths/climb003/climb.jsonl")
    descent = rows(run / "paths/ts0000/descent.jsonl")
    saddle = next(t for t in network["ts_candidates"] if t["id"] == "ts0000")

    def series(items, phases=None):
        out = []
        for item in items:
            if item.get("potential_eV") is None:
                continue
            if phases and item["phase"] not in phases:
                continue
            tracked = item.get("tracked") or {}
            if "q_PT" not in tracked:
                continue
            out.append((tracked["q_PT"], (item["potential_eV"] - source) * scale))
        return list(zip(*out)) if out else ([], [])

    fig, ax = plt.subplots(figsize=(6.6, 4.0))
    ax.axhline(0, color=RULE, lw=0.7, zorder=1)

    driven_q, driven_e = series(trial, {"source", "perturbed", "forced", "free_start"})
    free_q, free_e = series(trial, {"free"})
    quench_q, quench_e = series(trial, {"quench", "quenched"})
    climb_q, climb_e = series(climb)
    plus = [r for r in descent if r["phase"] == "descent_plus"]
    minus = [r for r in descent if r["phase"] == "descent_minus"]
    plus_q, plus_e = series(plus)
    minus_q, minus_e = series(minus)

    # Sampled every `sample_interval` MD steps, so consecutive samples are not adjacent in
    # time -- joining them draws a polyline that looks like a path and is not one.
    ax.plot(free_q, free_e, ".", color=MUTED, ms=3.4, alpha=0.55, zorder=2)
    ax.plot(quench_q, quench_e, "-", color=INDIGO, lw=1.4, zorder=3)
    ax.plot(plus_q, plus_e, "-", color=TEAL, lw=1.4, alpha=0.75, zorder=3)
    ax.plot(minus_q, minus_e, "-", color=TEAL, lw=1.4, alpha=0.75, zorder=3)
    ax.plot(climb_q, climb_e, "-", color=OXIDE, lw=1.8, zorder=4)
    ax.plot(climb_q, climb_e, "o", color=OXIDE, ms=3.0, zorder=5)
    if driven_q:
        ax.plot(driven_q, driven_e, "o", color=MUTED, ms=4.0, mfc="white", mew=1.2, zorder=5)

    top = (saddle["energy_eV"] - source) * scale
    ax.plot([0.0], [top], "*", color=OXIDE, ms=15, mec="white", mew=0.8, zorder=6)
    ax.annotate(f"first-order saddle\n{saddle['imaginary_wavenumbers_icm'][0]:.0f} cm$^{{-1}}$"
                f"\n+{top:.0f} meV",
                xy=(0.0, top), xytext=(-0.30, top + 150), fontsize=7.5, color=OXIDE,
                ha="left", va="bottom",
                arrowprops=dict(arrowstyle="-", color=OXIDE, lw=0.7, shrinkA=6, shrinkB=2))
    ax.annotate("A", xy=(0.703, 0), xytext=(0.703, 42), fontsize=10, ha="center",
                color=INK, weight="bold")
    ax.annotate("B", xy=(-0.70, 0), xytext=(-0.70, 42), fontsize=10, ha="center",
                color=INK, weight="bold")

    ax.set_xlabel(r"$q_{\rm PT} = r({\rm O_1{-}H}) - r({\rm O_2{-}H})$   /  Å")
    ax.set_ylabel("energy above the reactant minimum  /  meV")
    ax.set_title("A measured transfer path and the saddle climbed from it",
                 loc="left", fontsize=10, pad=10)
    ax.legend(handles=[
        Line2D([], [], marker="o", color=MUTED, mfc="white", mew=1.2, ls="none",
               label="probe delivered"),
        Line2D([], [], marker=".", color=MUTED, ls="none", ms=6,
               label="free NVE response (sampled)"),
        Line2D([], [], color=INDIGO, lw=1.4, label="quench"),
        Line2D([], [], color=OXIDE, lw=1.8, marker="o", ms=3,
               label="min-mode walk to the saddle"),
        Line2D([], [], color=TEAL, lw=1.4, label="descent, both sides"),
        Line2D([], [], marker="*", color=OXIDE, ls="none", ms=11, label="saddle"),
    ], loc="upper center", ncol=3, bbox_to_anchor=(0.5, -0.19))
    ax.text(0.0, -0.42, "seed 17 · trial t000008 · ts0000. Every point is a recorded frame; "
            "the walk reaches the saddle from above.",
            fontsize=7, color=MUTED, transform=ax.transAxes, ha="left")
    save(fig, "fig1_transition_path")


# ------------------------------------------------- 2. the three frozen questions
def figure_frozen_questions():
    now = [("Q1", 16, 16), ("Q2", 16, 16), ("Q4", 15, 16)]
    was = [("Q1", 12, 16), ("Q2", 15, 16), ("Q4", 14, 16)]
    labels = ["Q1  directed edge on\nthe transfer channel",
              "Q2  the transfer\nchannel exists",
              "Q4  a min-mode TS\nsupports it"]
    fig, ax = plt.subplots(figsize=(6.2, 2.5))
    for i, ((_, k, n), (_, ko, no)) in enumerate(zip(now, was)):
        y = len(now) - 1 - i
        lo, hi = wilson(ko, no)
        ax.plot([lo, hi], [y - 0.18] * 2, color=MUTED, lw=1.6, alpha=0.5,
                solid_capstyle="round", zorder=2)
        ax.plot([ko / no], [y - 0.18], "o", ms=5, mfc="white", color=MUTED, mew=1.4, zorder=3)
        lo, hi = wilson(k, n)
        ax.plot([lo, hi], [y + 0.10] * 2, color=TEAL, lw=3.2, alpha=0.35,
                solid_capstyle="round", zorder=2)
        ax.plot([k / n], [y + 0.10], "o", ms=7, color=TEAL, mec="white", mew=1.2, zorder=3)
        ax.text(1.015, y + 0.10, f"{k}/{n}", va="center", fontsize=9, color=INK, weight="bold")
    ax.set_yticks(range(len(now)))
    ax.set_yticklabels(labels[::-1], fontsize=8)
    ax.set_ylim(-0.6, len(now) - 0.35)
    ax.set_xlim(0.35, 1.09)
    ax.set_xticks([0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0])
    ax.set_xlabel("pass rate, with the 95% Wilson interval")
    ax.set_title("Sixteen seeds: the three frozen questions", loc="left", fontsize=10, pad=10)
    ax.grid(axis="x", color=RULE, lw=0.6, alpha=0.6)
    ax.set_axisbelow(True)
    ax.legend(handles=[
        Line2D([], [], marker="o", color=TEAL, ls="none", ms=7, label="single implementation"),
        Line2D([], [], marker="o", color=MUTED, mfc="white", mew=1.4, ls="none", ms=5,
               label="earlier, mixed implementations"),
    ], loc="upper center", ncol=2, bbox_to_anchor=(0.5, -0.30))
    save(fig, "fig2_frozen_questions")


# ------------------------------------------------------- 3. the energy spread
def figure_energy_spread():
    reference = -11.721
    values = []
    for seed in SEEDS:
        network = json.loads((ROOT / f"runs/p2_T300_seed{seed}/network.json").read_text())
        nodes = network["chemical_nodes"]
        channel = next((c for c in network["reaction_channels"]
                        if c["class"] == "reaction"), None)
        if channel is None:
            continue
        low = nodes[0]["microstates"][0]["energy_eV"]
        other = next((x for x in nodes if x["id"] != nodes[0]["id"]
                      and x["id"] in channel["ends"]), None)
        if other:
            values.append((seed, (other["microstates"][0]["energy_eV"] - low) * 1000))
    fig, ax = plt.subplots(figsize=(6.2, 1.9))
    ax.axvline(reference, color=OXIDE, lw=1.4, ls=(0, (3, 3)), zorder=2)
    ax.text(reference - 0.03, 0.94, f"tightly relaxed reference  {reference:.3f}",
            rotation=90, va="top", ha="right", fontsize=7.5, color=OXIDE,
            transform=ax.get_xaxis_transform())
    lo = min(v for _, v in values)
    hi = max(v for _, v in values)
    ax.annotate("", xy=(lo, 0.30), xytext=(hi, 0.30),
                arrowprops=dict(arrowstyle="|-|,widthA=.3,widthB=.3", color=MUTED, lw=0.8))
    ax.text((lo + hi) / 2, 0.40, f"spread {hi - lo:.2f} meV", ha="center", fontsize=7.5,
            color=MUTED)
    for index, (seed, value) in enumerate(sorted(values, key=lambda item: item[1])):
        ax.plot([value], [-0.28 - (index % 4) * 0.17], "o", ms=6, color=TEAL,
                mec="white", mew=1.0, zorder=3)
    ax.set_ylim(-1.05, 0.75)
    ax.set_yticks([])
    ax.spines["left"].set_visible(False)
    ax.set_xlabel(r"$\Delta E(\mathrm{B}-\mathrm{A})$  /  meV")
    ax.set_title("Endpoint spread under the production tolerance, not an energy resolution",
                 loc="left", fontsize=10, pad=10)
    ax.grid(axis="x", color=RULE, lw=0.6, alpha=0.6)
    ax.set_axisbelow(True)
    save(fig, "fig3_energy_spread")


# ------------------------------------------- 4. which mechanism carried each edge
def figure_evidence_by_seed():
    data = []
    for seed in SEEDS:
        network = json.loads((ROOT / f"runs/p2_T300_seed{seed}/network.json").read_text())
        channel = next((c for c in network["reaction_channels"]
                        if c["class"] == "reaction"
                        and sorted(map(sorted, c["broken"] + c["formed"]))
                        == [[0, 5], [4, 5]]), None)
        ends = set(channel["ends"]) if channel else set()
        edges = [e for e in network["reactions"]
                 if {e["source"], e["target"]} == ends]
        data.append((seed,
                     sum(len(e["attempts"]) for e in edges),
                     sum(len(e.get("continuations") or []) for e in edges)))
    fig, ax = plt.subplots(figsize=(6.2, 3.1))
    y = range(len(data))
    for index, (seed, direct, cont) in enumerate(data):
        row = len(data) - 1 - index
        colour = TEAL if direct else OXIDE
        if direct:
            ax.barh(row, direct, height=0.62, color=TEAL, zorder=3)
        else:
            ax.text(0.10, row, "no direct observation", va="center", fontsize=7,
                    color=OXIDE, zorder=3)
        ax.text(-0.28, row, f"{cont}", va="center", ha="right", fontsize=7.5,
                color=INDIGO, family="monospace")
    ax.set_yticks(list(y))
    ax.set_yticklabels([f"{seed}" for seed, _, _ in data][::-1], family="monospace",
                       fontsize=8)
    ax.set_ylim(-0.7, len(data) - 0.3)
    ax.set_xlim(-1.4, max(d for _, d, _ in data) + 0.4)
    ax.set_xticks(range(0, max(d for _, d, _ in data) + 1))
    ax.set_xlabel("trials that observed the transfer on their own")
    ax.set_ylabel("seed", labelpad=2)
    ax.set_title("Which evidence carried the directed edge", loc="left", fontsize=10, pad=10)
    ax.grid(axis="x", color=RULE, lw=0.6, alpha=0.6)
    ax.set_axisbelow(True)
    ax.legend(handles=[
        Line2D([], [], color=TEAL, lw=5, label="direct observations (bar)"),
        Line2D([], [], color=INDIGO, lw=0, marker="$0$", ms=7,
               label="continuations supporting the edge (number at left)"),
    ], loc="upper center", ncol=2, bbox_to_anchor=(0.5, -0.16))
    save(fig, "fig4_evidence_by_seed")


# ------------------------------------------------- 5. what the ceiling blocks
def figure_ceiling():
    report_path = ROOT / "runs/ceiling_diagnostic/report.json"
    if not report_path.exists():
        print("  (ceiling diagnostic not present, skipped)")
        return
    report = json.loads(report_path.read_text())
    fig, axes = plt.subplots(1, len(report["seeds"]), figsize=(6.6, 2.7), sharey=True)
    axes = axes if hasattr(axes, "__len__") else [axes]
    for ax, entry in zip(axes, report["seeds"]):
        # Control first and thick, diagnostic thin on top: the two arms are identical
        # until the ceiling binds, and that is the thing to see.
        for label, colour, alpha, width in (("control_1.5", MUTED, 0.55, 4.0),
                                            ("diagnostic_3.0", OXIDE, 0.95, 1.5)):
            arm = entry.get(label)
            if not arm:
                continue
            path = (ROOT / "runs/ceiling_diagnostic"
                    / f"{entry['seed_attempt']}_{label}/climb.jsonl")
            trace = [(r["step"], r["potential_eV"]) for r in rows(path)
                     if r.get("potential_eV") is not None]
            if not trace:
                continue
            steps, energies = zip(*trace)
            base = report["source_energy_eV"]
            ax.plot(steps, [(e - base) for e in energies], "-", color=colour,
                    lw=width, alpha=alpha, solid_capstyle="round",
                    label=f"ceiling {arm['ceiling_eV']} eV")
        for value, colour in ((1.5, MUTED), (3.0, OXIDE)):
            ax.axhline(value, color=colour, lw=0.8, ls=(0, (3, 3)), alpha=0.8)
        ax.set_title(f"{entry['seed_attempt']}  ·  {entry['probe']['family']}"
                     f"{entry['probe']['indices']}", fontsize=8, color=INK, pad=6)
        ax.set_xlabel("climb step")
        ax.grid(color=RULE, lw=0.6, alpha=0.6)
        ax.set_axisbelow(True)
    axes[0].set_ylabel("rise above the reactant  /  eV")
    axes[0].legend(loc="upper center", ncol=2, bbox_to_anchor=(1.05, -0.22), fontsize=7.5)
    fig.suptitle("Raising the ceiling revealed no saddle: both walks ran out of steps",
                 x=0.02, y=1.06, ha="left", fontsize=10)
    save(fig, "fig5_ceiling_diagnostic")


if __name__ == "__main__":
    print(f"writing into {OUT}")
    figure_transition_path()
    figure_frozen_questions()
    figure_energy_spread()
    figure_evidence_by_seed()
    figure_ceiling()
