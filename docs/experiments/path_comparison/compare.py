"""MLP 下降路径 vs QM IRC：同一反应、两条构造、两张势能面。

跑法： python docs/experiments/path_comparison/compare.py
"""

import json
from pathlib import Path

import numpy as np
from ase.io import read

ROOT = Path(__file__).resolve().parents[3]
MLP_XYZ = ROOT / "runs/b04_probe2_openceiling/paths/ts0000/descent.extxyz"
MLP_LOG = ROOT / "runs/b04_probe2_openceiling/paths/ts0000/descent.jsonl"
IRC_XYZ = ROOT / "docs/experiments/b04_qm/irc_IRC_Full_trj.xyz"


def bond_pair(atoms):
    """(短 C-Cl, 长 C-Cl)。有序，消掉两个氯的标号。"""
    s = atoms.get_chemical_symbols()
    c = s.index("C")
    cl = [i for i, x in enumerate(s) if x == "Cl"]
    return tuple(sorted((atoms.get_distance(c, cl[0]), atoms.get_distance(c, cl[1]))))


def arc(frames):
    p = np.array([a.positions for a in frames])
    step = np.linalg.norm(np.diff(p, axis=0), axis=(1, 2))
    return step, np.concatenate([[0.0], np.cumsum(step)])


def report():
    mlp, irc = read(MLP_XYZ, index=":"), read(IRC_XYZ, index=":")
    rows = [json.loads(line) for line in MLP_LOG.read_text().splitlines()]
    ms, ma = arc(mlp)
    qs, qa = arc(irc)

    print(f"MLP descent  {len(mlp):3d} 帧  弧长 {ma[-1]:.3f} A  每帧中位 {np.median(ms):.4f} A")
    print(f"QM IRC       {len(irc):3d} 帧  弧长 {qa[-1]:.3f} A  每帧中位 {np.median(qs):.4f} A")

    print("\nMLP 各帧段吃掉多少弧长")
    for a, b in ((1, 10), (10, 20), (20, 30), (30, 40), (40, 50), (50, len(ms))):
        seg = ms[a:b]
        if not len(seg):
            continue
        print(
            f"  {a:2d}-{b:<2d} {100 * seg.sum() / ma[-1]:5.1f}%  每帧中位 {np.median(seg):.4f} A"
        )

    # 能量停止变化之后还有多少帧
    for sign in (1, -1):
        e = np.array([r.get("potential_eV") or np.nan for r in rows if r.get("sign") == sign])
        if len(e) < 2:
            continue
        moving = np.where(np.abs(np.diff(e)) > 1e-4)[0]
        last = moving[-1] + 1 if len(moving) else 0
        print(
            f"  sign={sign:+d}: 能量最后一次变化在第 {last} 帧，其后 {len(e) - last - 1} 帧收敛尾巴"
        )

    print("\n同一长键处的短键（两条路径的投影对照）")
    tm = np.array([bond_pair(a) for a in mlp])
    tq = np.array([bond_pair(a) for a in irc])
    print(f"  {'长 C-Cl':>8s} {'IRC 短键':>9s} {'MLP 短键':>9s} {'差':>8s}")
    for target in (2.40, 2.50, 2.60, 2.80, 3.00):

        def near(t):
            i = int(np.argmin(abs(t[:, 1] - target)))
            return t[i, 0] if abs(t[i, 1] - target) < 0.25 else None

        a, b = near(tq), near(tm)
        if a and b:
            print(f"  {target:8.2f} {a:9.3f} {b:9.3f} {b - a:+8.3f}")


if __name__ == "__main__":
    report()
