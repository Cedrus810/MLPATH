import json
import os
import subprocess
import sys
from collections import Counter
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]


def drive(run_dir, counter, crash_after):
    path = os.pathsep.join(str(ROOT / p) for p in ("src", "../src", "tests"))
    return subprocess.run(
        [
            sys.executable,
            str(ROOT / "tests" / "d2_driver.py"),
            str(run_dir),
            str(counter),
            str(crash_after),
        ],
        env={**os.environ, "PYTHONPATH": path},
        capture_output=True,
        text=True,
        timeout=900,
    )


def events(run_dir):
    text = (run_dir / "journal.jsonl").read_text(encoding="utf-8")
    return [json.loads(line) for line in text.splitlines() if line.strip()]


def report(run_dir):
    return json.loads((run_dir / "report.json").read_text(encoding="utf-8"))


def paid(rep):
    return rep["accounts"]["search"]["spent"] + rep["accounts"]["certify"]["spent"]


@pytest.fixture(scope="module")
def reference(tmp_path_factory):
    base = tmp_path_factory.mktemp("ref")
    done = drive(base / "run", base / "count", -1)
    assert done.returncode == 0, done.stderr
    return base / "run", (base / "count").stat().st_size


def test_uninterrupted_run_external_count_matches_every_ledger(reference):
    run_dir, external = reference
    rep = report(run_dir)
    assert external == rep["ledger_total"] == paid(rep)
    assert rep["unreconciled"] == []


def crash_point(ref_dir, phase):
    if phase == "L1":
        return 5  # scans of the first wave are in flight
    first = next(
        e for e in events(ref_dir) if e["event"] == "settle" and e["task"].startswith("L4|")
    )
    return first["cumulative_spent"] - first["spent"] // 2  # inside the first certification


@pytest.mark.parametrize("phase", ["L1", "L4"])
def test_kill_and_resume_pays_once_and_finds_the_same_saddles(reference, tmp_path, phase):
    ref_dir, ref_external = reference
    run_dir, counter = tmp_path / "run", tmp_path / "count"
    first = drive(run_dir, counter, crash_point(ref_dir, phase))
    assert first.returncode == 137, first.stderr
    cut = events(run_dir)
    settled = {(e["task"], e["attempt"]) for e in cut if e["event"] == "settle"}
    open_kinds = {
        e["kind"]
        for e in cut
        if e["event"] == "dispatch" and (e["task"], e["attempt"]) not in settled
    }
    assert phase in open_kinds

    second = drive(run_dir, counter, -1)
    assert second.returncode == 0, second.stderr
    log, rep, ref = events(run_dir), report(run_dir), report(ref_dir)

    # same science as the uninterrupted run
    def energies(r):
        return sorted(round(c["saddle_energy_eV"], 8) for c in r["c2"])

    assert energies(rep) == energies(ref)
    assert rep["release_status"] == ref["release_status"]
    assert rep["certify_status"] == ref["certify_status"]

    # every task settles once; an interrupted attempt is retried under a larger number
    assert max(Counter(e["task"] for e in log if e["event"] == "settle").values()) == 1
    reconciles = [e for e in log if e["event"] == "reconcile"]
    assert reconciles
    for r in reconciles:
        assert any(
            e["event"] == "settle" and e["task"] == r["task"] and e["attempt"] > r["attempt"]
            for e in log
        )

    # cost: nothing refunded, nothing paid twice
    wasted = sum(r["spent"] for r in reconciles)
    assert paid(rep) == paid(ref) + wasted
    external = counter.stat().st_size
    assert external <= paid(rep)
    assert external - ref_external <= wasted
