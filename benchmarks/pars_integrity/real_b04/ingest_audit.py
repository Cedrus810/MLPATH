"""Import every new-format PRRS run under runs/ with the strict 2026-10-05 checks."""
import json, sys
from pathlib import Path
from pars.adapters.prrs_run import PrrsRun
from pars.config import ParsConfig
from pars.paths.ingest import ingest_run, skipped_attempts

rows = []
for network in sorted(Path("runs").rglob("network.json")):
    run_dir = network.parent
    row = {"run": str(run_dir)}
    try:
        payload = json.loads(network.read_text())
        if "reaction_channels" not in payload:
            row["result"] = "old_format_skip"; rows.append(row); continue
        if not (run_dir / "attempts.jsonl").is_file():
            row["result"] = "no_ledger_skip"; rows.append(row); continue
        run = PrrsRun(run_dir)
        row["state"] = run.electronic_state
        row["uncertainty_status"] = sorted({r.get("uncertainty_status") for r in run.attempts()}, key=str)
        traces = ingest_run(run, ParsConfig())
        row["traces"] = len(traces)
        row["skipped"] = skipped_attempts(run, ParsConfig())
        row["result"] = "ok"
    except Exception as error:
        row["result"] = f"REJECTED {type(error).__name__}: {error}"
    rows.append(row)
json.dump(rows, open(sys.argv[1], "w"), indent=1, default=str)
for row in rows:
    print(row["run"], "|", row["result"][:160], "|", row.get("state"), "| traces", row.get("traces"), "| unc", row.get("uncertainty_status"))
