"""Run directory -> DiscoveryTrace list (design section 6.1).

Selection of which trials become traces follows the run's own
classification, which stays authoritative on the PRRS side. Every view filters
the same checked records (result.json cross-checked against attempts.jsonl),
so a ledger/result conflict fails ingestion instead of silently skipping.
"""

from __future__ import annotations

from ..adapters.prrs_run import PrrsRun
from ..config import ParsConfig
from ..schema import DiscoveryTrace


def _accepted(record: dict, config: ParsConfig) -> bool:
    # only completed trials are confirmed A->B paths: ts_candidate stopped on a
    # saddle, unconfirmed did not hold its new topology through the free tail
    return (
        record.get("status") == "completed"
        and record.get("classification") in config.include_classifications
    )


def ingest_run(run_dir, config: ParsConfig) -> list[DiscoveryTrace]:
    """Load one trace per completed attempt whose classification passes the
    configured filter. Order follows attempts.jsonl (chronological).
    ``run_dir`` may also be an already opened PrrsRun."""
    run = run_dir if isinstance(run_dir, PrrsRun) else PrrsRun(run_dir)
    return [
        run.load_trace(record["attempt_id"])
        for record in run.checked_attempts()
        if _accepted(record, config)
    ]


def skipped_attempts(run: PrrsRun, config: ParsConfig) -> dict[str, int]:
    """Attempts ingestion left out, counted by ``status/classification``."""
    counts: dict[str, int] = {}
    for record in run.checked_attempts():
        if not _accepted(record, config):
            key = f"{record.get('status')}/{record.get('classification')}"
            counts[key] = counts.get(key, 0) + 1
    return dict(sorted(counts.items()))


def channel_members_not_ingested(run: PrrsRun, config: ParsConfig) -> dict[str, list[str]]:
    """Per channel, the observing attempts PRRS counted but PARS did not take
    as paths (e.g. unconfirmed crossings). Empty channels are omitted."""
    accepted = {r["attempt_id"] for r in run.checked_attempts() if _accepted(r, config)}
    missing: dict[str, list[str]] = {}
    for channel in run.network["reaction_channels"]:
        left_out = sorted(set(channel.get("observed_by", ())) - accepted)
        if left_out:
            missing[channel["id"]] = left_out
    return dict(sorted(missing.items()))
