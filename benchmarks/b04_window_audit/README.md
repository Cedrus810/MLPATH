# b04 response-horizon audit, 2026-10-04

This is a real POLAR-1 CPU replay of **scan 0** of the legacy `b04_frozen2`
run, not a whole-network acceptance run. Source/PES historical identity is
limited by the old manifest; the current checkpoint contents and numerical
environment are frozen in [the report](results/2026-10-04/report.json).

The seven published sample endpoints at the original 200 steps (50 fs)
reproduced within absolute 1e-5 Angstrom before any alternate horizon ran.
Then exactly the same source/probe/seed/amplitude set was run at 100 steps
(25 fs) and 400 steps (100 fs), with the other configuration unchanged.
Float64 CPU, two Torch threads, the explicitly hashed local
`/home/ruigengji/MLP/mace/MACE-POLAR-1-M.model`, 21 trials, approximately
2429.57 seconds. The source manifest/network/attempts hashes were unchanged.

| Horizon | Returned-to-departed window | Unresolved samples |
| --- | --- | --- |
| 25 fs | none certified | 0.80 (unconfirmed transient) |
| 50 fs | [0.75, 0.80] | none |
| 100 fs | [0.30, 0.60] | 0.75 (unconfirmed transient) |

At 100 fs, amplitudes 0.80 and 1.00 return to the source, whereas they
reached the other microstate at 50 fs. Overall `not_certified`, with
`label_changes=true`. Transient endpoints are not promoted to confirmed
departures: this stricter rule explains why the 25 fs certified-window
summary differs from the old geometry-only exploratory summary.
No static basin-boundary, infinite-time stability or whole-network claim.

Replaying requires the original source corpus and the local checkpoint.
Do not recreate the source run or overwrite its results. From the repository
root, create JSON factory kwargs with `device="cpu"`,
`default_dtype="float64"`, and metadata with `local_file_sha256` mapping the
absolute checkpoint path to its current content SHA256. Use a new directory:

```sh
PYTHONPATH=src:. /home/ruigengji/miniforge3/envs/openmm_dev/bin/python -m prrs audit-windows \
  --run runs/b04_frozen2 --output runs/my_b04_window_audit \
  --calculator benchmarks.b04_sn2_chloride.calculator:polar1 \
  --calculator-kwargs /path/to/cpu_kwargs.json \
  --calculator-metadata /path/to/checkpoint_metadata.json \
  --response-steps 100 400 --scans 0
```

The command returns 2 for this non-certified result. To match the recorded
thread count, call `replay_window_audit` from a Python process after
`torch.set_num_threads(2)`. All model loading remains explicit and local.
The complete repair/evidence boundary is in
[MAJOR_FIXES_2026-10-04](../../docs/MAJOR_FIXES_2026-10-04.md).
