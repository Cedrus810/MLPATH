"""python d2_driver.py RUN_DIR COUNTER CRASH_AFTER [hf|v14]: run or resume the TwoSaddle D2
bench (``hf`` freezes rule D2-HF-v1, ``v14`` rule D2-HF-v1.4 with parallel waves). Each model evaluation appends one byte to COUNTER first; evaluation number
CRASH_AFTER of this process calls os._exit(137) (no finally blocks, like SIGKILL).
CRASH_AFTER < 0 never crashes."""

import sys
from pathlib import Path

from d2_bench import HF, Tick, backends, init
from explore.d2 import run


def main():
    run_dir, counter, crash_after = Path(sys.argv[1]), Path(sys.argv[2]), int(sys.argv[3])
    mode = sys.argv[4] if len(sys.argv) > 4 else ""
    if not (run_dir / "manifest.json").exists():
        if mode == "v14":
            # D2-HF-v1.4, every same-graph ridge sampled, parallel waves of up to 4
            from explore.d2 import HSBracket, StartCluster

            init(
                run_dir,
                hf=HF,
                start_cluster=StartCluster(1e4, 1),
                same_graph_ridge_last=True,
                hs_bracket=HSBracket(2e-3, 0.025),
                same_graph_sample=True,
                hf_parallel=True,
                audit_fraction=1.0,
                max_active=2,
                wave_size=4,
            )
        else:
            init(run_dir, hf=HF if mode == "hf" else None)
    run(run_dir, backends(Tick(counter, crash_after)))


if __name__ == "__main__":
    main()
