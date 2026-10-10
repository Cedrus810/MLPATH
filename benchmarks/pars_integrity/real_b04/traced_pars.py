"""`python -m pars ...` with every thread's stack dumped to stderr every TRACE_EVERY seconds.

Diagnostic for the probe3 t4 refinement that ran ~67 min without output and then vanished
(2026-10-05 on kasuga01, 2026-10-08 on kasuga02): the periodic dumps show where it spends
its time, and PYTHONFAULTHANDLER covers a hard crash. A SIGKILL leaves neither; then the
exit code in gpu_batch.sh's RUN_EXIT line is the only record.
"""
import faulthandler, os, runpy, sys

faulthandler.dump_traceback_later(int(os.environ["TRACE_EVERY"]), repeat=True)
sys.argv[0] = "pars"
runpy.run_module("pars", run_name="__main__", alter_sys=True)
