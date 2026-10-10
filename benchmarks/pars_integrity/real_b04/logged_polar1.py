"""polar1 with a record of the electronic state of every structure it evaluates."""
import atexit, collections, json, os
from benchmarks.b04_sn2_chloride.calculator import polar1

SEEN = collections.Counter()


def make():
    calc = polar1()
    original = calc.calculate

    def calculate(atoms=None, *args, **kwargs):
        SEEN[str((atoms.info.get("charge"), atoms.info.get("spin")))] += 1
        return original(atoms, *args, **kwargs)

    calc.calculate = calculate
    atexit.register(_write_seen)
    return calc


def _write_seen():
    # Add to the counts already on disk rather than replace them: a resumed refine that finds
    # every representative done evaluates nothing, and on 2026-10-08 two such reruns overwrote
    # the first batch's counts with {}.
    path = os.environ["SEEN_OUT"]
    total = collections.Counter()
    if os.path.exists(path):
        total.update(json.loads(open(path).read()))
    total.update(SEEN)
    open(path, "w").write(json.dumps(total))
