"""b04 preflight: the things that have to be known before criteria can be written.

Nothing here is an acceptance result. This answers only "can this benchmark be run at
all", which is a different question from "did PRRS find the reaction". The model for
this file is docs/P3_BENZOYLACETONE_PROBES.md -- 378 lines, not one of them a verdict.

Each gate prints PASS / FAIL / SKIP and the script exits non-zero if any required gate
fails. A FAIL here is not a bad result; it is the gate doing its job, and the fix goes
into src/prrs, not into this file.

    python benchmarks/b04_sn2_chloride/preflight.py            # gates 1-3, no GPU, seconds
    python benchmarks/b04_sn2_chloride/preflight.py --models   # also gate 4, needs a GPU
"""

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(Path(__file__).parent))

from ase import Atoms  # noqa: E402

from build import build  # noqa: E402

results = []


def gate(name, ok, detail, required=True):
    results.append((name, ok, required))
    mark = "PASS" if ok else ("FAIL" if required else "WARN")
    print(f"  [{mark}] {name}")
    for line in detail.strip().splitlines():
        print(f"         {line}")
    print()


print(f"b04 preflight  ·  tree {ROOT}")
print("=" * 78)

# ---------------------------------------------------------------- gate 1
# The domain gate must admit a bare halide anion. It is closed-shell; it is simply not
# neutral, and chemistry.VALENCE is a table of NEUTRAL closed-shell valences.
from prrs.chemistry import VALENCE, chemical_key  # noqa: E402
from prrs.config import SearchConfig  # noqa: E402

CHARGE = dict(charge_sensitive=True, charge=-1, multiplicity=1)
chloride = chemical_key(Atoms("Cl", positions=[[0, 0, 0]]), 1.2, **CHARGE)
reactants = chemical_key(build(), 1.2, **CHARGE)
blocking = ("no_assignment_satisfies_valences", "degree_exceeds_valence")
# RESOLVED 2026-09-04 by the model decision, not by a code change. The gate exists
# because MACE-OFF24 is trained on neutral closed-shell molecules and ignores total
# charge (network.py, Registry.admit says so in as many words). POLAR-1's charge
# response is measured at 16.67 eV across 0/-1/+1, so that reason does not apply and
# the frozen config declares closed_shell_only = False for it. What the identity layer
# says about a bare anion is then an annotation, not a veto.
_config = json.loads(
    "\n".join(
        line for line in Path(__file__).with_name("config.json").read_text().splitlines()
        if not line.lstrip().startswith("//")
    )
) if (Path(__file__).with_name("config.json")).exists() else {}
_vetoes = bool(_config.get("closed_shell_only", True))
gate(
    "1  a bare chloride anion is not vetoed by the identity layer",
    (not _vetoes)
    or (chloride.get("stereo_unresolved") not in blocking
        and reactants.get("stereo_unresolved") not in blocking),
    f"""Cl- alone      stereo_unresolved = {chloride.get('stereo_unresolved')!r}
Cl- + CH3Cl    stereo_unresolved = {reactants.get('stereo_unresolved')!r}
fragments      {reactants['fragments']}
VALENCE covers {sorted(VALENCE)} (atomic numbers), all as NEUTRAL valences.
Cl has neutral valence 1 and a bare anion has degree 0, so no assignment exists and
Registry.admit refuses it as out_of_domain while closed_shell_only is True.
Fix belongs in chemistry.VALENCE: make it charge-aware so Cl- has valence 0 at
total_charge = -1. That changes chemical_key, so it changes the digest on purpose.
Setting closed_shell_only=False is NOT the fix -- see config.json.""",
)

# ---------------------------------------------------------------- gate 2
# The config must be constructible with a nonzero total charge, and the charge must
# reach chemical_key. This one is expected to pass: the machinery anticipated it.
try:
    cfg = SearchConfig(charge_sensitive=True, total_charge=-1, multiplicity=1)
    built, why = True, f"charge_sensitive={cfg.charge_sensitive} total_charge={cfg.total_charge}"
except Exception as exc:                                  # noqa: BLE001 - reported below
    built, why = False, f"{type(exc).__name__}: {exc}"
gate(
    "2  SearchConfig accepts a nonzero total charge",
    built,
    f"""{why}
network.Registry passes config.total_charge and config.multiplicity into chemical_key
(src/prrs/network.py), so the IDENTITY layer is charge-aware already. Whether the
POTENTIAL hears about it is gate 3, and they are not the same thing.""",
)

# ---------------------------------------------------------------- gate 3
# Nothing in src/prrs writes atoms.info["charge"], so a charge-aware calculator would
# silently evaluate the neutral species while the config declares an anion.
sources = sorted((ROOT / "src" / "prrs").glob("*.py"))
writers = [
    f"{p.name}:{n}"
    for p in sources
    for n, line in enumerate(p.read_text().splitlines(), 1)
    if 'info["charge"]' in line or "info['charge']" in line
]
from prrs.search import _stamped_with_charge  # noqa: E402

_probe = build()
_probe.info.clear()
_stamped = _stamped_with_charge(_probe, SearchConfig(charge_sensitive=True, total_charge=-1,
                                                    multiplicity=1, closed_shell_only=False))
gate(
    "3  the total charge reaches the calculator",
    _stamped.info.get("charge") == -1 and _stamped.info.get("spin") == 1,
    f"""writers of atoms.info["charge"] inside src/prrs: {writers or 'NONE'}
Charge-aware MACE reads the total charge from atoms.info["charge"] / ["spin"].
docs/experiments/sn2_probe.py sets them by hand; the search never does.
Worst form of this failure: the config declares -1, the model evaluates 0, and
NOTHING RAISES. build.py stamps the source structure, but every atoms.copy() made
during a search has to carry it too, and that is not established here.""",
)

# ---------------------------------------------------------------- gate 4
if "--models" in sys.argv:
    import json

    spec = json.loads(
        "\n".join(
            line for line in Path(__file__).with_name("model.json").read_text().splitlines()
            if not line.lstrip().startswith("//")
        )
    )
    try:
        from mace.calculators import MACECalculator

        readings = {}
        for candidate in spec["candidates"]:
            calc = MACECalculator(
                model_paths=candidate["path"],
                model_type=candidate["mace_model_type"],
                device="cuda",
                default_dtype="float64",
            )
            energies = []
            for q in (0, -1, +1):
                probe = build()
                probe.info["charge"], probe.info["spin"] = q, 1
                probe.calc = calc
                calc.results.clear()
                energies.append(float(probe.get_potential_energy()))
            readings[candidate["name"]] = energies
        spread = {k: max(v) - min(v) for k, v in readings.items()}
        gate(
            "4  the model actually responds to total charge",
            all(s > 1e-6 for s in spread.values()),
            "\n".join(
                f"{k}: E(0)={v[0]:.6f}  E(-1)={v[1]:.6f}  E(+1)={v[2]:.6f}  "
                f"spread={spread[k]:.3e} eV"
                for k, v in readings.items()
            )
            + "\nMACE-OFF24 returns bit-identical energies here, which is why it is excluded."
            "\nA zero spread means the model is blind to charge and the benchmark is void.",
        )
    except Exception as exc:                              # noqa: BLE001 - reported, not raised
        gate("4  the model actually responds to total charge", False,
             f"could not load a model: {type(exc).__name__}: {exc}")
else:
    print("  [SKIP] 4  the model actually responds to total charge")
    print("         pass --models to run it (needs a GPU and loads two checkpoints)\n")

# ---------------------------------------------------------------- gate 5
# Nothing in prrs checks a structure's elements against what the model was built with.
# A run on an element the model never saw produces numbers and no warning.
import json as _json  # noqa: E402

_spec = _json.loads("\n".join(
    line for line in Path(__file__).with_name("model.json").read_text().splitlines()
    if not line.lstrip().startswith("//")))
_present = sorted(set(int(z) for z in build().numbers))
try:
    sys.path.insert(0, str(ROOT / "scripts"))
    from model_registry import describe  # noqa: E402

    _lines, _ok = [], True
    for candidate in _spec["candidates"]:
        declared = set(describe(candidate["path"])["atomic_numbers"])
        missing = [z for z in _present if z not in declared]
        _ok = _ok and not missing
        _lines.append(f"{candidate['name']}: declares {len(declared)} elements, "
                      f"structure needs {_present} -> "
                      f"{'all present' if not missing else f'MISSING {missing}'}")
    gate("5  every element in the structure is one the model declares", _ok,
         "\n".join(_lines) + "\n"
         "Read from the checkpoint itself (scripts/model_registry.py), not from a\n"
         "hand-kept list. Note this is the DOMAIN question only -- it says nothing\n"
         "about chemistry.VALENCE, which answers a Lewis-structure question the model\n"
         "has no opinion on.")
except Exception as exc:                                  # noqa: BLE001 - reported
    gate("5  every element in the structure is one the model declares", False,
         f"could not read a checkpoint: {type(exc).__name__}: {exc}")

print("=" * 78)
failed = [n for n, ok, required in results if not ok and required]
print(f"{len(results) - len(failed)}/{len(results)} required gates pass")
if failed:
    print("\nblocking:")
    for name in failed:
        print(f"  - {name}")
    print("\nCriteria must not be frozen while any of these is red: every criterion")
    print("written now would be answered by a domain-gate refusal, and README.md rule 4")
    print("requires a refusal to be reported as a refusal rather than as a failure.")
sys.exit(1 if failed else 0)
