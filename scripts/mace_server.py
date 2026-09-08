"""Persistent MACE energy/gradient server for ORCA's ExtOpt interface.

Why this exists: ORCA calls the ExtOpt program once per gradient, and a numerical
Hessian on N atoms is 6N of those. Loading a MACE checkpoint and initialising CUDA
costs about 5 s, while evaluating 9 atoms costs milliseconds -- measured 213 s for
roughly 40 calls, essentially all of it reloading. The model is loaded once here and
stays resident.

Launched the way opi.external_methods.server.OpiServer launches a server:

    python scripts/mace_server.py -b 127.0.0.1:8888

The model is named by environment, the same variables scripts/mace_engrad.py uses,
so a server and a direct call cannot disagree about which checkpoint is in play:

    PRRS_ENGRAD_MODEL, PRRS_ENGRAD_TYPE, PRRS_ENGRAD_DEVICE

Wire format is a pickled dict over multiprocessing.connection, matching OPI's Client:

    {"type": "engrad", "numbers": [...], "positions": [[x,y,z],...],
     "charge": int, "spin": int, "gradient": bool}
      -> {"energy_eV": float, "forces_eV_A": [[fx,fy,fz],...] | None}
    {"type": "ping"}      -> {"model": <path>, "device": <device>, "calls": int}
    {"type": "shutdown"}  -> {"ok": True}

The calculator is NOT sent over the wire. A MACECalculator holds CUDA tensors and
pickling it is a fragile path; naming the checkpoint is not.
"""

import os
import sys
from multiprocessing.connection import Listener


def parse_bind(argv):
    for i, a in enumerate(argv):
        if a == "-b" and i + 1 < len(argv):
            host, _, port = argv[i + 1].rpartition(":")
            return host or "127.0.0.1", int(port)
    return "127.0.0.1", 8888


def main(argv):
    host, port = parse_bind(argv)
    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    sys.path.insert(0, os.path.join(root, "src"))

    model = os.environ["PRRS_ENGRAD_MODEL"]
    device = os.environ.get("PRRS_ENGRAD_DEVICE", "cuda")
    from prrs import torch_guard

    torch_guard.install(trusted=[model])
    from ase import Atoms
    from mace.calculators import MACECalculator

    calc = MACECalculator(
        model_paths=model,
        model_type=os.environ.get("PRRS_ENGRAD_TYPE", "MACE"),
        device=device,
        default_dtype="float64",
    )
    calls = 0
    listener = Listener((host, port))
    print(f"mace_server ready on {host}:{port}  model={model}  device={device}", flush=True)

    while True:
        with listener.accept() as conn:
            try:
                message = conn.recv()
            except EOFError:
                continue
            kind = message.get("type")
            if kind == "shutdown":
                conn.send({"ok": True})
                break
            if kind == "ping":
                conn.send({"model": model, "device": device, "calls": calls})
                continue
            if kind != "engrad":
                # setup_calculator is OPI's own message; we load from the path instead,
                # so acknowledge it rather than failing a caller that sends it.
                conn.send({"ok": True, "note": "server loads its model from the environment"})
                continue
            atoms = Atoms(numbers=message["numbers"], positions=message["positions"])
            atoms.info["charge"] = message.get("charge", 0)
            atoms.info["spin"] = message.get("spin", 1)
            atoms.calc = calc
            calc.results.clear()
            energy = float(atoms.get_potential_energy())
            forces = atoms.get_forces().tolist() if message.get("gradient", True) else None
            calls += 1
            conn.send({"energy_eV": energy, "forces_eV_A": forces})
    listener.close()
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
