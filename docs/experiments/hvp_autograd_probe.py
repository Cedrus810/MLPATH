"""Does exact double-backward Hv work on MACE-OFF24_medium, and does it agree with FD?"""
import time, numpy as np, torch
from ase.build import molecule
from mace.calculators import MACECalculator

MODEL = "/home/kasuga/.cache/mace/MACE-OFF24_medium.model"
DEV = "cuda"
rng = np.random.default_rng(0)

atoms = molecule("CH3CH2OH")
atoms.rattle(0.05, seed=1)          # generic point, not a stationary one
calc = MACECalculator(model_paths=MODEL, device=DEV, default_dtype="float64")
atoms.calc = calc
E0 = atoms.get_potential_energy()
F0 = atoms.get_forces().copy()
n = len(atoms)
print(f"atoms={n}  E0={E0:.9f} eV  fmax={np.abs(F0).max():.4f} eV/A")

model = calc.models[0]

def graph_gradient():
    """gradient of E wrt positions, with create_graph=True so it can be differentiated again"""
    batch = calc._atoms_to_batch(atoms)
    bd = batch.clone().to_dict()
    pos = bd["positions"]
    pos.requires_grad_(True)
    out = model(bd, compute_force=False, compute_stress=False, training=True)
    energy = out["energy"].sum()
    grad = torch.autograd.grad(energy, pos, create_graph=True)[0]
    return pos, grad, float(energy.detach())

try:
    pos, grad, E_graph = graph_gradient()
except Exception as exc:
    print("DOUBLE BACKWARD SETUP FAILED:", type(exc).__name__, exc)
    raise SystemExit(1)

print(f"energy via direct model call: {E_graph:.9f}  dE vs calculator = {E_graph - E0:.3e}")
print(f"force  via direct model call: max|dF| vs calculator = "
      f"{np.abs(-grad.detach().cpu().numpy() - F0).max():.3e} eV/A")

def hvp_exact(v, retain=True):
    return torch.autograd.grad((grad * v).sum(), pos, retain_graph=retain)[0]

def hvp_fd(v_np, eps):
    ref = atoms.positions.copy()
    out = []
    for sign in (1, -1):
        atoms.set_positions(ref + sign * eps * v_np)
        out.append(-atoms.get_forces().copy())
    atoms.set_positions(ref)
    return (out[0] - out[1]) / (2 * eps)

# three probe directions: random, and the two largest-force-component unit vectors
vs = []
for _ in range(2):
    w = rng.normal(size=(n, 3)); vs.append(w / np.linalg.norm(w))
w = np.zeros((n, 3)); w[int(np.abs(F0).sum(1).argmax()), 0] = 1.0; vs.append(w)

torch.cuda.reset_peak_memory_stats()
print("\n--- exact double-backward Hv vs central-difference Hv (eV/A^2) ---")
print(f"{'v':>3} {'|Hv_exact|inf':>14} {'eps=1e-3':>12} {'eps=1e-2':>12} {'eps=5e-2':>12}   (max abs dev)")
exact_store = []
for i, v_np in enumerate(vs):
    v_t = torch.as_tensor(v_np, dtype=pos.dtype, device=pos.device)
    try:
        hv = hvp_exact(v_t).detach().cpu().numpy()
    except Exception as exc:
        print("DOUBLE BACKWARD FAILED:", type(exc).__name__, exc); raise SystemExit(1)
    exact_store.append((v_np, hv))
    devs = [np.abs(hvp_fd(v_np, e) - hv).max() for e in (1e-3, 1e-2, 5e-2)]
    print(f"{i:>3} {np.abs(hv).max():>14.6f} " + " ".join(f"{d:>12.3e}" for d in devs))

# symmetry: v^T H w must equal w^T H v
(v0, hv0), (v1, hv1) = exact_store[0], exact_store[1]
print(f"\nsymmetry  v0.H v1 = {np.sum(v0*hv1):.12f}   v1.H v0 = {np.sum(v1*hv0):.12f}"
      f"   diff = {abs(np.sum(v0*hv1)-np.sum(v1*hv0)):.3e}")

# cross-check against the model's own full analytic Hessian
H = calc.get_hessian(atoms).reshape(3*n, 3*n)
print(f"asymmetry of get_hessian: max|H-H^T| = {np.abs(H-H.T).max():.3e}")
for i, (v_np, hv) in enumerate(exact_store):
    print(f"  v{i}: max|H@v - Hv_exact| = {np.abs(H @ v_np.ravel() - hv.ravel()).max():.3e}")

# cost
v_t = torch.as_tensor(vs[0], dtype=pos.dtype, device=pos.device)
def timeit(fn, reps=20):
    fn(); torch.cuda.synchronize(); t = time.perf_counter()
    for _ in range(reps): fn()
    torch.cuda.synchronize(); return (time.perf_counter() - t) / reps

t_force = timeit(lambda: (atoms.calc.results.clear(), atoms.get_forces()))
def one_exact():
    p, g, _ = graph_gradient()
    torch.autograd.grad((g * v_t).sum(), p)
t_exact = timeit(one_exact, reps=10)
print(f"\ncost per force eval        {t_force*1e3:8.2f} ms")
print(f"cost per exact Hv (full)   {t_exact*1e3:8.2f} ms   = {t_exact/t_force:.2f} force evals")
print(f"cost per FD Hv (2 forces)  {2*t_force*1e3:8.2f} ms")
print(f"peak CUDA mem {torch.cuda.max_memory_allocated()/2**20:.1f} MiB")

# full Hessian both ways, for the gold-standard comparison the doc relies on
t0 = time.perf_counter(); calc.get_hessian(atoms); t_H = time.perf_counter() - t0
print(f"full analytic Hessian      {t_H*1e3:8.2f} ms   = {t_H/t_force:.1f} force evals "
      f"(FD full Hessian would be {6*n} = {6*n*t_force*1e3:.0f} ms)")
