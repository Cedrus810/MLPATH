"""How rough is MACE-OFF24's second derivative, and does a finite-amplitude plateau exist?

Reference-free diagnostics only (no DFT available):
  D1  magnitude of the 6 trivial (translation/rotation) eigenvalues at a tight minimum
      -- these are EXACTLY zero for any exactly invariant PES, so their size is a hard
      error floor on the analytic Hessian.
  D2  kappa_a(v) vs a over decades, per mode: is there a plateau, and is the a->0 limit
      representative of the finite-a range PRRS actually uses?
  D3  roughness of H itself: ||H(R+d w)v - H(R)v|| / d, i.e. third-derivative magnitude.
  D4  float32 vs float64: same model, same geometry, compare spectra.
"""
import time, numpy as np, torch
from ase.build import molecule
from ase.optimize import LBFGS
from mace.calculators import MACECalculator

MODEL = "/home/kasuga/.cache/mace/MACE-OFF24_medium.model"
np.set_printoptions(precision=6, suppress=False, linewidth=140)

def build(dtype):
    calc = MACECalculator(model_paths=MODEL, device="cuda", default_dtype=dtype)
    return calc

def trivial_basis(pos, masses):
    root = np.sqrt(masses)
    modes = []
    for axis in range(3):
        m = np.zeros_like(pos); m[:, axis] = root; modes.append(m.ravel())
    com = pos - np.average(pos, axis=0, weights=masses)
    for axis in range(3):
        e = np.zeros(3); e[axis] = 1.0
        modes.append((np.cross(com, e) * root[:, None]).ravel())
    basis = []
    for m in modes:
        for k in basis: m = m - (m @ k) * k
        n = np.linalg.norm(m)
        if n > 1e-8: basis.append(m / n)
    return np.array(basis)

def spectrum(calc, atoms):
    n = len(atoms)
    H = np.asarray(calc.get_hessian(atoms), dtype=float).reshape(3*n, 3*n)
    H = 0.5*(H + H.T)
    masses = atoms.get_masses()
    root = np.sqrt(np.repeat(masses, 3))
    Hm = H / np.outer(root, root)
    B = trivial_basis(atoms.positions, masses)
    P = np.eye(3*n) - B.T @ B
    lam_raw, vec_raw = np.linalg.eigh(Hm)                      # no projection
    lam, vec = np.linalg.eigh(P @ Hm @ P)                      # projected
    triv = np.sort([float(b @ Hm @ b) for b in B])             # Rayleigh quotient on trivial modes
    return Hm, lam, vec, lam_raw, triv

# ---------- tight minimum in float64 ----------
atoms = molecule("CH3CH2OH")
calc64 = build("float64")
atoms.calc = calc64
LBFGS(atoms, logfile=None).run(fmax=1e-4, steps=500)
E0 = atoms.get_potential_energy()
n = len(atoms); masses = atoms.get_masses()
print(f"minimum: fmax = {np.abs(atoms.get_forces()).max():.3e} eV/A   E = {E0:.9f} eV")
print("model submodules:", sorted({type(m).__name__ for m in calc64.models[0].modules()
                                   if 'Cutoff' in type(m).__name__ or 'Transform' in type(m).__name__
                                   or 'ZBL' in type(m).__name__}))

Hm, lam, vec, lam_raw, triv = spectrum(calc64, atoms)
print("\n=== D1: trivial-mode error floor (eV/(A^2 amu)) ===")
print("Rayleigh quotients on the 6 exact-zero modes:", np.array(triv))
print("unprojected spectrum, 8 smallest |.|:", np.array(sorted(lam_raw, key=abs))[:8])
phys = lam[np.abs(lam) > 0]
print(f"projected physical spectrum, 6 lowest: {lam[:6]}")
print(f"|trivial|max = {np.abs(triv).max():.3e}   vs lowest physical |lambda| = "
      f"{np.abs(lam[np.abs(lam)>1e-12]).min():.3e}   vs eigenvalue_tol = 1e-3")

# ---------- D2: kappa_a vs a, in mass-weighted coordinates ----------
order = np.argsort(lam)
picks = [("lowest", order[0]), ("2nd lowest", order[1]),
         ("mid", order[len(order)//2]), ("stiffest", order[-1])]
inv_root = 1.0/np.sqrt(np.repeat(masses, 3))
grid = np.array([1e-4, 3e-4, 1e-3, 3e-3, 1e-2, 2e-2, 5e-2, 1e-1, 2e-1, 3e-1])

print("\n=== D2: kappa_a(v) in eV/(A^2 amu), mass-weighted amplitude a in A*sqrt(amu) ===")
for name, idx in picks:
    u = vec[:, idx]
    lam_u = float(u @ Hm @ u)
    dR = (u * inv_root).reshape(-1, 3)          # cartesian displacement per unit mass-wt amp
    ref = atoms.positions.copy()
    rowE, rowF = [], []
    for a in grid:
        vals = {}
        for sign in (1, -1):
            atoms.set_positions(ref + sign * a * dR)
            vals[sign] = (atoms.get_potential_energy(), atoms.get_forces().copy())
        atoms.set_positions(ref)
        rowE.append((vals[1][0] + vals[-1][0] - 2*E0) / a**2)
        gmw = -(vals[1][1] - vals[-1][1]).ravel() * inv_root
        rowF.append(float(gmw @ u) / (2*a))
    print(f"\n  mode={name:<11} lambda(analytic) = {lam_u:+.6f}")
    print("    a      " + "".join(f"{a:>11.0e}" for a in grid))
    print("    kap_F  " + "".join(f"{v:>11.5f}" for v in rowF))
    print("    kap_E  " + "".join(f"{v:>11.5f}" for v in rowE))
    print("    dev_F% " + "".join(f"{100*(v-lam_u)/abs(lam_u):>11.3f}" for v in rowF))

# ---------- D3: roughness of H (third derivative) ----------
print("\n=== D3: how fast does H change? ||H(R+d w)v - H(R)v|| / d ===")
rng = np.random.default_rng(0)
w = rng.normal(size=(n, 3)); w /= np.linalg.norm(w)
v = vec[:, order[0]] * inv_root; v = (v/np.linalg.norm(v)).reshape(-1, 3)
def exact_hv(pos, vec3):
    saved = atoms.positions.copy(); atoms.set_positions(pos)
    batch = calc64._atoms_to_batch(atoms); bd = batch.clone().to_dict()
    p = bd["positions"]; p.requires_grad_(True)
    out = calc64.models[0](bd, compute_force=False, compute_stress=False, training=True)
    g = torch.autograd.grad(out["energy"].sum(), p, create_graph=True)[0]
    vt = torch.as_tensor(vec3, dtype=p.dtype, device=p.device)
    hv = torch.autograd.grad((g*vt).sum(), p)[0].detach().cpu().numpy()
    atoms.set_positions(saved); return hv
base = exact_hv(atoms.positions.copy(), v)
print(f"  ||Hv|| = {np.linalg.norm(base):.4f} eV/A^2")
print(f"  {'delta (A)':>10} {'||dHv||/delta':>15} {'rel. per A':>12}")
for d in (1e-5, 1e-4, 1e-3, 1e-2, 1e-1):
    hv = exact_hv(atoms.positions + d*w, v)
    nrm = np.linalg.norm(hv - base)/d
    print(f"  {d:>10.0e} {nrm:>15.4f} {nrm/np.linalg.norm(base):>12.4f}")

# ---------- D4: float32 vs float64 ----------
print("\n=== D4: float32 vs float64 spectra ===")
calc32 = build("float32")
a32 = atoms.copy(); a32.calc = calc32
_, lam32, _, _, triv32 = spectrum(calc32, a32)
print(f"  float32 fmax at the float64 minimum: {np.abs(a32.get_forces()).max():.3e} eV/A")
print(f"  float64 6 lowest: {lam[:6]}")
print(f"  float32 6 lowest: {lam32[:6]}")
print(f"  max |dlambda| over full spectrum: {np.abs(np.sort(lam)-np.sort(lam32)).max():.3e}")
print(f"  float32 |trivial|max = {np.abs(triv32).max():.3e}  (float64: {np.abs(triv).max():.3e})")
