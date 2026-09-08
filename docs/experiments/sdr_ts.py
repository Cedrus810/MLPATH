"""Same reliability diagnostics AT A SADDLE POINT, on PHYSICAL modes.

The first probe had two flaws: the "lowest" modes it picked were the projector's own
zero modes (rigid-body), not physical ones, and a minimum of a small well-covered
molecule is where an MLP is at its best. Curvature reliability matters at the saddle
and along the UNSTABLE mode. The syn OH-torsion barrier top is reachable as a
dihedral-constrained minimum.
"""
import numpy as np, torch
from ase.build import molecule
from ase.optimize import LBFGS
from ase.constraints import FixInternals
from mace.calculators import MACECalculator

MODEL = "/home/kasuga/.cache/mace/MACE-OFF24_medium.model"
calc = MACECalculator(model_paths=MODEL, device="cuda", default_dtype="float64")

atoms = molecule("CH3CH2OH")
sym = atoms.get_chemical_symbols()
d = atoms.get_all_distances()
o = sym.index("O")
h_o = min([i for i in range(len(atoms)) if sym[i] == "H"], key=lambda i: d[o, i])
c1 = min([i for i in range(len(atoms)) if sym[i] == "C"], key=lambda i: d[o, i])
c2 = [i for i in range(len(atoms)) if sym[i] == "C" and i != c1][0]
print(f"torsion C{c2}-C{c1}-O{o}-H{h_o}")

atoms.calc = calc
atoms.set_dihedral(c2, c1, o, h_o, 0.0, indices=[h_o])
atoms.set_constraint(FixInternals(dihedrals_deg=[[0.0, [c2, c1, o, h_o]]]))
LBFGS(atoms, logfile=None).run(fmax=1e-4, steps=800)
atoms.set_constraint()
F = atoms.get_forces(); E0 = atoms.get_potential_energy()
print(f"saddle candidate: fmax = {np.abs(F).max():.3e} eV/A  dihedral = "
      f"{atoms.get_dihedral(c2, c1, o, h_o):.3f} deg  E = {E0:.9f} eV")

n = len(atoms); masses = atoms.get_masses()
root3 = np.sqrt(np.repeat(masses, 3)); inv_root = 1.0/root3

def trivial_basis(pos, m):
    r = np.sqrt(m); modes = []
    for ax in range(3):
        z = np.zeros_like(pos); z[:, ax] = r; modes.append(z.ravel())
    com = pos - np.average(pos, axis=0, weights=m)
    for ax in range(3):
        e = np.zeros(3); e[ax] = 1.0
        modes.append((np.cross(com, e) * r[:, None]).ravel())
    b = []
    for mo in modes:
        for k in b: mo = mo - (mo @ k) * k
        nn = np.linalg.norm(mo)
        if nn > 1e-8: b.append(mo/nn)
    return np.array(b)

H = np.asarray(calc.get_hessian(atoms), dtype=float).reshape(3*n, 3*n)
H = 0.5*(H + H.T)
Hm = H/np.outer(root3, root3)
B = trivial_basis(atoms.positions, masses)
print(f"\nD1 trivial-mode floor: max|Rayleigh| = {max(abs(b @ Hm @ b) for b in B):.3e}")

lam_all, vec_all = np.linalg.eigh(Hm)
overlap = np.linalg.norm(B @ vec_all, axis=0)          # how much each mode is rigid-body
phys = [i for i in range(3*n) if overlap[i] < 0.5]
lam, vec = lam_all[phys], vec_all[:, phys]
o2 = np.argsort(lam); lam, vec = lam[o2], vec[:, o2]
print(f"physical modes: {len(lam)}   6 lowest: {lam[:6]}")
print(f"negative (tol 1e-3): {[f'{v:+.6f}' for v in lam if v < -1e-3]}")
scale = 1.0545718e-34*1e10/np.sqrt(1.602176634e-19*1.66053907e-27)/100.0/2.99792458e8*1e2
icm = 1.0/(2*np.pi)  # placeholder, report eigenvalues only

grid = np.array([1e-3, 3e-3, 1e-2, 2e-2, 5e-2, 1e-1, 1.5e-1, 2e-1, 3e-1])
picks = [("UNSTABLE lam1", 0), ("lowest positive", int(np.argmax(lam > 1e-3))),
         ("2nd positive", int(np.argmax(lam > 1e-3))+1)]
ref = atoms.positions.copy()
print("\n=== D2 at the saddle: kappa_a on PHYSICAL modes, eV/(A^2 amu) ===")
store = {}
for name, idx in picks:
    u = vec[:, idx]; lam_u = float(u @ Hm @ u)
    dR = (u*inv_root).reshape(-1, 3)
    kF, kE = [], []
    for a in grid:
        vals = {}
        for s in (1, -1):
            atoms.set_positions(ref + s*a*dR)
            vals[s] = (atoms.get_potential_energy(), atoms.get_forces().copy())
        atoms.set_positions(ref)
        kE.append((vals[1][0]+vals[-1][0]-2*E0)/a**2)
        kF.append(float((-(vals[1][1]-vals[-1][1]).ravel()*inv_root) @ u)/(2*a))
    store[name] = (lam_u, np.array(kF), np.array(kE))
    print(f"\n  {name:<16} lambda = {lam_u:+.6f}")
    print("    a      " + "".join(f"{x:>10.0e}" for x in grid))
    print("    kap_F  " + "".join(f"{x:>10.5f}" for x in kF))
    print("    kap_E  " + "".join(f"{x:>10.5f}" for x in kE))
    print("    dev_F% " + "".join(f"{100*(x-lam_u)/abs(lam_u):>10.2f}" for x in kF))

# regression: kap_F(a) = k0 + c a^2 fitted on the band, vs the analytic lambda
print("\n=== C-estimator check: least squares kap_F = k0 + c a^2 over bands ===")
for name, (lam_u, kF, kE) in store.items():
    for lo, hi in [(1e-3, 2e-2), (1e-2, 1e-1), (5e-2, 3e-1)]:
        m = (grid >= lo) & (grid <= hi)
        if m.sum() < 3: continue
        A = np.vstack([np.ones(m.sum()), grid[m]**2]).T
        coef, res, *_ = np.linalg.lstsq(A, kF[m], rcond=None)
        pred = A @ coef
        dof = max(1, m.sum()-2)
        s2 = float(((kF[m]-pred)**2).sum()/dof)
        cov = s2*np.linalg.inv(A.T @ A)
        print(f"  {name:<16} band [{lo:.0e},{hi:.0e}]  k0 = {coef[0]:+.6f} +- {np.sqrt(cov[0,0]):.2e}"
              f"   c = {coef[1]:+.4f}   k0-lambda = {coef[0]-lam_u:+.3e}")

# D3 roughness along the unstable mode
u = vec[:, 0]; v3 = (u*inv_root); v3 = (v3/np.linalg.norm(v3)).reshape(-1, 3)
def exact_hv(pos, vec3):
    saved = atoms.positions.copy(); atoms.set_positions(pos)
    bd = calc._atoms_to_batch(atoms).clone().to_dict()
    p = bd["positions"]; p.requires_grad_(True)
    out = calc.models[0](bd, compute_force=False, compute_stress=False, training=True)
    g = torch.autograd.grad(out["energy"].sum(), p, create_graph=True)[0]
    vt = torch.as_tensor(vec3, dtype=p.dtype, device=p.device)
    hv = torch.autograd.grad((g*vt).sum(), p)[0].detach().cpu().numpy()
    atoms.set_positions(saved); return hv
base = exact_hv(ref.copy(), v3)
w = np.random.default_rng(0).normal(size=(n, 3)); w /= np.linalg.norm(w)
print(f"\n=== D3 along the UNSTABLE mode: ||Hv|| = {np.linalg.norm(base):.4f} eV/A^2 ===")
for dl in (1e-5, 1e-4, 1e-3, 1e-2, 1e-1):
    hv = exact_hv(ref + dl*w, v3)
    r = np.linalg.norm(hv-base)/dl
    print(f"  delta={dl:.0e}  ||dHv||/delta = {r:9.4f}  relative per A = {r/np.linalg.norm(base):8.4f}")

# D4 at the saddle: does float32 preserve the sign and size of lambda_1?
c32 = MACECalculator(model_paths=MODEL, device="cuda", default_dtype="float32")
a32 = atoms.copy(); a32.calc = c32
H32 = np.asarray(c32.get_hessian(a32), dtype=float).reshape(3*n, 3*n)
H32 = 0.5*(H32+H32.T)/np.outer(root3, root3)
l32 = np.linalg.eigh(H32)[0]
neg64 = sorted([v for v in lam_all if v < -1e-3]); neg32 = sorted([v for v in l32 if v < -1e-3])
print(f"\n=== D4 at the saddle ===")
print(f"  float64 negatives: {[f'{v:+.6f}' for v in neg64]}")
print(f"  float32 negatives: {[f'{v:+.6f}' for v in neg32]}")
print(f"  max |dlambda| = {np.abs(np.sort(lam_all)-np.sort(l32)).max():.3e}")
print(f"  float32 trivial floor = {max(abs(b @ H32 @ b) for b in B):.3e}")
