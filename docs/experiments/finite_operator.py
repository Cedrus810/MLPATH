"""Is there an OPERATOR at finite a, i.e. do "eigenvalues of L_a" exist?

Expansions (F = 4th derivative of E, fully symmetric; T = 3rd):
  D_a(v) = [g(x+av) - g(x-av)]/(2a) = H v + (a^2/6) F[v,v,v] + O(a^4)
  => <v1, D_a v2> - <v2, D_a v1> = (a^2/6)( F[v1,v2,v2,v2] - F[v2,v1,v1,v1] ) != 0
  kappa_a^E(v) = [E(x+av)+E(x-av)-2E(x)]/a^2 = <v,Hv> + (a^2/12) F[v,v,v,v] + O(a^4)
     -> a QUARTIC form on the sphere, not a quadratic one.

T1  asymmetry of the force-form bilinear map      (expect ~ a^2)
T2  non-quadraticity of the energy form            (expect ~ a^2)
T3  does the bottom mode ROTATE with a?  scan theta in the (v1,v2) plane per a.
"""
import numpy as np, torch
from ase.build import molecule
from ase.optimize import LBFGS
from ase.constraints import FixInternals
from mace.calculators import MACECalculator

calc = MACECalculator(model_paths="/home/kasuga/.cache/mace/MACE-OFF24_medium.model",
                      device="cuda", default_dtype="float64")
atoms = molecule("CH3CH2OH"); sym = atoms.get_chemical_symbols(); d = atoms.get_all_distances()
o = sym.index("O")
h = min([i for i in range(len(atoms)) if sym[i]=="H"], key=lambda i: d[o,i])
c1 = min([i for i in range(len(atoms)) if sym[i]=="C"], key=lambda i: d[o,i])
c2 = [i for i in range(len(atoms)) if sym[i]=="C" and i!=c1][0]
atoms.calc = calc
atoms.set_dihedral(c2,c1,o,h,0.0,indices=[h])
atoms.set_constraint(FixInternals(dihedrals_deg=[[0.0,[c2,c1,o,h]]]))
LBFGS(atoms, logfile=None).run(fmax=1e-4, steps=800)
atoms.set_constraint()
E0 = atoms.get_potential_energy(); ref = atoms.positions.copy()
n = len(atoms); masses = atoms.get_masses()
root3 = np.sqrt(np.repeat(masses,3)); inv = 1.0/root3
print(f"saddle: fmax={np.abs(atoms.get_forces()).max():.2e}  E={E0:.9f}")

def triv(pos,m):
    r=np.sqrt(m); mo=[]
    for ax in range(3):
        z=np.zeros_like(pos); z[:,ax]=r; mo.append(z.ravel())
    com=pos-np.average(pos,axis=0,weights=m)
    for ax in range(3):
        e=np.zeros(3); e[ax]=1.0; mo.append((np.cross(com,e)*r[:,None]).ravel())
    b=[]
    for x in mo:
        for k in b: x=x-(x@k)*k
        if np.linalg.norm(x)>1e-8: b.append(x/np.linalg.norm(x))
    return np.array(b)

H = np.asarray(calc.get_hessian(atoms),dtype=float).reshape(3*n,3*n); H=0.5*(H+H.T)
Hm = H/np.outer(root3,root3); B = triv(ref,masses)
la,va = np.linalg.eigh(Hm); ov=np.linalg.norm(B@va,axis=0)
ph=[i for i in range(3*n) if ov[i]<0.5]; la,va=la[ph],va[:,ph]
so=np.argsort(la); la,va=la[so],va[:,so]
v1,v2 = va[:,0], va[:,1]; vS = va[:,-1]
print(f"lambda_1={la[0]:+.6f}  lambda_2={la[1]:+.6f}  lambda_max={la[-1]:+.4f}")

def E_at(u_mw, a):
    atoms.set_positions(ref + a*(u_mw*inv).reshape(-1,3))
    e = atoms.get_potential_energy(); atoms.set_positions(ref); return e
def g_mw_at(u_mw, a):
    atoms.set_positions(ref + a*(u_mw*inv).reshape(-1,3))
    g = -atoms.get_forces().ravel()*inv; atoms.set_positions(ref); return g
def D_a(u_mw, a):                     # force-form finite response, mass-weighted
    return (g_mw_at(u_mw,a) - g_mw_at(u_mw,-a))/(2*a)
def kapE(u_mw, a):
    return (E_at(u_mw,a) + E_at(u_mw,-a) - 2*E0)/a**2

grid = np.array([1e-2, 2e-2, 5e-2, 1e-1, 2e-1, 3e-1])

print("\n=== T1: is the force-form bilinear map symmetric? ===")
print(f"{'a':>8} {'<v1,D_a v2>':>14} {'<v2,D_a v1>':>14} {'asymmetry':>12} {'asym/a^2':>10}")
for a in grid:
    x12 = float(v1 @ D_a(v2,a)); x21 = float(v2 @ D_a(v1,a))
    print(f"{a:>8.0e} {x12:>14.8f} {x21:>14.8f} {x12-x21:>12.2e} {(x12-x21)/a**2:>10.4f}")
print("  (a symmetric operator would give identical columns at every a)")

print("\n=== T2: is the energy form a quadratic form? ===")
print("  quadratic form => kappa(cos t v1 + sin t v2) = c2 cos^2 t + 2 c1 sin t cos t + c0 sin^2 t exactly")
th = np.linspace(-0.6, 0.6, 13)
for a in [1e-2, 1e-1, 3e-1]:
    u = np.array([np.cos(t)*v1 + np.sin(t)*v2 for t in th])
    k = np.array([kapE(x, a) for x in u])
    A = np.vstack([np.cos(th)**2, 2*np.sin(th)*np.cos(th), np.sin(th)**2]).T
    coef,_,_,_ = np.linalg.lstsq(A, k, rcond=None)
    res = k - A@coef
    print(f"  a={a:.0e}  best quadratic fit residual: rms={np.sqrt((res**2).mean()):.3e} "
          f"max={np.abs(res).max():.3e}   (|kappa| range {k.min():+.4f}..{k.max():+.4f})")

print("\n=== T3: does the bottom direction rotate with a? argmin over theta ===")
def scan(plane_v, label):
    print(f"\n  plane (v1, {label}):")
    print(f"{'a':>8} {'theta* (deg)':>13} {'kappa at theta*':>16} {'kappa at 0':>12} {'gain':>11}")
    for a in grid:
        ts = np.linspace(-0.5, 0.5, 41)
        ks = np.array([kapE(np.cos(t)*v1 + np.sin(t)*plane_v, a) for t in ts])
        i = int(np.argmin(ks))
        # parabolic refine
        if 0 < i < len(ts)-1:
            y0,y1,y2 = ks[i-1],ks[i],ks[i+1]
            shift = 0.5*(y0-y2)/(y0-2*y1+y2) if (y0-2*y1+y2)!=0 else 0.0
            tstar = ts[i] + shift*(ts[1]-ts[0]); kstar = y1 - 0.25*(y0-y2)*shift
        else:
            tstar, kstar = ts[i], ks[i]
        k0 = kapE(v1, a)
        print(f"{a:>8.0e} {np.degrees(tstar):>13.3f} {kstar:>16.6f} {k0:>12.6f} {k0-kstar:>11.2e}")
scan(v2, "v2")
scan(vS, "v_stiffest")
