"""Same claims on an ANALYTIC potential with known H, T, F -- no MLP involved.
E(x) = 1/2 x.Hx + 1/6 T[x,x,x] + 1/24 F[x,x,x,x], random symmetric H,T,F.
Checks:
  (1) asymmetry  <v1,D_a v2> - <v2,D_a v1>  ==  (a^2/6)(F[v1,v2,v2,v2]-F[v2,v1,v1,v1])
  (2) kappa_a^E(v) == v.Hv + (a^2/12) F[v,v,v,v]
  (3) grad_v kappa_a^E(v) == 2 D_a(v)
"""
import math
import numpy as np
rng = np.random.default_rng(1); N = 8

def sym(t, k):
    out = np.zeros_like(t)
    from itertools import permutations
    for p in permutations(range(k)): out += np.transpose(t, p)
    return out/math.factorial(k)

H = sym(rng.normal(size=(N,N)), 2)
T = sym(rng.normal(size=(N,)*3), 3)
F = sym(rng.normal(size=(N,)*4), 4)

def E(x):
    return (0.5*np.einsum('ij,i,j->',H,x,x) + np.einsum('ijk,i,j,k->',T,x,x,x)/6
            + np.einsum('ijkl,i,j,k,l->',F,x,x,x,x)/24)
def g(x):
    return (H@x + 0.5*np.einsum('ijk,j,k->i',T,x,x)
            + np.einsum('ijkl,j,k,l->i',F,x,x,x)/6)
def D_a(v,a): return (g(a*v)-g(-a*v))/(2*a)
def kapE(v,a): return (E(a*v)+E(-a*v)-2*E(np.zeros(N)))/a**2

v1 = rng.normal(size=N); v1/=np.linalg.norm(v1)
v2 = rng.normal(size=N); v2-= (v2@v1)*v1; v2/=np.linalg.norm(v2)

print(f"{'a':>8} {'measured asym':>16} {'predicted asym':>16} {'ratio':>9} | "
      f"{'kapE err vs formula':>20} | {'||grad - 2 D_a||':>17}")
for a in (1e-3, 1e-2, 1e-1, 3e-1, 1.0):
    asym = float(v1 @ D_a(v2,a) - v2 @ D_a(v1,a))
    pred = (a**2/6)*(np.einsum('ijkl,i,j,k,l->',F,v1,v2,v2,v2)
                     - np.einsum('ijkl,i,j,k,l->',F,v2,v1,v1,v1))
    kap = kapE(v1,a)
    kap_formula = v1@H@v1 + (a**2/12)*np.einsum('ijkl,i,j,k,l->',F,v1,v1,v1,v1)
    # numerical grad of kapE wrt v (unconstrained), compare to 2 D_a(v)
    eps = 1e-6; num = np.zeros(N)
    for i in range(N):
        e = np.zeros(N); e[i] = eps
        num[i] = (kapE(v1+e,a)-kapE(v1-e,a))/(2*eps)
    print(f"{a:>8.0e} {asym:>16.3e} {pred:>16.3e} {asym/pred if pred else np.nan:>9.5f} | "
          f"{kap-kap_formula:>20.3e} | {np.linalg.norm(num-2*D_a(v1,a)):>17.3e}")
