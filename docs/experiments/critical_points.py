"""A QUADRATIC form on the sphere has exactly N critical directions (its eigenvectors,
up to sign). Solve grad_v kappa_a^E(v) = mu v, ||v||=1 by Newton from many random
starts and COUNT them. More than N == no spectrum. Non-orthogonal == no eigenbasis.
"""
import math, numpy as np
rng = np.random.default_rng(1); N = 8
def symm(t,k):
    from itertools import permutations
    o=np.zeros_like(t)
    for p in permutations(range(k)): o += np.transpose(t,p)
    return o/math.factorial(k)
H = symm(rng.normal(size=(N,N)),2); F = symm(rng.normal(size=(N,)*4),4)
lam, V = np.linalg.eigh(H)

def kap(v,a):  return v@H@v + (a**2/12)*np.einsum('ijkl,i,j,k,l->',F,v,v,v,v)
def grad(v,a): return 2*H@v + (a**2/3)*np.einsum('ijkl,j,k,l->i',F,v,v,v)
def jac(v,a,mu):
    A = 2*H + (a**2)*np.einsum('ijkl,k,l->ij',F,v,v) - mu*np.eye(N)
    J = np.zeros((N+1,N+1)); J[:N,:N]=A; J[:N,N]=-v; J[N,:N]=2*v
    return J

def all_critical(a, starts=3000):
    found=[]
    for _ in range(starts):
        v = rng.normal(size=N); v/=np.linalg.norm(v); mu = float(v@grad(v,a))
        ok=False
        for _ in range(200):
            r = np.concatenate([grad(v,a)-mu*v, [v@v-1.0]])
            if np.linalg.norm(r) < 1e-12: ok=True; break
            try: step = np.linalg.solve(jac(v,a,mu), -r)
            except np.linalg.LinAlgError: break
            v = v + step[:N]; mu = mu + step[N]
            nv = np.linalg.norm(v)
            if not np.isfinite(nv) or nv > 1e6: break
        if not ok: continue
        v = v/np.linalg.norm(v)
        if not any(abs(abs(v@u)-1) < 1e-6 for u in found): found.append(v.copy())
    return found

print(f"N = {N}: a quadratic form has exactly {N} critical directions (up to sign)")
print(f"H eigenvalues: {np.array2string(lam, precision=4)}\n")
print(f"{'a':>6} {'# critical dirs':>16} {'lowest kappa':>13} {'|v1.v2| of 2 lowest':>20} {'verdict':>32}")
for a in (0.0, 0.1, 0.3, 1.0, 2.0, 3.0):
    cd = all_critical(a)
    ks = np.array([kap(v,a) for v in cd]); o = np.argsort(ks)
    dot = abs(cd[o[0]] @ cd[o[1]]) if len(cd) > 1 else float('nan')
    verdict = "spectrum (matches eigenvectors)" if len(cd)==N and dot<1e-8 else \
              ("MORE THAN N -> no spectrum" if len(cd)>N else "non-orthogonal -> no eigenbasis")
    print(f"{a:>6.1f} {len(cd):>16} {ks.min():>13.5f} {dot:>20.6f} {verdict:>32}")
