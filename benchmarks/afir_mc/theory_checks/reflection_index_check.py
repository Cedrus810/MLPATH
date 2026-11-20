import numpy as np
rng = np.random.default_rng(1)
bad = 0; stable_by_index = {0:0,1:0,2:0,3:0}
for trial in range(40000):
    d = rng.integers(3, 9); nneg = rng.integers(0, min(4, d))
    ev = np.concatenate([-10**rng.uniform(-3, 1, nneg), 10**rng.uniform(-3, 2, d-nneg)])
    t = rng.normal(size=d); t /= np.linalg.norm(t)
    A = (np.eye(d) - 2*np.outer(t, t)) @ np.diag(ev)
    nunst = (np.linalg.eigvals(A).real < 0).sum()
    if nunst == 0: stable_by_index[nneg] += 1
    if nneg >= 2 and nunst < nneg - 1: bad += 1
print("index>=2 with fewer than n-1 unstable eigs:", bad)
print("count of linearly stable fixed points by Hessian index:", stable_by_index)
# Idea 2: explicit descent step limited by stiffest mode -> e^{mu t} gain after n steps
for ratio in (10, 50):
    mu, lst = 1.0, ratio
    dt = 1.0/lst            # stable explicit step ~ 1/lambda_stiff
    for n in (4, 10, 30):
        gain = (1+mu*dt)**n   # reaction-mode amplification per side
        th0 = np.radians(75); th = np.degrees(np.arctan(np.tan(th0)/gain))
        print(f"stiff/mu={ratio:3d} steps={n:3d} gain={gain:.2f} 75deg->{th:.1f}deg (soft mode, lambda_soft->0)")
