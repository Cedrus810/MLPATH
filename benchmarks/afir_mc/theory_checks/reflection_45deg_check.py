import numpy as np
rng = np.random.default_rng(0)
def crit(lams, w, mu=1.0):
    # smallest tau1^2 at which A=(I-2tt^T)H stays Hurwitz-positive, scanning tau1^2 from 1 down
    lo = None
    for t1sq in np.linspace(0.999, 0.001, 999):
        t = np.concatenate([[np.sqrt(t1sq)], np.sqrt((1-t1sq)*w)])
        H = np.diag(np.concatenate([[-mu], lams]))
        A = (np.eye(len(t)) - 2*np.outer(t, t)) @ H
        if np.linalg.eigvals(A).real.min() <= 0:
            return t1sq
    return 0.0
viol = 0; worst = 0
for trial in range(20000):
    d = rng.integers(2, 8)
    lams = 10**rng.uniform(-3, 2, d-1)
    w = rng.dirichlet(np.ones(d-1))
    t1sq = rng.uniform(0.5001, 1)
    t = np.concatenate([[np.sqrt(t1sq)], np.sqrt((1-t1sq)*w)])
    H = np.diag(np.concatenate([[-1.0], lams]))
    A = (np.eye(d) - 2*np.outer(t, t)) @ H
    m = np.linalg.eigvals(A).real.min()
    if m <= 0: viol += 1
print("violations of tau1^2>1/2 => stable:", viol, "/ 20000")
print("2D crit (lam=0.01):", crit(np.array([0.01]), np.array([1.0])))
print("3D crit, lam=(0.001,100), w=(.5,.5):", crit(np.array([1e-3, 100.]), np.array([.5, .5])))
print("3D crit, lam=(0.01,0.01), degenerate:", crit(np.array([1e-2, 1e-2]), np.array([.5, .5])))
# constrained (HF) condition for same 3D example at tau1^2=0.6
t = np.array([np.sqrt(.6), np.sqrt(.2), np.sqrt(.2)]); print("tau^T H^-1 tau (HF needs <0):", -.6/1 + .2/1e-3 + .2/100)
