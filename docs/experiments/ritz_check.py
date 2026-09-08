"""Check doc A.6: is theta_2 > 0 really no proof that lambda_2 > 0?
Build a saddle-like spectrum (one clearly negative mode, one slightly negative mode)
and run plain Lanczos, watching theta_2 and its residual."""
import numpy as np
rng = np.random.default_rng(0)
N = 300
lam = np.sort(np.concatenate([[-0.50, -0.004], rng.uniform(0.01, 5.0, N-2)]))
Q, _ = np.linalg.qr(rng.normal(size=(N, N)))
H = Q @ np.diag(lam) @ Q.T
H = 0.5*(H+H.T)
print(f"true lambda_1={lam[0]:+.6f}  lambda_2={lam[1]:+.6f}  (nu_true = {int((lam<0).sum())})")

def lanczos(H, m, v0):
    n = len(v0); V = np.zeros((n, m+1)); alpha = np.zeros(m); beta = np.zeros(m+1)
    V[:, 0] = v0/np.linalg.norm(v0)
    for j in range(m):
        w = H @ V[:, j]
        alpha[j] = V[:, j] @ w
        w = w - alpha[j]*V[:, j] - (beta[j]*V[:, j-1] if j > 0 else 0)
        w -= V[:, :j+1] @ (V[:, :j+1].T @ w)   # full reorthogonalization
        beta[j+1] = np.linalg.norm(w)
        if beta[j+1] < 1e-13: m = j+1; break
        V[:, j+1] = w/beta[j+1]
    T = np.diag(alpha[:m]) + np.diag(beta[1:m], 1) + np.diag(beta[1:m], -1)
    theta, S = np.linalg.eigh(T)
    U = V[:, :m] @ S
    res = [np.linalg.norm(H @ U[:, k] - theta[k]*U[:, k]) for k in range(min(3, m))]
    return theta, res

v0 = rng.normal(size=N)
print(f"{'m':>4} {'theta_1':>10} {'|r_1|':>9} {'theta_2':>10} {'|r_2|':>9}  "
      f"{'naive nu (theta<0)':>18}  {'certified (theta_2-|r_2|>0)':>27}")
for m in (5, 10, 15, 20, 30, 40, 60, 80):
    theta, res = lanczos(H, m, v0)
    t1, t2, r1, r2 = theta[0], theta[1], res[0], res[1]
    naive = int((theta < 0).sum())
    cert = "yes -> would claim nu=1" if t2 - r2 > 0 else "no  -> undecided"
    print(f"{m:>4} {t1:>+10.5f} {r1:>9.2e} {t2:>+10.5f} {r2:>9.2e}  {naive:>18}  {cert:>27}")
