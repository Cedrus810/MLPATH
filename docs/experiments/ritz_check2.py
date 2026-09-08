"""Which criterion survives? Test three spectra with a shared Lanczos driver."""
import numpy as np
rng = np.random.default_rng(0)
N = 300

def build(l1, l2):
    lam = np.sort(np.concatenate([[l1, l2], rng.uniform(0.05, 5.0, N-2)]))
    Q, _ = np.linalg.qr(np.random.default_rng(7).normal(size=(N, N)))
    H = Q @ np.diag(lam) @ Q.T
    return 0.5*(H+H.T), lam

def lanczos(H, m, v0):
    n = len(v0); V = np.zeros((n, m+1)); alpha = np.zeros(m); beta = np.zeros(m+1)
    V[:, 0] = v0/np.linalg.norm(v0)
    for j in range(m):
        w = H @ V[:, j]; alpha[j] = V[:, j] @ w
        w = w - alpha[j]*V[:, j] - (beta[j]*V[:, j-1] if j > 0 else 0)
        w -= V[:, :j+1] @ (V[:, :j+1].T @ w)
        beta[j+1] = np.linalg.norm(w)
        if beta[j+1] < 1e-13: m = j+1; break
        V[:, j+1] = w/beta[j+1]
    T = np.diag(alpha[:m]) + np.diag(beta[1:m], 1) + np.diag(beta[1:m], -1)
    theta, S = np.linalg.eigh(T); U = V[:, :m] @ S
    res = [np.linalg.norm(H @ U[:, k] - theta[k]*U[:, k]) for k in range(min(3, m))]
    return theta[:3], res

TOL = 1e-3
for name, (l1, l2) in [("true saddle, clean gap", (-0.50, +0.050)),
                       ("true saddle, tight gap", (-0.50, +0.004)),
                       ("2nd-order saddle",       (-0.50, -0.004))]:
    H, lam = build(l1, l2)
    nu = int((lam < -TOL).sum())
    print(f"\n=== {name}: lambda_1={lam[0]:+.4f} lambda_2={lam[1]:+.4f}  nu_true={nu} ===")
    print(f"{'m':>4} {'theta_1':>10} {'theta_2':>10} {'|r_2|':>9}  "
          f"{'A: th2-|r2|>-tol':>17} {'B: th2-10|r2|>-tol':>19} {'C: count stable':>16}")
    prev = None
    v0 = np.random.default_rng(3).normal(size=N)
    for m in (5, 10, 15, 20, 30, 40, 60, 80, 120):
        theta, res = lanczos(H, m, v0)
        t2, r2 = theta[1], res[1]
        cnt = int((theta < -TOL).sum())
        A = "PASS nu=1" if t2 - r2 > -TOL else "-"
        B = "PASS nu=1" if t2 - 10*r2 > -TOL else "-"
        C = f"{cnt} (was {prev})" if prev is not None else f"{cnt}"
        wrong = " <-- WRONG" if (A.startswith("PASS") and nu != 1) else ""
        wrongB = " <-- WRONG" if (B.startswith("PASS") and nu != 1) else ""
        print(f"{m:>4} {theta[0]:>+10.5f} {t2:>+10.5f} {r2:>9.2e}  {A:>17}{wrong} {B:>19}{wrongB} {C:>16}")
        prev = cnt
