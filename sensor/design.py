"""Element-set selection rule of PLAN.md Amendment 1 (A1.2).

J(R) = min over design-grid axis pairs with separation > X of the symmetrised Poisson KL divergence
D(n, n') = 0.5 * sum_j (lam_j - lam'_j) ln(lam_j / lam'_j),  lam_j = N s_j(n).
"""
from __future__ import annotations

import numpy as np

from sensor.poisson import rotvec

X_DEG = 10.0
N_DESIGN = 1e8


def build_rots(params, K):
    """R_1 = I; remaining K-1 rotations from rotation vectors (3 parameters each)."""
    p = np.asarray(params, float).reshape(K - 1, 3)
    return [np.eye(3)] + [rotvec(v) for v in p]


def s_matrix(phi, phibar, rots, grid, eps=1.0):
    s = np.empty((len(grid), len(rots)))
    for j, R in enumerate(rots):
        s[:, j] = 1.0 + eps * (phi(grid @ R) - phibar) / phibar
    return s


def separation_mask(grid, x_deg=X_DEG):
    return np.abs(grid @ grid.T) < np.cos(np.deg2rad(x_deg))


def min_kl(s, mask, N=N_DESIGN, chunk=400):
    best = np.inf
    for a in range(0, len(s), chunk):
        sa = s[a : a + chunk, None, :]
        d = 0.5 * N * ((sa - s[None]) * np.log(sa / s[None])).sum(-1)
        d = np.where(mask[a : a + chunk], d, np.inf)
        best = min(best, d.min())
    return float(best)


def objective(params, K, phi, phibar, grid, mask, eps=1.0, N=N_DESIGN):
    rots = build_rots(params, K)
    s = s_matrix(phi, phibar, rots, grid, eps)
    if np.any(s <= 0):
        return 0.0
    return min_kl(s, mask, N)


def mirror_dmin(phi_eval, rots, dirs, x_deg=X_DEG, N=N_DESIGN):
    """PLAN A1.5: min over directions and the three sign flips of D(n, s*n), ang(n, s*n) >= X.
    phi_eval: callable (B,3)->(B,) giving the exact yield (eigensolver)."""
    flips = [np.array([-1, -1, 1.0]), np.array([1, -1, -1.0]), np.array([-1, 1, -1.0])]
    from model.sphere import fibonacci_sphere  # noqa: F401
    phibar = float(np.mean(phi_eval(_ref_grid())))
    best, per_flip = np.inf, []
    for f in flips:
        sd = dirs * f
        ang = np.degrees(np.arccos(np.clip(np.abs((dirs * sd).sum(1)), -1, 1)))
        keep = ang >= x_deg
        a = np.stack([1 + (phi_eval(dirs[keep] @ R) - phibar) / phibar for R in rots], 1)
        b = np.stack([1 + (phi_eval(sd[keep] @ R) - phibar) / phibar for R in rots], 1)
        d = 0.5 * N * ((a - b) * np.log(a / b)).sum(1)
        per_flip.append(float(d.min()) if keep.any() else np.inf)
        best = min(best, per_flip[-1])
    return best, per_flip


_REF = None


def _ref_grid():
    global _REF
    if _REF is None:
        from model.sphere import fibonacci_sphere
        _REF = fibonacci_sphere(3000)
    return _REF


# ---------------------------------------------------------------------- Amendment 1.1: refined objective J_ref
def _pair_rows(s, rows, mask, N):
    sa = s[rows][:, None, :]
    d = 0.5 * N * ((sa - s[None]) * np.log(sa / s[None])).sum(-1)
    return np.where(mask[rows], d, np.inf)


def refined_jref(phi, phibar, rots, grid, mask=None, N=N_DESIGN, n_worst=40, n_cand=3, x_deg=X_DEG):
    """Stage B objective of Amendment 1.1: continuous minimisation of D(n, n') over n' for the worst design directions."""
    from scipy.optimize import minimize
    from sensor.poisson import tangent_basis

    if mask is None:
        mask = separation_mask(grid, x_deg)
    s = s_matrix(phi, phibar, rots, grid)
    if np.any(s <= 0):
        return 0.0
    row_min = np.empty(len(grid))
    for a in range(0, len(grid), 400):
        rows = np.arange(a, min(a + 400, len(grid)))
        row_min[rows] = _pair_rows(s, rows, mask, N).min(1)
    worst = np.argsort(row_min)[:n_worst]
    D = _pair_rows(s, worst, mask, N)
    xr = np.deg2rad(x_deg)

    def s_point(v):
        return np.array([1.0 + (phi((v @ R)[None])[0] - phibar) / phibar for R in rots])

    best = np.inf
    for k, i in enumerate(worst):
        n = grid[i]
        sn = s[i]
        e1, e2 = tangent_basis(n[None])
        e1, e2 = e1[0], e2[0]

        def f(x):
            psi, chi = x
            v = np.cos(psi) * n + np.sin(psi) * (np.cos(chi) * e1 + np.sin(chi) * e2)
            sv = s_point(v)
            if np.any(sv <= 0):
                return 1e30
            return float(0.5 * N * ((sn - sv) * np.log(sn / sv)).sum())

        for c in np.argsort(D[k])[:n_cand]:
            v = grid[c] if grid[c] @ n >= 0 else -grid[c]
            psi0 = float(np.arccos(np.clip(v @ n, -1, 1)))
            chi0 = float(np.arctan2(v @ e2, v @ e1))
            r = minimize(f, [max(psi0, xr), chi0], method="L-BFGS-B", bounds=[(xr, np.pi / 2), (None, None)])
            best = min(best, r.fun, f([max(psi0, xr), chi0]))
    return float(best)
