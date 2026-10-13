"""Photon-limited K-element sensor, Cramer-Rao bound and Poisson maximum-likelihood axis estimator (PLAN.md sec. 4)."""
from __future__ import annotations

import numpy as np

from model.sphere import fibonacci_hemisphere
from model.surrogate import tangent_project


def rot_x(a):
    c, s = np.cos(a), np.sin(a)
    return np.array([[1, 0, 0], [0, c, -s], [0, s, c]])


def rotvec(r):
    """Rotation matrix from rotation vector (Rodrigues)."""
    th = np.linalg.norm(r)
    if th < 1e-12:
        return np.eye(3)
    k = r / th
    K = np.array([[0, -k[2], k[1]], [k[2], 0, -k[0]], [-k[1], k[0], 0]])
    return np.eye(3) + np.sin(th) * K + (1 - np.cos(th)) * K @ K


CYC = np.array([[0, 0, 1], [1, 0, 0], [0, 1, 0]], float)  # cyclic axis permutation, det = +1
ELEMENT_SETS = {
    "K1": [np.eye(3)],
    "K2": [np.eye(3), rot_x(np.pi / 2)],
    "K3orth": [np.eye(3), CYC, CYC @ CYC],
}


class RenPhi:
    """Ren et al. phenomenological signal Phi(zeta) = a - b exp(-(zeta-90deg)^2/(2c^2)), zeta = angle to the flavin z axis."""

    def __init__(self, a=0.5, b=0.00127, c_deg=26.0):
        self.a, self.b, self.c = a, b, np.deg2rad(c_deg)

    def value_grad(self, n):
        n = np.atleast_2d(n)
        u = np.clip(n[:, 2], -1 + 1e-12, 1 - 1e-12)
        z = np.arccos(u)
        g = np.exp(-((z - np.pi / 2) ** 2) / (2 * self.c**2))
        val = self.a - self.b * g
        dz = self.b * g * (z - np.pi / 2) / self.c**2  # dPhi/dzeta
        du = dz * (-1 / np.sqrt(1 - u * u))            # dPhi/du
        grad = np.zeros_like(n)
        grad[:, 2] = du
        return val, grad

    def __call__(self, n):
        return self.value_grad(n)[0]

    @property
    def mean(self):
        from model.sphere import fibonacci_sphere
        return float(self(fibonacci_sphere(200000)).mean())


def tangent_basis(n):
    """(B,3) unit vectors -> e1, e2 orthonormal tangent bases (deterministic)."""
    a = np.where((np.abs(n[:, 0]) < 0.9)[:, None], np.array([1.0, 0, 0]), np.array([0, 1.0, 0]))
    e1 = np.cross(n, a)
    e1 /= np.linalg.norm(e1, axis=1, keepdims=True)
    e2 = np.cross(n, e1)
    return e1, e2


class Sensor:
    """K identical molecular elements at rotations R_j; mean counts lam_j = N s_j, s_j = 1 + eps*delta_j + beta."""

    def __init__(self, phi, rots, eps=1.0, beta=0.0, rho=0.0, sigma_g=0.0):
        self.phi, self.rots = phi, [np.asarray(r, float) for r in rots]
        self.eps, self.beta, self.rho, self.sigma_g = eps, beta, rho, sigma_g
        self.phibar = phi.mean
        self.K = len(rots)

    # shape-only quantities (independent of N)
    def delta_grad(self, n):
        n = np.atleast_2d(n)
        d = np.empty((len(n), self.K))
        g = np.empty((len(n), self.K, 3))
        for j, R in enumerate(self.rots):
            nj = n @ R  # rows of R^T n
            v, gr = self.phi.value_grad(nj)
            d[:, j] = (v - self.phibar) / self.phibar
            g[:, j] = tangent_project(gr @ R.T, n) / self.phibar
        return d, g

    def s(self, d):
        return 1.0 + self.eps * d + self.beta

    def fisher(self, n, N):
        """Fisher information in tangent coordinates (B,2,2), including the gain nuisance if sigma_g > 0."""
        n = np.atleast_2d(n)
        d, g = self.delta_grad(n)
        s = self.s(d)
        e1, e2 = tangent_basis(n)
        gt = np.stack([(g * e1[:, None]).sum(-1), (g * e2[:, None]).sum(-1)], axis=-1)  # (B,K,2)
        w = 1.0 / (s + self.rho**2)
        J = N * self.eps**2 * np.einsum("bk,bki,bkj->bij", w, gt, gt)
        if self.sigma_g > 0:
            Jnx = N * self.eps * np.einsum("bk,bk,bki->bi", s, w, gt)
            Jxx = N * (s * s * w).sum(1) + 1.0 / self.sigma_g**2
            J = J - Jnx[:, :, None] * Jnx[:, None, :] / Jxx[:, None, None]
        return J

    def crb(self, n, N):
        """Returns (rms bound sqrt(tr J^-1), median-error bound sqrt(ln2 tr J^-1), singular flag)."""
        J = self.fisher(n, N)
        ev = np.linalg.eigvalsh(J)
        singular = ev[:, 0] <= 1e-12 * np.maximum(ev[:, 1], 1e-300)
        with np.errstate(divide="ignore", invalid="ignore"):
            tr = np.where(singular, np.inf, (1.0 / ev).sum(1))
        return np.sqrt(tr), np.sqrt(np.log(2) * tr), singular

    # ------------------------------------------------------------------ data and estimator
    def lam(self, n, N):
        d, _ = self.delta_grad(n)
        return N * self.s(d)

    def simulate(self, n_true, N, rng):
        lam = self.lam(n_true[None], N)[0]
        if self.sigma_g > 0:
            lam = lam * np.exp(self.sigma_g * rng.normal())
        c = rng.poisson(lam).astype(float)
        if self.rho > 0:
            c = c + rng.normal(0.0, self.rho * np.sqrt(N), size=c.shape)
        return c

    def _loglik(self, c, lam, N):
        """Log-likelihood for rows of lam (B,K); Poisson (rho = 0) or Gaussian with variance lam + rho^2 N."""
        lam = np.maximum(lam, 1e-12)
        if self.rho > 0:
            v = lam + self.rho**2 * N
            return -0.5 * (((c - lam) ** 2) / v + np.log(v)).sum(1)
        if self.sigma_g > 0:  # profile over x = ln g with a N(0, sigma_g^2) prior (Newton)
            C, Lam = c.sum(), lam.sum(1)
            x = np.log(np.maximum(C, 1e-9) / Lam)
            for _ in range(12):
                f1 = C - np.exp(x) * Lam - x / self.sigma_g**2
                f2 = -np.exp(x) * Lam - 1.0 / self.sigma_g**2
                x = x - f1 / f2
            return (c * (np.log(lam) + x[:, None])).sum(1) - np.exp(x) * Lam - x**2 / (2 * self.sigma_g**2)
        return (c * np.log(lam) - lam).sum(1)

    def make_grid(self, n_grid=4000):
        G = fibonacci_hemisphere(n_grid)
        d, _ = self.delta_grad(G)
        self._G, self._Gs = G, self.s(d)
        return G

    def mle(self, c, N, n_refine=4, iters=25, sep_deg=8.0):
        """Grid search, then Fisher-scoring refinement of up to n_refine grid candidates that are > sep_deg apart.

        (CHANGELOG 2026-10-13: the first version refined only the single best grid point, which at very large N can
        sit in the wrong basin because the likelihood peak is far narrower than the grid spacing.)"""
        if not hasattr(self, "_G"):
            self.make_grid()
        ll = self._loglik(c, N * self._Gs, N)
        best, cos_sep = [], np.cos(np.deg2rad(sep_deg))
        for b in np.argsort(-ll)[:400]:
            if all(abs(self._G[b] @ self._G[a]) < cos_sep for a in best):
                best.append(b)
            if len(best) == n_refine:
                break
        out_n, out_ll = None, -np.inf
        for b in best:
            n = self._G[b].copy()
            cur = self._loglik(c, self.lam(n[None], N), N)[0]
            for _ in range(iters):
                d, g = self.delta_grad(n[None])
                s = self.s(d)[0]
                lam = N * s
                e1, e2 = tangent_basis(n[None])
                gt = np.stack([g[0] @ e1[0], g[0] @ e2[0]], axis=-1)  # (K,2)
                w = 1.0 / (s + self.rho**2)
                if self.rho > 0:
                    resid = (c - lam) / (lam + self.rho**2 * N)
                else:
                    resid = c / lam - 1.0
                    if self.sigma_g > 0:
                        resid = c / (lam) - 1.0  # approximate score (prior ignored in refinement; checked by logL)
                score = N * self.eps * (resid[:, None] * gt).sum(0)
                J = N * self.eps**2 * (w[:, None, None] * gt[:, :, None] * gt[:, None, :]).sum(0)
                J = J + 1e-9 * np.trace(J) * np.eye(2) + 1e-30 * np.eye(2)
                step = np.linalg.solve(J, score)
                t = 1.0
                improved = False
                for _ in range(8):
                    cand = n + t * (step[0] * e1[0] + step[1] * e2[0])
                    cand /= np.linalg.norm(cand)
                    new = self._loglik(c, self.lam(cand[None], N), N)[0]
                    if new > cur + 1e-12:
                        n, cur, improved = cand, new, True
                        break
                    t *= 0.5
                if not improved:
                    break
            if cur > out_ll:
                out_n, out_ll = n, cur
        return out_n


def axis_error_deg(n_hat, n_true):
    return np.degrees(np.arccos(np.clip(abs(float(n_hat @ n_true)), -1, 1)))


def hemisphere_directions(n, rng):
    z = rng.uniform(0, 1, n)
    ph = rng.uniform(0, 2 * np.pi, n)
    r = np.sqrt(1 - z * z)
    return np.stack([r * np.cos(ph), r * np.sin(ph), z], axis=1)
