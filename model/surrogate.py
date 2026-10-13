"""Even-degree real spherical-harmonic surrogate of the singlet yield Phi(n) (Phi(n) = Phi(-n))."""
from __future__ import annotations

import numpy as np
from scipy.special import sph_harm_y

from .sphere import fibonacci_sphere


def sh_basis(vecs: np.ndarray, L: int) -> np.ndarray:
    v = np.atleast_2d(vecs)
    theta = np.arccos(np.clip(v[:, 2], -1, 1))
    phi = np.arctan2(v[:, 1], v[:, 0])
    cols = []
    for l in range(0, L + 1, 2):
        for m in range(0, l + 1):
            Y = sph_harm_y(l, m, theta, phi)
            cols.append(Y.real)
            if m > 0:
                cols.append(Y.imag)
    return np.stack(cols, axis=1)


class SHSurrogate:
    def __init__(self, coef: np.ndarray, L: int):
        self.coef = coef
        self.L = L

    @classmethod
    def fit(cls, func, L: int, n_fit: int = 6000):
        g = fibonacci_sphere(n_fit)
        A = sh_basis(g, L)
        coef, *_ = np.linalg.lstsq(A, func(g), rcond=None)
        return cls(coef, L)

    def __call__(self, vecs):
        return sh_basis(vecs, self.L) @ self.coef

    @property
    def mean(self) -> float:
        """Sphere average: only Y00 = 1/sqrt(4 pi) has a non-zero mean."""
        return float(self.coef[0] / np.sqrt(4 * np.pi))


def choose_degree(func, degrees=(8, 12, 16, 20, 24, 32), tol=1e-7, seed=12345, n_holdout=1000):
    """Smallest even degree whose held-out absolute error is below tol. Returns (surrogate, table)."""
    rng = np.random.default_rng(seed)
    h = rng.normal(size=(n_holdout, 3))
    h /= np.linalg.norm(h, axis=1, keepdims=True)
    truth = func(h)
    table = []
    chosen = None
    for L in degrees:
        s = SHSurrogate.fit(func, L)
        err = float(np.max(np.abs(s(h) - truth)))
        table.append((L, err))
        if chosen is None and err < tol:
            chosen = s
    return chosen, table


class PolySurrogate:
    """Homogeneous even-degree-L Cartesian polynomial restricted to the unit sphere.

    Spans exactly the even spherical harmonics up to degree L, but evaluates (value + analytic gradient) much
    faster than scipy's sph_harm_y.  Used from Stage 2 onward; verified against the SH surrogate and the
    eigenbasis solver in tests/test_sensor.py.
    """

    def __init__(self, coef, L, expo):
        self.coef, self.L, self.expo = coef, L, expo

    @staticmethod
    def exponents(L):
        return np.array([(a, b, L - a - b) for a in range(L + 1) for b in range(L + 1 - a)])

    @classmethod
    def fit(cls, func, L=8, n_fit=4000):
        expo = cls.exponents(L)
        g = fibonacci_sphere(n_fit)
        A = cls._design(g, expo, L)
        coef, *_ = np.linalg.lstsq(A, func(g), rcond=None)
        return cls(coef, L, expo)

    @classmethod
    def _design(cls, n, expo, L):
        px, py, pz = (cls._pw(n[:, i], L) for i in range(3))
        return px[:, expo[:, 0]] * py[:, expo[:, 1]] * pz[:, expo[:, 2]]

    @staticmethod
    def _pw(v, L):
        P = np.ones((len(v), L + 1))
        for p in range(1, L + 1):
            P[:, p] = P[:, p - 1] * v
        return P

    def value_grad(self, n):
        """Returns (value (B,), cartesian gradient (B,3)); gradient is that of the polynomial extension
        (use tangent_project for the sphere gradient)."""
        n = np.atleast_2d(n)
        L, e = self.L, self.expo
        px, py, pz = (self._pw(n[:, i], L) for i in range(3))
        mx, my, mz = px[:, e[:, 0]], py[:, e[:, 1]], pz[:, e[:, 2]]
        c = self.coef
        val = (mx * my * mz) @ c
        ex, ey, ez = e[:, 0], e[:, 1], e[:, 2]
        dx = (ex * px[:, np.maximum(ex - 1, 0)]) * my * mz
        dy = mx * (ey * py[:, np.maximum(ey - 1, 0)]) * mz
        dz = mx * my * (ez * pz[:, np.maximum(ez - 1, 0)])
        return val, np.stack([dx @ c, dy @ c, dz @ c], axis=1)

    def __call__(self, n):
        return self.value_grad(n)[0]

    @property
    def mean(self):
        """Exact sphere mean: mean of x^a y^b z^c is (a-1)!!(b-1)!!(c-1)!!/(a+b+c+1)!! for even a,b,c, else 0."""
        def df(m):  # double factorial with (-1)!! = 1
            out = 1
            while m > 1:
                out *= m
                m -= 2
            return out
        w = np.array([df(a - 1) * df(b - 1) * df(c - 1) / df(a + b + c + 1) if (a % 2 == b % 2 == c % 2 == 0) else 0.0
                      for a, b, c in self.expo])
        return float(w @ self.coef)


def tangent_project(grad, n):
    return grad - (grad * n).sum(-1, keepdims=True) * n
