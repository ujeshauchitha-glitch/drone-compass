"""Sphere grids, extrema search and anisotropy."""
from __future__ import annotations

import numpy as np
from scipy.optimize import minimize


def fibonacci_sphere(n: int) -> np.ndarray:
    i = np.arange(n) + 0.5
    z = 1 - 2 * i / n
    phi = np.pi * (1 + 5**0.5) * i
    r = np.sqrt(1 - z * z)
    return np.stack([r * np.cos(phi), r * np.sin(phi), z], axis=1)


def fibonacci_hemisphere(n: int) -> np.ndarray:
    """n points with z >= 0 (axis representatives)."""
    i = np.arange(n) + 0.5
    z = 1 - i / n
    phi = np.pi * (1 + 5**0.5) * i
    r = np.sqrt(1 - z * z)
    return np.stack([r * np.cos(phi), r * np.sin(phi), z], axis=1)


def angles_to_vec(theta, phi):
    return np.array([np.sin(theta) * np.cos(phi), np.sin(theta) * np.sin(phi), np.cos(theta)])


def sphere_extrema(func, n_grid=2000, n_polish=5):
    """func: (N,3)->(N,). Returns dict with grid and polished max/min and their locations."""
    g = fibonacci_sphere(n_grid)
    v = func(g)
    res = {"grid_max": v.max(), "grid_min": v.min(), "mean": v.mean()}
    for sign, key in ((+1.0, "max"), (-1.0, "min")):
        order = np.argsort(-sign * v)[:n_polish]
        best = None
        for o in order:
            th0 = np.arccos(np.clip(g[o, 2], -1, 1))
            ph0 = np.arctan2(g[o, 1], g[o, 0])
            f = lambda x: -sign * float(func(angles_to_vec(x[0], x[1])[None])[0])
            r = minimize(f, [th0, ph0], method="Nelder-Mead", options={"xatol": 1e-9, "fatol": 1e-14, "maxiter": 400})
            if best is None or r.fun < best.fun:
                best = r
        res[key] = -sign * best.fun
        res[key + "_at"] = angles_to_vec(*best.x)
    res["aniso"] = res["max"] - res["min"]
    return res
