"""Gate S1 checks (a)-(f) from PLAN.md section 3.4.  Thresholds are copied from PLAN.md and are not tunable."""
from __future__ import annotations

import numpy as np

from .sphere import fibonacci_sphere, sphere_extrema
from .spin import H16, RECALLED, RPModel, isotropic_version, primary_model, yield_from_hamiltonian
from .surrogate import choose_degree

TOL_A = 1e-9
TOL_B = 1e-9
TOL_D = 1e-9
TOL_E = 1e-8
TOL_F_GRID = 1e-6
TOL_F_SH = 1e-7
ANISO_MIN = 1e-3


def rand_dirs(n, seed):
    rng = np.random.default_rng(seed)
    v = rng.normal(size=(n, 3))
    return v / np.linalg.norm(v, axis=1, keepdims=True)


def check_a(B0=0.05, k=1e6, tensors=H16, n=400):
    m = isotropic_version(tensors)
    y = m.yield_eig(rand_dirs(n, 1), B0, k)
    spread = float(y.max() - y.min())
    return {"spread": spread, "pass": spread < TOL_A}


def check_b(tensors=H16, k=1e6):
    m = primary_model(tensors)
    y = yield_from_hamiltonian(np.zeros((m.dim, m.dim), complex), m.PS, m.M, k)
    return {"yield": y, "dev": abs(y - 1), "pass": abs(y - 1) < TOL_B}


def check_c(B0=0.05, k=1e6, tensors=H16, n=2000):
    y = primary_model(tensors).yield_eig(fibonacci_sphere(n), B0, k)
    return {"min": float(y.min()), "max": float(y.max()), "pass": bool(y.min() >= 0 and y.max() <= 1)}


def check_d(B0=0.05, k=1e6, tensors=H16, n=500):
    m = primary_model(tensors)
    v = rand_dirs(n, 2)
    d = float(np.max(np.abs(m.yield_eig(v, B0, k) - m.yield_eig(-v, B0, k))))
    # also sign flip of B0 (same as flipping n)
    return {"max_diff": d, "pass": d < TOL_D}


def check_e(B0=0.05, k=1e6, tensors=H16, n=200):
    m = primary_model(tensors)
    v = rand_dirs(n, 3)
    d = float(np.max(np.abs(m.yield_eig(v, B0, k) - m.yield_liouville(v, B0, k))))
    return {"max_diff": d, "pass": d < TOL_E}


def check_f(B0=0.05, k=1e6, tensors=H16):
    m = primary_model(tensors)
    f = lambda v: m.yield_eig(v, B0, k)
    a2 = sphere_extrema(f, 2000)["aniso"]
    a8 = sphere_extrema(f, 8000)["aniso"]
    sur, table = choose_degree(f)
    sh_err = None if sur is None else dict(table)[sur.L]
    return {
        "aniso_2000": a2,
        "aniso_8000": a8,
        "grid_change": abs(a2 - a8),
        "sh_degree": None if sur is None else sur.L,
        "sh_holdout_err": sh_err,
        "sh_table": table,
        "pass": bool(abs(a2 - a8) < TOL_F_GRID and sur is not None),
    }
