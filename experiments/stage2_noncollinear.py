"""Amendment 1, A1.5: does the mirror ambiguity depend on all tensors sharing one principal frame?
All orientations are ASSUMED; principal values of N5/N10 are the verified H16 values."""
from __future__ import annotations

import csv
import json
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from model.sphere import fibonacci_sphere, sphere_extrema  # noqa: E402
from model.spin import H16, RPModel  # noqa: E402
from sensor.design import N_DESIGN, X_DEG, mirror_dmin  # noqa: E402
from sensor.poisson import ELEMENT_SETS, rotvec  # noqa: E402

RES = ROOT / "results"
B0, K_RATE = 0.05, 1e6
FLIPS = [np.array([-1, -1, 1.0]), np.array([1, -1, -1.0]), np.array([-1, 1, -1.0])]


def rot(axis, deg):
    a = np.asarray(axis, float)
    return rotvec(a / np.linalg.norm(a) * np.deg2rad(deg))


def axial(a_iso, d):
    return np.diag([a_iso - d, a_iso - d, a_iso + 2 * d])


def make(kind, alpha):
    R1 = rot([1, 1, 0], alpha)
    A5 = np.diag(H16["N5"])
    A10 = R1 @ np.diag(H16["N10"]) @ R1.T
    if kind == "NC":
        return RPModel([(1, A5), (1, A10)])
    R2 = rot([0, 1, 1], alpha)
    return RPModel([(1, A5), (1, A10), (0.5, axial(0.8, 0.1))], [(0.5, R2 @ axial(0.5, 0.5) @ R2.T)])


def main():
    t0 = time.time()
    rng = np.random.default_rng(np.random.SeedSequence(20261013, spawn_key=(6, 0)))
    v = rng.normal(size=(500, 3)); v /= np.linalg.norm(v, axis=1, keepdims=True)
    dirs = v[:400]
    rows = []
    for kind, alphas in (("NC", (0, 5, 15, 30, 60)), ("NC-R2", (0, 15, 30, 60))):
        for a in alphas:
            m = make(kind, a)
            f = lambda x: m.yield_eig(x, B0, K_RATE)
            base = f(v)
            viol = max(float(np.max(np.abs(f(v * fl) - base))) for fl in FLIPS)
            dmin, per = mirror_dmin(f, ELEMENT_SETS["K3orth"], dirs, X_DEG, N_DESIGN)
            n_mirror = np.inf if dmin < 1e-9 else 10.0 * N_DESIGN / dmin
            ext = sphere_extrema(f, 600)
            rows.append(dict(model=kind, alpha_deg=a, sym_violation=viol, Dmin_nats_at_1e8=dmin,
                             Dmin_per_flip=[round(x, 6) for x in per], N_mirror_photons=n_mirror,
                             ambiguity_survives=bool(dmin < 10.0), aniso=ext["aniso"]))
            print(kind, a, f"viol={viol:.2e} Dmin={dmin:.3g} N_mirror={n_mirror:.3g} survives={dmin < 10} aniso={ext['aniso']:.4f}",
                  f"{time.time() - t0:.0f}s", flush=True)
    json.dump(rows, open(RES / "stage2_noncollinear.json", "w"), indent=1)
    with open(RES / "stage2_noncollinear.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)


if __name__ == "__main__":
    main()
