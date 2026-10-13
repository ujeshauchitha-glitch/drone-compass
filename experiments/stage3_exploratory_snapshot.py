"""EXPLORATORY (A2.3b): global sign/180-degree ambiguity of an axis-only magnetic sensor in a snapshot yaw initialisation.
TRIAD-style: tilt from the accelerometer; heading from the magnetic vector (arm B) or from the axis with polarity chosen so the
vertical field component points down (arm C).  Fraction of draws with |yaw error| > 90 deg."""
from __future__ import annotations

import csv
import sys
from itertools import product
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from nav import attitude as A  # noqa: E402

RES = ROOT / "results"
N = 20000


def triad(f_body, m_body, g_nav_up, m_nav):
    """Rotation body->nav from (up, m) pairs (TRIAD)."""
    def frame(a, b):
        t1 = a / np.linalg.norm(a)
        t2 = np.cross(a, b); t2 /= np.linalg.norm(t2)
        return np.stack([t1, t2, np.cross(t1, t2)], axis=1)
    return frame(g_nav_up, m_nav) @ frame(f_body, m_body).T


def main():
    rng = np.random.default_rng(np.random.SeedSequence(A.SEED, spawn_key=(12,)))
    sig_acc = float(np.hypot(0.05, 0.3))
    sig_mag = 0.5
    rows = []
    for dip, srp, sd in product((15.0, 45.0, 66.0), (1.0, 5.0, 10.0), (0.0, 2.0, 10.0)):
        be = A.earth_field(dip); bn = np.linalg.norm(be); un = be / bn
        up = np.array([0.0, 0.0, 1.0])
        fail_B = fail_C = 0
        for _ in range(N):
            yaw = rng.uniform(-np.pi, np.pi)
            Rt = A.Rz(np.array(yaw)) @ A.Ry(np.array(rng.normal(0, 0.05))) @ A.Rx(np.array(rng.normal(0, 0.05)))
            f = Rt.T @ (A.G * up) + sig_acc * rng.standard_normal(3)
            dist = sd / np.sqrt(3) * rng.standard_normal(3) + 0.0
            btrue = Rt.T @ be + dist
            # arm B
            mB = btrue + sig_mag * rng.standard_normal(3)
            Rb = triad(f, mB, up, be)
            eb = (np.degrees(np.arctan2((Rt @ Rb.T)[1, 0], (Rt @ Rb.T)[0, 0])))
            fail_B += abs(eb) > 90
            # arm C: axis only, tangent noise, random sign, polarity chosen so the vertical component points down
            u = btrue / np.linalg.norm(btrue)
            a = np.array([1.0, 0, 0]) if abs(u[0]) < 0.9 else np.array([0, 1.0, 0])
            e1 = np.cross(u, a); e1 /= np.linalg.norm(e1); e2 = np.cross(u, e1)
            um = u + np.deg2rad(srp) * (rng.standard_normal() * e1 + rng.standard_normal() * e2)
            um *= rng.choice([-1.0, 1.0]) / np.linalg.norm(um)
            fh = f / np.linalg.norm(f)
            um = um if (um @ fh) < 0 else -um           # field points down: component along 'up' negative
            Rc = triad(f, um, up, un)
            ec = np.degrees(np.arctan2((Rt @ Rc.T)[1, 0], (Rt @ Rc.T)[0, 0]))
            fail_C += abs(ec) > 90
        rows.append(dict(dip=dip, sigma_rp_deg=srp, sigma_d=sd, frac_yaw_gt_90deg_vector=fail_B / N, frac_yaw_gt_90deg_axis_only=fail_C / N))
        print(rows[-1], flush=True)
    with open(RES / "stage3_exploratory_snapshot.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)


if __name__ == "__main__":
    main()
