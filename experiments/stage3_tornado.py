"""Stage 4 sensitivity (tornado), Amendment 3 A3.4.  One-at-a-time variation of every ASSUMED parameter around the base;
arms B and C (5 sigma_rp) for BOTH the pre-registered (uncalibrated) and the calibrated filter.  <= 2 workers."""
from __future__ import annotations

import csv
import json
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from nav import attitude as A  # noqa: E402
from nav.attitude_cal import run_filter_cal  # noqa: E402

RES = ROOT / "results"
M = 200
SIG_RP = (0.3, 1.0, 2.0, 5.0, 10.0)
BASE = dict(bias_dph=10.0, arw=0.3, sigma_d=2.0, dip=45.0, sigma_mag=0.5, soft_iron_scale=0.01, a_scale=1.0, vib_sigma=0.3, init_yaw_deg=10.0, kappa=1.0)
VARIATIONS = [("gyro bias instability (deg/h)", "bias_dph", (1.0, 100.0)), ("angle random walk (deg/sqrt(h))", "arw", (0.1, 1.0)),
              ("disturbance sigma_d (uT)", "sigma_d", (0.0, 10.0)), ("dip angle (deg)", "dip", (15.0, 66.0)),
              ("MEMS noise sigma_mag (uT)", "sigma_mag", (0.1, 2.0)), ("soft-iron std", "soft_iron_scale", (0.0, 0.02)),
              ("linear-acceleration scale", "a_scale", (0.0, 2.0)), ("accelerometer vibration sigma (m/s^2)", "vib_sigma", (0.1, 1.0)),
              ("initial yaw sigma (deg)", "init_yaw_deg", (2.0, 30.0)), ("sensor-specific disturbance factor kappa", "kappa", (0.5,))]
_RND = None
_SC = {}


def _init():
    global _RND
    _RND = A.Randoms(M)


def scen(a):
    if a not in _SC:
        _SC[a] = A.Scenario(a_scale=a)
    return _SC[a]


def points():
    pts = [("base", None, None, dict(BASE))]
    for name, key, vals in VARIATIONS:
        for v in vals:
            p = dict(BASE); p[key] = v
            pts.append((name, key, v, p))
    return pts


def run_task(t):
    pid, filt, arm, srp, p = t
    cfg = A.default_cfg(bias_dph=p["bias_dph"], arw=p["arw"], sigma_d=p["sigma_d"], dip=p["dip"], sigma_mag=p["sigma_mag"],
                        soft_iron=True, soft_iron_scale=p["soft_iron_scale"], vib_sigma=p["vib_sigma"], init_yaw_deg=p["init_yaw_deg"],
                        kappa=p["kappa"] if arm == "C" else 1.0, sigma_rp_deg=srp if srp else 0.0)
    fn = run_filter_cal if filt == "cal" else A.run_filter
    e = fn(arm, cfg, scen(p["a_scale"]), _RND)["yaw_err"]
    return dict(pid=pid, filt=filt, arm=arm, srp=srp, e60=e[:, 600].tolist(), e600=e[:, 6000].tolist())


def main():
    t0 = time.time()
    pts = points()
    tasks = []
    for pid, (name, key, v, p) in enumerate(pts):
        for filt in ("uncal", "cal"):
            tasks.append((pid, filt, "B", None, p))
            for srp in SIG_RP:
                tasks.append((pid, filt, "C", srp, p))
    print(f"{len(pts)} points, {len(tasks)} filter runs; estimated ~{len(tasks) * 15 / 2 / 60:.0f} min with 2 workers", flush=True)
    res = {}
    with ProcessPoolExecutor(max_workers=2, initializer=_init) as ex:
        for i, r in enumerate(ex.map(run_task, tasks, chunksize=1)):
            res[(r["pid"], r["filt"], r["arm"], r["srp"])] = r
            if (i + 1) % 20 == 0:
                print(f"{i + 1}/{len(tasks)} done, {time.time() - t0:.0f}s", flush=True)
    wrap = lambda x: (x + np.pi) % (2 * np.pi) - np.pi
    rm = lambda e: float(np.sqrt(np.mean(wrap(np.array(e)) ** 2)))
    rows = []
    for pid, (name, key, v, p) in enumerate(pts):
        for filt in ("uncal", "cal"):
            b = res[(pid, filt, "B", None)]
            reds = []
            for srp in SIG_RP:
                c = res[(pid, filt, "C", srp)]
                r60 = (rm(b["e60"]) - rm(c["e60"])) / rm(b["e60"]); r600 = (rm(b["e600"]) - rm(c["e600"])) / rm(b["e600"])
                reds.append((srp, r60, r600))
            best_min = max(min(r60, r600) for _, r60, r600 in reds)
            best_600 = max(r600 for _, _, r600 in reds)
            rows.append(dict(parameter=name, value=v if v is not None else "base", filter=filt, metric_best_min60_600=best_min,
                             metric_best_reduction_600=best_600, rmse_B_60_deg=np.degrees(rm(b["e60"])), rmse_B_600_deg=np.degrees(rm(b["e600"])),
                             reductions_by_sigma_rp=json.dumps([(s, round(a, 4), round(c, 4)) for s, a, c in reds])))
    with open(RES / "stage3_tornado.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
    print("done", f"{time.time() - t0:.0f}s")


if __name__ == "__main__":
    main()
