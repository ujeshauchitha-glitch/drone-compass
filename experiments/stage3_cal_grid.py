"""Calibrated decision study (Amendment 3, A3.3): 27 settings x (B + 5 C + 5 D) = 297 configurations, <= 2 workers."""
from __future__ import annotations

import csv
import gzip
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from itertools import product
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from nav import attitude as A  # noqa: E402
from nav.attitude_cal import run_filter_cal  # noqa: E402

RES = ROOT / "results"
M = 200
SIG_RP = (0.3, 1.0, 2.0, 5.0, 10.0)
SIG_D = (0.0, 2.0, 10.0)
MID = (10.0, 0.3)
_SC = _RND = None


def _init():
    global _SC, _RND
    _SC, _RND = A.Scenario(), A.Randoms(M)


def settings():
    s = []
    for smag in (0.1, 0.5, 2.0):
        s.append(dict(setting="S1", gyro=MID, dip=45.0, smag=smag))
    for dip in (15.0, 66.0):
        s.append(dict(setting="S2", gyro=MID, dip=dip, smag=0.5))
    for g in ((1.0, 0.1), (1.0, 1.0), (100.0, 0.1), (100.0, 1.0)):
        s.append(dict(setting="S3", gyro=g, dip=45.0, smag=0.5))
    return s


def enumerate_configs():
    cfgs = []
    for st, sd in product(settings(), SIG_D):
        base = dict(setting=st["setting"], bias_dph=st["gyro"][0], arw=st["gyro"][1], sigma_d=sd, dip=st["dip"], sigma_mag=st["smag"])
        cfgs.append(dict(base, arm="B", sigma_rp_deg=float("nan")))
        for srp in SIG_RP:
            cfgs.append(dict(base, arm="C", sigma_rp_deg=srp))
            cfgs.append(dict(base, arm="D", sigma_rp_deg=srp))
    for i, c in enumerate(cfgs):
        c["id"] = i
    return cfgs


def run_one(c):
    cfg = A.default_cfg(bias_dph=c["bias_dph"], arw=c["arw"], sigma_d=c["sigma_d"], dip=c["dip"], sigma_mag=c["sigma_mag"],
                        sigma_rp_deg=0.0 if np.isnan(c["sigma_rp_deg"]) else c["sigma_rp_deg"], soft_iron=True)
    r = run_filter_cal(c["arm"], cfg, _SC, _RND, nees_steps=(600, 3000, 6000))
    e = r["yaw_err"]
    return dict(id=c["id"], e60=e[:, 600].tolist(), e600=e[:, 6000].tolist(),
                rmse_60=A.yaw_rmse(e, 60.0), rmse_300=A.yaw_rmse(e, 300.0), rmse_600=A.yaw_rmse(e, 600.0),
                nees_60=float(r["nees"][600].mean()), nees_300=float(r["nees"][3000].mean()), nees_600=float(r["nees"][6000].mean()))


def main():
    t0 = time.time()
    cfgs = enumerate_configs()
    print(f"{len(cfgs)} configurations; estimated runtime ~{len(cfgs) * 27 / 2 / 60:.0f} min with 2 workers", flush=True)
    res = {}
    with ProcessPoolExecutor(max_workers=2, initializer=_init) as ex:
        for i, r in enumerate(ex.map(run_one, cfgs, chunksize=1)):
            res[r["id"]] = r
            if (i + 1) % 20 == 0:
                print(f"{i + 1}/{len(cfgs)} done, {time.time() - t0:.0f}s", flush=True)
    rows = []
    for c in cfgs:
        row = dict(c)
        row.update({k: res[c["id"]][k] for k in ("rmse_60", "rmse_300", "rmse_600", "nees_60", "nees_300", "nees_600")})
        rows.append(row)
    with open(RES / "stage3_cal_summary.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
    with gzip.open(RES / "stage3_cal_yaw_errors.csv.gz", "wt", newline="") as fh:
        w = csv.writer(fh); w.writerow(["id", "t_s", "run", "yaw_err_rad"])
        for c in cfgs:
            for t_s, key in ((60, "e60"), (600, "e600")):
                for run, v in enumerate(res[c["id"]][key]):
                    w.writerow([c["id"], t_s, run, f"{v:.9g}"])
    print("done", f"{time.time() - t0:.0f}s")


if __name__ == "__main__":
    main()
