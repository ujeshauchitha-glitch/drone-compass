"""Stage 3 main grid (PLAN.md 5.5, Amendment 2 A2.6) + EXPLORATORY kappa runs.  <= 2 worker processes.
Writes results/stage3_summary.csv and results/stage3_yaw_errors.csv.gz (yaw error at t = 60 s and 600 s, every run)."""
from __future__ import annotations

import csv
import gzip
import json
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from itertools import product
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from nav import attitude as A  # noqa: E402

RES = ROOT / "results"
M = 200
GYRO = [(b, a) for b in (1.0, 10.0, 100.0) for a in (0.1, 0.3, 1.0)]
SIG_D = (0.0, 2.0, 10.0)
SIG_RP = (0.3, 1.0, 2.0, 5.0, 10.0)
MID = (10.0, 0.3)

_SC = _RND = None


def _init():
    global _SC, _RND
    _SC, _RND = A.Scenario(), A.Randoms(M)


def enumerate_configs():
    cfgs = []

    def add(group, arm, gyro, sd, dip=45.0, smag=0.5, srp=2.0, kappa=1.0):
        cfgs.append(dict(group=group, arm=arm, bias_dph=gyro[0], arw=gyro[1], sigma_d=sd, dip=dip, sigma_mag=smag,
                         sigma_rp_deg=srp if arm in "CD" else float("nan"), kappa=kappa))

    for g in GYRO:
        add("base", "A", g, 0.0)
    for g, sd in product(GYRO, SIG_D):
        add("base", "B", g, sd)
        for srp in SIG_RP:
            add("base", "C", g, sd, srp=srp)
            add("base", "D", g, sd, srp=srp)
    for (dip, smag), sd in product(((15.0, 0.5), (66.0, 0.5), (45.0, 0.1), (45.0, 2.0)), SIG_D):
        add("oat", "B", MID, sd, dip=dip, smag=smag)
        for srp in SIG_RP:
            add("oat", "C", MID, sd, dip=dip, smag=smag, srp=srp)
            add("oat", "D", MID, sd, dip=dip, smag=smag, srp=srp)
    for kappa, sd, srp in product((0.5, 0.1, 0.0), (2.0, 10.0), (0.3, 1.0)):
        add("kappa_EXPLORATORY", "C", MID, sd, srp=srp, kappa=kappa)
        add("kappa_EXPLORATORY", "D", MID, sd, srp=srp, kappa=kappa)
    for i, c in enumerate(cfgs):
        c["id"] = i
    return cfgs


def run_one(c):
    cfg = A.default_cfg(bias_dph=c["bias_dph"], arw=c["arw"], sigma_d=c["sigma_d"], dip=c["dip"], sigma_mag=c["sigma_mag"],
                        sigma_rp_deg=0.0 if np.isnan(c["sigma_rp_deg"]) else c["sigma_rp_deg"], kappa=c["kappa"], soft_iron=True)
    r = A.run_filter(c["arm"], cfg, _SC, _RND, nees_steps=(600, 3000, 6000))
    e = r["yaw_err"]
    rm = {t: A.yaw_rmse(e, t) for t in (60.0, 300.0, 600.0)}
    ts = np.arange(0, 601, 10.0)
    series = [A.yaw_rmse(e, t) for t in ts]
    return dict(id=c["id"], rmse_60=rm[60.0], rmse_300=rm[300.0], rmse_600=rm[600.0],
                nees_60=float(r["nees"][600].mean()), nees_300=float(r["nees"][3000].mean()), nees_600=float(r["nees"][6000].mean()),
                tilt_rmse_final=float(np.sqrt(np.mean(r["tilt_err"] ** 2))), series=series,
                e60=e[:, 600].tolist(), e600=e[:, 6000].tolist())


def main(only_group=None):
    t0 = time.time()
    cfgs = enumerate_configs()
    if only_group:
        cfgs = [c for c in cfgs if c["group"] == only_group]
    print(f"{len(cfgs)} configurations; estimated runtime ~{len(cfgs) * 9.5 / 2 / 60:.0f} min with 2 workers", flush=True)
    results = {}
    with ProcessPoolExecutor(max_workers=2, initializer=_init) as ex:
        for i, res in enumerate(ex.map(run_one, cfgs, chunksize=1)):
            results[res["id"]] = res
            if (i + 1) % 20 == 0:
                print(f"{i + 1}/{len(cfgs)} done, {time.time() - t0:.0f}s", flush=True)
    rows = []
    for c in cfgs:
        r = results[c["id"]]
        row = dict(c)
        row.update({k: r[k] for k in ("rmse_60", "rmse_300", "rmse_600", "nees_60", "nees_300", "nees_600", "tilt_rmse_final")})
        rows.append(row)
    suffix = f"_{only_group}" if only_group else ""
    with open(RES / f"stage3_summary{suffix}.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
    with gzip.open(RES / f"stage3_yaw_errors{suffix}.csv.gz", "wt", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["id", "t_s", "run", "yaw_err_rad"])
        for c in cfgs:
            r = results[c["id"]]
            for t_s, key in ((60, "e60"), (600, "e600")):
                for run, v in enumerate(r[key]):
                    w.writerow([c["id"], t_s, run, f"{v:.9g}"])
    with open(RES / f"stage3_rmse_series{suffix}.json", "w") as fh:
        json.dump({str(c["id"]): results[c["id"]]["series"] for c in cfgs}, fh)
    print("done", f"{time.time() - t0:.0f}s")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else None)
