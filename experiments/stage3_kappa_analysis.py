"""EXPLORATORY kappa analysis with the corrected model (CHANGELOG): supersedes the earlier kappa table."""
from __future__ import annotations

import csv
import gzip
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
RES = ROOT / "results"
M, N_BOOT = 200, 2000


def fnum(x):
    try:
        return float(x)
    except ValueError:
        return float("nan")


def read(sfile, efile):
    rows = list(csv.DictReader(open(RES / sfile)))
    for r in rows:
        for k in list(r):
            if k not in ("group", "arm"):
                r[k] = fnum(r[k])
        r["id"] = int(r["id"])
    errs = defaultdict(dict)
    with gzip.open(RES / efile, "rt") as fh:
        rd = csv.reader(fh); next(rd)
        for cid, t, run, v in rd:
            errs[int(cid)].setdefault(int(t), []).append(float(v))
    return rows, {k: {t: np.array(v) for t, v in d.items()} for k, d in errs.items()}


def main():
    krows, kerrs = read("stage3_summary_kappa_EXPLORATORY.csv", "stage3_yaw_errors_kappa_EXPLORATORY.csv.gz")
    mrows, merrs = read("stage3_summary.csv", "stage3_yaw_errors.csv.gz")
    B = {(r["bias_dph"], r["arw"], r["sigma_d"], r["dip"], r["sigma_mag"]): r for r in mrows if r["arm"] == "B"}
    rng = np.random.default_rng(np.random.SeedSequence(20261013, spawn_key=(11,)))
    idx = rng.integers(0, M, (N_BOOT, M))
    wrap = lambda x: (x + np.pi) % (2 * np.pi) - np.pi

    def red(eb, ex):
        eb, ex = wrap(eb), wrap(ex)
        bb, bx = np.sqrt(np.mean(eb[idx] ** 2, 1)), np.sqrt(np.mean(ex[idx] ** 2, 1))
        r = (bb - bx) / bb
        return float((np.sqrt(np.mean(eb**2)) - np.sqrt(np.mean(ex**2))) / np.sqrt(np.mean(eb**2))), float(np.percentile(r, 2.5))

    out = []
    for r in krows:
        b = B[(r["bias_dph"], r["arw"], r["sigma_d"], r["dip"], r["sigma_mag"])]
        r60, c60 = red(merrs[b["id"]][60], kerrs[r["id"]][60]); r600, c600 = red(merrs[b["id"]][600], kerrs[r["id"]][600])
        out.append(dict(arm=r["arm"], kappa=r["kappa"], sigma_d=r["sigma_d"], sigma_rp_deg=r["sigma_rp_deg"], reduction_60=r60, ci_lo_60=c60,
                        reduction_600=r600, ci_lo_600=c600, criterion_met=bool(r60 >= 0.2 and c60 > 0 and r600 >= 0.2 and c600 > 0)))
    with open(RES / "stage3_kappa_exploratory_corrected.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(out[0].keys())); w.writeheader(); w.writerows(out)
    for kappa in (0.5, 0.1, 0.0):
        for sd in (2.0, 10.0):
            sub = [o for o in out if o["kappa"] == kappa and o["sigma_d"] == sd]
            print(f"kappa={kappa} sigma_d={sd}: criterion met in {sum(o['criterion_met'] for o in sub)} of {len(sub)}",
                  [(o['arm'], o['sigma_rp_deg'], round(o['reduction_60'], 3), round(o['reduction_600'], 3)) for o in sub])


if __name__ == "__main__":
    main()
