"""EXPLORATORY (CHANGELOG, logged before running): separate sensor noise from the unmodelled soft-iron systematic.
E1: filters whose R includes the soft-iron variance (si_R_um = 0.5 uT per axis) for both magnetic sensors.
E2: no soft iron at all (pre-registered R).  Mid gyro, dip 45 deg.  Same paired-bootstrap decision statistic as A2.3."""
from __future__ import annotations

import csv
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
SIG_RP = (0.3, 1.0, 2.0, 5.0, 10.0)
_SC = _RND = None


def _init():
    global _SC, _RND
    _SC, _RND = A.Scenario(), A.Randoms(M)


def configs():
    out = []
    for name, soft, si, smags in (("E1_tuned_R", True, 0.5, (0.1, 0.5, 2.0)), ("E2_no_soft_iron", False, 0.0, (0.5,))):
        for smag, sd in product(smags, (0.0, 2.0, 10.0)):
            out.append(dict(exp=name, arm="B", soft_iron=soft, si_R_um=si, sigma_mag=smag, sigma_d=sd, sigma_rp_deg=float("nan")))
            for srp in SIG_RP:
                out.append(dict(exp=name, arm="C", soft_iron=soft, si_R_um=si, sigma_mag=smag, sigma_d=sd, sigma_rp_deg=srp))
                out.append(dict(exp=name, arm="D", soft_iron=soft, si_R_um=si, sigma_mag=smag, sigma_d=sd, sigma_rp_deg=srp))
    for i, c in enumerate(out):
        c["id"] = i
    return out


def run_one(c):
    cfg = A.default_cfg(bias_dph=10.0, arw=0.3, sigma_d=c["sigma_d"], dip=45.0, sigma_mag=c["sigma_mag"],
                        sigma_rp_deg=0.0 if np.isnan(c["sigma_rp_deg"]) else c["sigma_rp_deg"], soft_iron=c["soft_iron"], si_R_um=c["si_R_um"])
    r = A.run_filter(c["arm"], cfg, _SC, _RND, nees_steps=(600, 3000, 6000))
    e = r["yaw_err"]
    return dict(id=c["id"], e60=e[:, 600], e600=e[:, 6000], nees=[float(r["nees"][k].mean()) for k in (600, 3000, 6000)])


def main():
    t0 = time.time()
    cfgs = configs()
    print(f"{len(cfgs)} configurations (~{len(cfgs) * 10 / 2 / 60:.0f} min with 2 workers)", flush=True)
    with ProcessPoolExecutor(max_workers=2, initializer=_init) as ex:
        res = {r["id"]: r for r in ex.map(run_one, cfgs, chunksize=1)}
    rng = np.random.default_rng(np.random.SeedSequence(20261013, spawn_key=(11,)))
    idx = rng.integers(0, M, (2000, M))
    base = {(c["exp"], c["sigma_mag"], c["sigma_d"]): c for c in cfgs if c["arm"] == "B"}
    rows = []
    wrap = lambda x: (x + np.pi) % (2 * np.pi) - np.pi

    def stats(eb, ex):
        eb, ex = wrap(eb), wrap(ex)
        rb, rx = np.sqrt(np.mean(eb**2)), np.sqrt(np.mean(ex**2))
        red = (np.sqrt(np.mean(eb[idx] ** 2, 1)) - np.sqrt(np.mean(ex[idx] ** 2, 1))) / np.sqrt(np.mean(eb[idx] ** 2, 1))
        return rb, rx, (rb - rx) / rb, float(np.percentile(red, 2.5))

    for c in cfgs:
        if c["arm"] == "B":
            continue
        b = base[(c["exp"], c["sigma_mag"], c["sigma_d"])]
        r, rb = res[c["id"]], res[b["id"]]
        s60, s600 = stats(rb["e60"], r["e60"]), stats(rb["e600"], r["e600"])
        rows.append(dict(exp=c["exp"], arm=c["arm"], sigma_mag=c["sigma_mag"], sigma_d=c["sigma_d"], sigma_rp_deg=c["sigma_rp_deg"],
                         rmse_B_60_deg=np.degrees(s60[0]), rmse_X_60_deg=np.degrees(s60[1]), reduction_60=s60[2], ci_lo_60=s60[3],
                         rmse_B_600_deg=np.degrees(s600[0]), rmse_X_600_deg=np.degrees(s600[1]), reduction_600=s600[2], ci_lo_600=s600[3],
                         criterion_met=bool(all(s[2] >= 0.2 and s[3] > 0 for s in (s60, s600))),
                         nees_B_600=res[b["id"]]["nees"][2], nees_X_600=r["nees"][2]))
    with open(RES / "stage3_exploratory_tuned.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
    nb = {(k[0], k[1], k[2]): res[v["id"]]["nees"] for k, v in base.items()}
    json.dump({f"{k[0]}|sigma_mag={k[1]}|sigma_d={k[2]}": v for k, v in nb.items()}, open(RES / "stage3_exploratory_tuned_nees_B.json", "w"), indent=1)
    for exp in ("E1_tuned_R", "E2_no_soft_iron"):
        sub = [r for r in rows if r["exp"] == exp]
        print(exp, "criterion met:", sum(r["criterion_met"] for r in sub), "of", len(sub))
    print("done", f"{time.time() - t0:.0f}s")


if __name__ == "__main__":
    main()
