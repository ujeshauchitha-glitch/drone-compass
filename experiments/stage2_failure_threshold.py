"""EXPLORATORY addition (CHANGELOG): smallest N with ambiguity failure rate <= 5 % for the selected set."""
from __future__ import annotations

import csv
import json
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from experiments.stage2_sprime import FAIL_DEG, mc, sets  # noqa: E402
from model.spin import primary_model  # noqa: E402
from model.surrogate import PolySurrogate  # noqa: E402
from sensor.poisson import Sensor, hemisphere_directions  # noqa: E402

RES = ROOT / "results"
SEED = 20261013
B0, K_RATE = 0.05, 1e6
N_GRID = np.logspace(2, 12, 21)


def main():
    t0 = time.time()
    design = json.load(open(RES / "stage2_design.json"))
    meta = json.load(open(RES / "stage2_requirements_meta.json"))
    phi = PolySurrogate.fit(lambda v: primary_model().yield_eig(v, B0, K_RATE), 8)
    rots = sets(design)["K3selected"]
    eps_list = [("1", 1.0), ("0.1", 0.1)] + [(k, v) for k, v in meta["eps_derived"].items()]
    dirs = hemisphere_directions(300, np.random.default_rng(np.random.SeedSequence(SEED, spawn_key=(8, 0))))
    rows, summary = [], []
    for ei, (en, eps) in enumerate(eps_list):
        sen = Sensor(phi, rots, eps=eps)
        sen.make_grid(4000)
        fr = []
        for i, N in enumerate(N_GRID):
            rng = np.random.default_rng(np.random.SeedSequence(SEED, spawn_key=(8, 1, ei, i)))
            err = mc(sen, dirs, N, rng)
            fr.append(float((err > FAIL_DEG).mean()))
            rows.append(dict(eps_name=en, eps=eps, N=float(N), fail_rate_gt_10deg=fr[-1], mc_median_deg=float(np.median(err)),
                             mc_p90_deg=float(np.quantile(err, .9))))
        fr = np.array(fr)
        ok = [i for i in range(len(N_GRID)) if np.all(fr[i:] <= 0.05)]
        summary.append(dict(eps_name=en, eps=eps, N_fail_le_5pct=float(N_GRID[ok[0]]) if ok else None,
                            N_fail_le_5pct_rule="smallest grid N from which the failure rate stays <= 0.05 up to 1e12 (300 trials/point)"))
        print(en, eps, summary[-1]["N_fail_le_5pct"], f"{time.time() - t0:.0f}s", flush=True)
    with open(RES / "stage2_failure_curve.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
    json.dump(summary, open(RES / "stage2_failure_threshold.json", "w"), indent=1)


if __name__ == "__main__":
    main()
