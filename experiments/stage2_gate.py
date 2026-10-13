"""Gate S2 (PLAN.md 4.7) evaluated exactly as pre-registered on K3orth, plus a clearly labelled EXPLORATORY
diagnostic on a fixed generic (non-orthogonal) K=3 set.  Writes results/stage2_gate_*.json."""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from model.spin import primary_model  # noqa: E402
from model.surrogate import PolySurrogate  # noqa: E402
from sensor.poisson import ELEMENT_SETS, Sensor, axis_error_deg, hemisphere_directions, rotvec  # noqa: E402

RES = ROOT / "results"
SEED = 20261013
N_GRID = np.logspace(5, 11, 13)
N_TRIALS = 300
B0, K_RATE = 0.05, 1e6

ELEMENT_SETS["K3generic_EXPLORATORY"] = [np.eye(3), rotvec(np.array([0.3, 0.5, 0.2])), rotvec(np.array([-0.4, 0.2, 0.6]))]


def run_config(name, phi, rots, eps=1.0):
    sen = Sensor(phi, rots, eps=eps)
    sen.make_grid(4000)
    dirs = hemisphere_directions(N_TRIALS, np.random.default_rng(np.random.SeedSequence(SEED, spawn_key=(2, 0))))
    rows = []
    for i, N in enumerate(N_GRID):
        rng = np.random.default_rng(np.random.SeedSequence(SEED, spawn_key=(2, 1, i)))
        err = np.array([axis_error_deg(sen.mle(sen.simulate(n, N, rng), N), n) for n in dirs])
        _, smed, sing = sen.crb(dirs, N)
        smed_deg = np.degrees(smed)
        rows.append(dict(N=float(N), mc_median=float(np.median(err)), mc_mean=float(err.mean()),
                         mc_p90=float(np.quantile(err, 0.9)), crb_median_of_sigmed=float(np.median(smed_deg)),
                         ratio=float(np.median(err) / np.median(smed_deg)),
                         frac_err_gt_5deg=float((err > 5).mean()), n_singular=int(sing.sum())))
    # gate evaluation
    crb_med = np.array([r["crb_median_of_sigmed"] for r in rows])
    ok_idx = np.where(crb_med <= 1.0)[0]
    big = ok_idx[-3:]
    ratios = [rows[j]["ratio"] for j in big]
    g1 = bool(len(big) == 3 and all(0.5 <= r <= 2.0 for r in ratios))
    sel = [j for j in range(len(rows)) if 0.1 <= crb_med[j] <= 3.0]
    slope = None
    if len(sel) >= 3:
        x = np.log10([rows[j]["N"] for j in sel]); y = np.log10([rows[j]["mc_median"] for j in sel])
        slope = float(np.polyfit(x, y, 1)[0])
    g2 = bool(slope is not None and abs(slope + 0.5) <= 0.1)
    return dict(name=name, eps=eps, rows=rows, large_N=[rows[j]["N"] for j in big], large_N_ratios=ratios,
                gate_S2_1_pass=g1, slope=slope, slope_points=len(sel), gate_S2_2_pass=g2)


def main():
    t0 = time.time()
    m = primary_model()
    phi = PolySurrogate.fit(lambda v: m.yield_eig(v, B0, K_RATE), 8)
    out = {}
    for name in ("K3orth", "K3generic_EXPLORATORY"):
        out[name] = run_config(name, phi, ELEMENT_SETS[name])
        o = out[name]
        print(name, "gate1", o["gate_S2_1_pass"], "ratios", np.round(o["large_N_ratios"], 3), "slope", o["slope"],
              "gate2", o["gate_S2_2_pass"], f"{time.time() - t0:.0f}s", flush=True)
    out["note"] = ("Gate S2 conditions 1-2 are pre-registered for K3orth only.  K3generic_EXPLORATORY is a diagnostic "
                   "added after K3orth failed and is NOT a gate result.")
    json.dump(out, open(RES / "stage2_gate.json", "w"), indent=1)


if __name__ == "__main__":
    main()
