"""Gate S2' (Amendment 1, A1.3-A1.4): evaluate the selected element set on FRESH seeds and FRESH directions.
Also: ambiguity failure rates (A1.4) and the Ren et al. comparison (condition 3)."""
from __future__ import annotations

import csv
import json
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from model.sphere import fibonacci_hemisphere  # noqa: E402
from model.spin import primary_model  # noqa: E402
from model.surrogate import PolySurrogate  # noqa: E402
from sensor.design import build_rots  # noqa: E402
from sensor.poisson import ELEMENT_SETS, RenPhi, Sensor, axis_error_deg, hemisphere_directions, rotvec  # noqa: E402

RES = ROOT / "results"
SEED = 20261013
N_GRID = np.logspace(5, 11, 13)
N_TRIALS = 300
B0, K_RATE = 0.05, 1e6
FAIL_DEG = 10.0


def sets(design):
    gen = [np.eye(3), rotvec(np.array([0.3, 0.5, 0.2])), rotvec(np.array([-0.4, 0.2, 0.6]))]
    return {"K3orth": ELEMENT_SETS["K3orth"], "K3generic_EXPLORATORY": gen,
            "K3selected": build_rots(design["selected"]["3"]["params"], 3)}


def mc(sen, dirs, N, rng):
    return np.array([axis_error_deg(sen.mle(sen.simulate(n, N, rng), N), n) for n in dirs])


def gate_from_rows(rows):
    crb_med = np.array([r["crb_median_of_sigmed"] for r in rows])
    ok = np.where(crb_med <= 1.0)[0]
    big = ok[-3:]
    ratios = [rows[j]["ratio"] for j in big]
    g1 = bool(len(big) == 3 and all(0.5 <= r <= 2.0 for r in ratios))
    sel = [j for j in range(len(rows)) if 0.1 <= crb_med[j] <= 3.0]
    slope = None
    if len(sel) >= 3:
        x = np.log10([rows[j]["N"] for j in sel]); y = np.log10([rows[j]["mc_median"] for j in sel])
        slope = float(np.polyfit(x, y, 1)[0])
    g2 = bool(slope is not None and abs(slope + 0.5) <= 0.1)
    return g1, g2, [rows[j]["N"] for j in big], ratios, slope, len(sel)


def main():
    t0 = time.time()
    design = json.load(open(RES / "stage2_design.json"))
    m = primary_model()
    phi = PolySurrogate.fit(lambda v: m.yield_eig(v, B0, K_RATE), 8)
    S = sets(design)

    # ---------------- Gate S2' : selected set, eps = 1, Poisson only, fresh directions and seeds (prefix 4)
    sen = Sensor(phi, S["K3selected"], eps=1.0)
    sen.make_grid(4000)
    dirs = hemisphere_directions(N_TRIALS, np.random.default_rng(np.random.SeedSequence(SEED, spawn_key=(4, 0))))
    rows = []
    for i, N in enumerate(N_GRID):
        rng = np.random.default_rng(np.random.SeedSequence(SEED, spawn_key=(4, 1, i)))
        err = mc(sen, dirs, N, rng)
        _, smed, sing = sen.crb(dirs, N)
        sd = np.degrees(smed)
        rows.append(dict(N=float(N), mc_median=float(np.median(err)), mc_mean=float(err.mean()), mc_p90=float(np.quantile(err, .9)),
                         crb_median_of_sigmed=float(np.median(sd)), ratio=float(np.median(err) / np.median(sd)),
                         frac_err_gt_10deg=float((err > FAIL_DEG).mean()), n_singular=int(sing.sum())))
    g1, g2, bigN, ratios, slope, npts = gate_from_rows(rows)
    out = {"gate_S2prime": dict(set="K3selected", eps=1.0, rows=rows, large_N=bigN, large_N_ratios=ratios, cond1_pass=g1,
                                slope=slope, slope_points=npts, cond2_pass=g2)}
    print("S2' cond1", g1, np.round(ratios, 3), "cond2", g2, slope, f"{time.time() - t0:.0f}s", flush=True)

    # ---------------- ambiguity failure rates (prefix 5), three sets x eps in {1, 0.1}
    fdirs = hemisphere_directions(N_TRIALS, np.random.default_rng(np.random.SeedSequence(SEED, spawn_key=(5, 0))))
    fail_rows = []
    for si, (name, rots) in enumerate(S.items()):
        for ei, eps in enumerate((1.0, 0.1)):
            s2 = Sensor(phi, rots, eps=eps)
            s2.make_grid(4000)
            for i, N in enumerate(N_GRID):
                rng = np.random.default_rng(np.random.SeedSequence(SEED, spawn_key=(5, 1, si, ei, i)))
                err = mc(s2, fdirs, N, rng)
                fail_rows.append(dict(set=name, eps=eps, N=float(N), fail_rate_gt_10deg=float((err > FAIL_DEG).mean()),
                                      mc_median=float(np.median(err)), mc_p90=float(np.quantile(err, .9))))
        print("failure rates done:", name, f"{time.time() - t0:.0f}s", flush=True)
    with open(RES / "stage2_failure_rates.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(fail_rows[0].keys())); w.writeheader(); w.writerows(fail_rows)
    out["failure_rates_file"] = "results/stage2_failure_rates.csv"

    # ---------------- Ren et al. comparison (condition 3): CRB only, selected set
    grid = fibonacci_hemisphere(2000)
    ren = {}
    refs = {"default_a0.5_b0.00127_c26": (RenPhi(0.5, 0.00127, 26.0), {"O1": 2.29e8, "O2_text": 5.26e7}),
            "optimised_a0.9_b0.00508_c10": (RenPhi(0.9, 0.00508, 10.0), {"O2_table2": 3.94e5})}
    for nm, (rp, theirs) in refs.items():
        s3 = Sensor(rp, S["K3selected"], eps=1.0)
        _, smed, sing = s3.crb(grid, 1e8)
        sd = np.degrees(smed[np.isfinite(smed)])
        med, worst = float(np.median(sd)), float(sd.max())
        n_med5 = 1e8 * (med / 5.0) ** 2
        n_worst5 = 1e8 * (worst / 5.0) ** 2
        ren[nm] = dict(sigma_med_deg_at_1e8_median_over_directions=med, sigma_med_deg_at_1e8_worst=worst,
                       N_per_element_for_5deg_median=n_med5, N_total_K3_median=3 * n_med5, N_per_element_worst=n_worst5,
                       N_total_K3_worst=3 * n_worst5, ren_reported=theirs, singular_directions=int(sing.sum()))
    out["ren_comparison"] = ren
    out["ren_differences"] = [
        "Ren: information-theoretic lower bound on RMS head-direction error (data processing inequality); ours: Cramer-Rao bound on a 2-D field axis (median-error form).",
        "Ren: r receptor groups on a hemispherical retina with receptor orientations O1/O2, head rotation in the horizontal plane (1-D unknown); ours: K=3 elements, unknown 2-D axis.",
        "Ren: uniform prior on head direction and optimal decoder; ours: local frequentist bound, per direction.",
        "Ren: noise variance Phi(1-Phi)/p per group; ours: Poisson counts with relative-contrast ratio eps.",
        "Ren: N = r*p total absorbed photons over all groups; ours: N per element, total = 3N.",
        "Signal model in this row is Ren's phenomenological Gaussian dip (exactly axial); ours elsewhere is the spin-model yield.",
    ]
    out["runtime_s"] = time.time() - t0
    json.dump(out, open(RES / "stage2_sprime.json", "w"), indent=1)
    print("done", f"{out['runtime_s']:.0f}s")


if __name__ == "__main__":
    main()
