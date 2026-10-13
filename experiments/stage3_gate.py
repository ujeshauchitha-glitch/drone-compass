"""Gate S3 (PLAN.md 5.6 + Amendment 2 A2.5), evaluated exactly as pre-registered, plus labelled diagnostics."""
from __future__ import annotations

import json
import sys
import time
from itertools import product
from pathlib import Path

import numpy as np
from scipy.stats import chi2

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from nav import attitude as A  # noqa: E402
from nav.reference import ref_run  # noqa: E402

RES = ROOT / "results"
M = 200
STEPS = (600, 3000, 6000)
N_TESTS = 12


def nees_interval(alpha, dof=6, runs=M):
    return chi2.ppf(alpha / 2, dof * runs) / runs, chi2.ppf(1 - alpha / 2, dof * runs) / runs


def main():
    t0 = time.time()
    sc, rnd = A.Scenario(), A.Randoms(M)
    out = {"scenario": dict(max_yaw_rate_dps=sc.max_yaw_rate_dps, max_acc_mps2=sc.max_acc, max_tilt_deg=sc.max_tilt_deg,
                            yaw_at_300s_deg=float(np.degrees(sc.yaw_true[3000])), yaw_end_deg=float(np.degrees(sc.yaw_true[-1])))}
    lo95, hi95 = nees_interval(0.05)
    loB, hiB = nees_interval(0.05 / N_TESTS)
    out["nees_intervals"] = dict(uncorrected_95=[lo95, hi95], bonferroni_family_95=[loB, hiB], n_tests=N_TESTS)

    # ---------------- (i) NEES, nominal condition, mid gyro
    nom = A.default_cfg(soft_iron=False, sigma_d=0.0, sigma_rp_deg=2.0)
    rows = []
    sc_diag = A.Scenario(); sc_diag.fb = np.einsum("kji,j->ki", sc_diag.R, np.array([0, 0, A.G]))
    for label, scen in (("PRE-REGISTERED trajectory", sc), ("DIAGNOSTIC (no linear acceleration; not a gate)", sc_diag)):
        for arm in "ABCD":
            r = A.run_filter(arm, nom, scen, rnd, nees_steps=STEPS)
            for k in STEPS:
                v = float(r["nees"][k].mean())
                rows.append(dict(condition=label, arm=arm, t_s=k * A.DT, nees=v,
                                 in_uncorrected_95=bool(lo95 <= v <= hi95), in_bonferroni_95=bool(loB <= v <= hiB)))
            print(label[:12], arm, [round(x["nees"], 2) for x in rows[-3:]], f"{time.time() - t0:.0f}s", flush=True)
    out["nees_rows"] = rows
    pre = [r for r in rows if r["condition"].startswith("PRE")]
    out["gate_S3_i_pass_bonferroni"] = bool(all(r["in_bonferroni_95"] for r in pre))
    out["gate_S3_i_pass_uncorrected"] = bool(all(r["in_uncorrected_95"] for r in pre))
    out["gate_S3_i_failures"] = [(r["arm"], r["t_s"], round(r["nees"], 3)) for r in pre if not r["in_bonferroni_95"]]

    # NEES at the other eight gyro settings (reported, not gated), arms A and B
    other = []
    for bd, aw in product((1.0, 10.0, 100.0), (0.1, 0.3, 1.0)):
        for arm in "AB":
            r = A.run_filter(arm, A.default_cfg(soft_iron=False, bias_dph=bd, arw=aw), sc, rnd, nees_steps=STEPS)
            other.append(dict(bias_dph=bd, arw=aw, arm=arm, nees=[float(r["nees"][k].mean()) for k in STEPS]))
    out["nees_other_gyro_settings_not_gated"] = other
    print("other gyro settings done", f"{time.time() - t0:.0f}s", flush=True)

    # ---------------- (ii) arm A vs analytic gyro-drift expectation (nine gyro settings, t = 60 and 600 s)
    ii = []
    for bd, aw in product((1.0, 10.0, 100.0), (0.1, 0.3, 1.0)):
        cfg = A.default_cfg(soft_iron=False, bias_dph=bd, arw=aw)
        r = A.run_filter("A", cfg, sc, rnd)
        sb = bd / 3600 * A.DEG; arw = aw / 60 * A.DEG; sbw = 0.1 * sb / np.sqrt(600.0); s0 = 10 * A.DEG
        for t in (60.0, 600.0):
            exp_ = float(np.sqrt(sb**2 * t**2 + arw**2 * t + sbw**2 * t**3 / 3 + s0**2))
            sim = A.yaw_rmse(r["yaw_err"], t)
            ii.append(dict(bias_dph=bd, arw=aw, t_s=t, sim_rad=sim, analytic_rad=exp_, ratio=sim / exp_, within_factor_2=bool(0.5 <= sim / exp_ <= 2.0)))
    out["gate_S3_ii_rows"] = ii
    out["gate_S3_ii_pass"] = bool(all(x["within_factor_2"] for x in ii))
    print("gate ii", out["gate_S3_ii_pass"], "ratio range", round(min(x["ratio"] for x in ii), 3), round(max(x["ratio"] for x in ii), 3), f"{time.time() - t0:.0f}s", flush=True)

    # ---------------- (iii) independent non-vectorised filter on 5 runs, arms A-D
    iii = []
    cfg = A.default_cfg(soft_iron=True, sigma_d=2.0, sigma_rp_deg=2.0)   # a condition with all terms active
    for arm in "ABCD":
        batched = A.run_filter(arm, cfg, sc, rnd)["yaw_err"]
        for rr in range(5):
            ref = ref_run(arm, cfg, sc, rnd, rr)
            iii.append(dict(arm=arm, run=rr, max_abs_diff_rad=float(np.max(np.abs(ref - batched[rr])))))
        print("gate iii arm", arm, [round(x["max_abs_diff_rad"], 12) for x in iii[-5:]], f"{time.time() - t0:.0f}s", flush=True)
    out["gate_S3_iii_rows"] = iii
    out["gate_S3_iii_pass"] = bool(all(x["max_abs_diff_rad"] < 1e-6 for x in iii))
    out["gate_S3_pass"] = bool(out["gate_S3_i_pass_bonferroni"] and out["gate_S3_ii_pass"] and out["gate_S3_iii_pass"])
    out["runtime_s"] = time.time() - t0
    json.dump(out, open(RES / "stage3_gate.json", "w"), indent=1)
    print("GATE S3: i(bonferroni)", out["gate_S3_i_pass_bonferroni"], "i(uncorrected)", out["gate_S3_i_pass_uncorrected"],
          "ii", out["gate_S3_ii_pass"], "iii", out["gate_S3_iii_pass"], "failures", out["gate_S3_i_failures"], f"{out['runtime_s']:.0f}s")


if __name__ == "__main__":
    main()
