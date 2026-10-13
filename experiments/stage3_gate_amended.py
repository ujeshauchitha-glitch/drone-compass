"""Gate S3(i)' (Amendment 3, A3.1): condition (i) evaluated on trajectory T0 (linear acceleration = 0).
POST-HOC: the outcome was visible in the earlier diagnostic; this is a confirmation, not an independent test."""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
from scipy.stats import chi2

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from nav import attitude as A  # noqa: E402

RES = ROOT / "results"
M, STEPS, N_TESTS = 200, (600, 3000, 6000), 12


def interval(alpha):
    return chi2.ppf(alpha / 2, 6 * M) / M, chi2.ppf(1 - alpha / 2, 6 * M) / M


def main():
    sc0 = A.Scenario(a_scale=0.0)
    rnd = A.Randoms(M)
    lo95, hi95 = interval(0.05)
    loB, hiB = interval(0.05 / N_TESTS)
    nom = A.default_cfg(soft_iron=False, sigma_d=0.0, sigma_rp_deg=2.0)
    rows = []
    for arm in "ABCD":
        r = A.run_filter(arm, nom, sc0, rnd, nees_steps=STEPS)
        for k in STEPS:
            v = float(r["nees"][k].mean())
            rows.append(dict(arm=arm, t_s=k * A.DT, nees=v, in_uncorrected_95=bool(lo95 <= v <= hi95), in_bonferroni_95=bool(loB <= v <= hiB)))
        print(arm, [round(x["nees"], 2) for x in rows[-3:]], flush=True)
    out = dict(trajectory="T0 (linear acceleration = 0)", post_hoc=True, intervals=dict(uncorrected_95=[lo95, hi95], bonferroni_family_95=[loB, hiB]),
               rows=rows, gate_S3_i_prime_pass_bonferroni=bool(all(r["in_bonferroni_95"] for r in rows)),
               gate_S3_i_prime_pass_uncorrected=bool(all(r["in_uncorrected_95"] for r in rows)))
    g = json.load(open(RES / "stage3_gate.json"))
    out["S3_prime_overall_pass"] = bool(out["gate_S3_i_prime_pass_bonferroni"] and g["gate_S3_ii_pass"] and g["gate_S3_iii_pass"])
    json.dump(out, open(RES / "stage3_gate_amended.json", "w"), indent=1)
    print("S3(i)' Bonferroni:", out["gate_S3_i_prime_pass_bonferroni"], " uncorrected:", out["gate_S3_i_prime_pass_uncorrected"], " S3' overall:", out["S3_prime_overall_pass"])


if __name__ == "__main__":
    main()
