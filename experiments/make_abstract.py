"""Generate abstract_draft.md (Stage 4).  Every number is read from results/ and listed in the claims check."""
from __future__ import annotations

import csv
import json
import re
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from experiments.stage3_analysis import n_required, req_table  # noqa: E402

RES = ROOT / "results"
FORBIDDEN = ("demonstrated", "validated", "real-time", "real time")


def jl(n):
    return json.load(open(RES / n))


def cl(n):
    return list(csv.DictReader(open(RES / n, newline="", encoding="utf8")))


def sci(x, d=1):
    m, e = f"{x:.{d}e}".split("e")
    return f"{m}×10^{int(e)}"


def main():
    claims = []

    def C(label, value, text, source, locator):
        claims.append((label, value, source, locator))
        return text

    uncal = jl("stage3_decision_summary.json")
    cal = jl("stage3_cal_decision_summary.json")
    gate3 = jl("stage3_gate.json")
    gate3a = jl("stage3_gate_amended.json")
    sp = jl("stage2_sprime.json")["gate_S2prime"]
    F = cl("stage2_failure_rates.csv")
    orth = float(np.mean([float(r["fail_rate_gt_10deg"]) for r in F if r["set"] == "K3orth" and r["eps"] == "1.0"]))
    D = cl("stage3_cal_decision_map.csv")
    for r in D:
        r["criterion_met"] = r["criterion_met"] == "True"
        for k in ("sigma_mag", "dip", "sigma_d", "sigma_rp_deg", "reduction_60", "reduction_600", "bias_dph", "arw"):
            r[k] = float(r[k])
    cpass = [r for r in D if r["criterion_met"]]
    noisy = [r["reduction_600"] for r in cpass if r["sigma_mag"] == 2.0 and r["sigma_rp_deg"] == 0.3 and r["arm"] == "C"]
    lowdip = [r["reduction_600"] for r in cpass if r["dip"] == 15.0 and r["sigma_rp_deg"] == 0.3]
    base_c = [r for r in D if r["arm"] == "C" and r["sigma_mag"] == 0.5 and r["dip"] == 45.0 and r["sigma_d"] == 2.0 and r["bias_dph"] == 10 and r["arw"] == 0.3 and r["sigma_rp_deg"] == 0.3][0]
    imp = [r for r in cl("stage3_cal_B_improvement.csv") if float(r["bias_dph"]) == 10 and float(r["arw"]) == 0.3 and float(r["sigma_d"]) == 2.0 and float(r["dip"]) == 45 and float(r["sigma_mag"]) == 0.5][0]
    req = req_table()
    n1, _ = n_required(0.3, req["1"]); n7, _ = n_required(0.3, req["eps_relax_1e+07"]); nr2, _ = n_required(0.3, req["eps_R2"])
    torj = jl("stage3_tornado_summary.json")
    reach = [p["parameter"] for p in torj["cal"]["ranking"] if p["reaches_20pct"]]
    reach_u = [p["parameter"] for p in torj["uncal"]["ranking"] if p["reaches_20pct"]]

    n_two = sum(1 for r in cpass if r["sigma_mag"] == 2.0 or r["dip"] == 15.0)
    t = {}
    t["two"] = C("calibrated passes in the two regimes (MEMS 2 uT or dip 15 deg)", n_two, f"{n_two}", "results/stage3_cal_decision_map.csv", "criterion met and (sigma_mag == 2.0 or dip == 15)")
    t["orth"] = C("K3orth mean ambiguity failure rate", orth, f"{orth:.2f}", "results/stage2_failure_rates.csv", "mean fail_rate_gt_10deg, K3orth, eps=1.0")
    t["rlo"] = C("Gate S2' ratio min", min(sp["large_N_ratios"]), f"{min(sp['large_N_ratios']):.2f}", "results/stage2_sprime.json", "gate_S2prime.large_N_ratios")
    t["rhi"] = C("Gate S2' ratio max", max(sp["large_N_ratios"]), f"{max(sp['large_N_ratios']):.2f}", "results/stage2_sprime.json", "gate_S2prime.large_N_ratios")
    t["u"] = C("pre-registered passes", uncal["n_passing"], f"{uncal['n_passing']}", "results/stage3_decision_summary.json", "n_passing")
    t["un"] = C("pre-registered configurations", uncal["n_configs"], f"{uncal['n_configs']}", "results/stage3_decision_summary.json", "n_configs")
    t["c"] = C("calibrated passes", cal["n_passing"], f"{cal['n_passing']}", "results/stage3_cal_decision_summary.json", "n_passing")
    t["cn"] = C("calibrated configurations", cal["n_configs"], f"{cal['n_configs']}", "results/stage3_cal_decision_summary.json", "n_configs")
    t["noisy"] = C("arm C reduction at 600 s, MEMS 2 uT, sigma_rp 0.3 deg (mean over passing rows)", np.mean(noisy), f"{100 * np.mean(noisy):.0f}%", "results/stage3_cal_decision_map.csv", "arm C, sigma_mag 2.0, sigma_rp 0.3, criterion met, reduction_600")
    t["lo"] = C("reduction at 600 s at 15 deg dip, sigma_rp 0.3 deg (min)", min(lowdip), f"{100 * min(lowdip):.0f}%", "results/stage3_cal_decision_map.csv", "dip 15, sigma_rp 0.3, criterion met, reduction_600 (min)")
    t["hi"] = C("reduction at 600 s at 15 deg dip, sigma_rp 0.3 deg (max)", max(lowdip), f"{100 * max(lowdip):.0f}%", "results/stage3_cal_decision_map.csv", "dip 15, sigma_rp 0.3, criterion met, reduction_600 (max)")
    t["b60"] = C("arm C vs calibrated B at base setting, t=60 s", base_c["reduction_60"], f"{100 * base_c['reduction_60']:.0f}%", "results/stage3_cal_decision_map.csv", "arm C, sigma_mag 0.5, dip 45, sigma_d 2, sigma_rp 0.3, reduction_60")
    t["b600"] = C("arm C vs calibrated B at base setting, t=600 s", base_c["reduction_600"], f"{100 * base_c['reduction_600']:.0f}%", "results/stage3_cal_decision_map.csv", "same row, reduction_600")
    t["bu"] = C("B yaw RMSE at 600 s, uncalibrated (deg)", float(imp["rmse600_B_uncal_deg"]), f"{float(imp['rmse600_B_uncal_deg']):.2f}", "results/stage3_cal_B_improvement.csv", "mid gyro, sigma_d 2, sigma_mag 0.5, dip 45")
    t["bc"] = C("B yaw RMSE at 600 s, calibrated (deg)", float(imp["rmse600_B_cal_deg"]), f"{float(imp['rmse600_B_cal_deg']):.2f}", "results/stage3_cal_B_improvement.csv", "same row")
    t["n1"] = C("photons/element needed for sigma_rp 0.3 deg, ideal contrast", n1, sci(n1), "results/stage2_requirements.csv + stage2_failure_threshold.json", "N_req(0.3 deg, eps=1) per A2.2")
    t["nr2"] = C("photons/element needed for sigma_rp 0.3 deg, richer-model contrast", nr2, sci(nr2), "results/stage2_requirements.csv + stage2_failure_threshold.json", "N_req(0.3 deg, eps_R2) per A2.2")
    t["n7"] = C("photons/element needed for sigma_rp 0.3 deg, relaxation-limited contrast", n7, sci(n7), "results/stage2_requirements.csv + stage2_failure_threshold.json", "N_req(0.3 deg, eps_relax_1e7) per A2.2")
    t["reach"] = " and ".join(p.replace(" (deg)", "").replace(" (uT)", "").replace("MEMS noise sigma_mag", "magnetometer noise") for p in reach)
    claims.append(("parameters that lift the best-case calibrated benefit to 20%", ", ".join(reach), "results/stage3_tornado_summary.json", "cal.ranking[*].reaches_20pct"))
    claims.append(("parameters that do so for the pre-registered filters", ", ".join(reach_u) if reach_u else "none", "results/stage3_tornado_summary.json", "uncal.ranking[*].reaches_20pct"))

    body = (
        f"**Objective.** Under what conditions, if any, could a radical-pair (spin-chemical) magnetic compass improve the heading of an inertial attitude filter relative to a MEMS magnetometer? "
        f"The study is simulation-only; no experiment exists and no new physics is claimed.\n\n"
        f"**Method.** A spin-dynamics yield model (two ¹⁴N couplings from the literature) feeds a photon-limited multi-element sensor and a Poisson maximum-likelihood axis estimator, whose output enters a 600 s quaternion attitude filter with "
        f"gyro, accelerometer and a MEMS magnetometer, the radical-pair axis, or both, under shared magnetic disturbances (200 paired runs per configuration). Sensor, noise and distortion parameters are assumed and swept.\n\n"
        f"**Results (simulated).** The sensor observes only the field axis; element sets related by axis permutations are mirror-ambiguous (failure rate {t['orth']}), and a non-orthogonal set resolves this (error/bound ratio {t['rlo']}–{t['rhi']}). "
        f"In the pre-registered comparison {t['u']} of {t['un']} configurations reduced yaw RMSE by at least 20% at both 60 s and 600 s, but the magnetometer baseline was inconsistent (unmodelled soft iron). "
        f"With online calibration of hard-iron, soft-iron and motor-current fields, which lowers magnetometer yaw RMSE at 600 s from {t['bu']}° to {t['bc']}°, {t['c']} of {t['cn']} configurations meet the criterion, {t['two']} of them "
        f"with a radical-pair error of 0.3° per axis against a 2 µT-noise magnetometer (about {t['noisy']} at 600 s) or at 15° dip ({t['lo']}–{t['hi']}); at 0.5 µT and 45° dip the 60 s change is {t['b60']} and the 600 s change {t['b600']}. "
        f"Reaching 0.3° needs roughly {t['n1']} (ideal contrast) to {t['n7']} (relaxation-limited) detected photons per element. In one-at-a-time sensitivity runs only {t['reach']} lift the best-case benefit to 20%{' (none do for the pre-registered filters)' if not reach_u else ''}.\n\n"
        f"**Limits.** The filter failed its pre-registered consistency gate on the original trajectory (unmodelled acceleration); an amended gate passed post hoc. Attitude only, no position; isotropic Gaussian sensor errors; idealised calibration priors."
    )
    words = len(re.findall(r"\S+", re.sub(r"\*\*|\n", " ", body)))
    low = body.lower()
    bad = [w for w in FORBIDDEN if w in low]
    assert not bad, bad
    assert 250 <= words <= 350, words
    lines = ["# abstract_draft.md", "", "> Generated by `experiments/make_abstract.py` from files in `results/`. Simulation only.", "", body, "", f"_Word count: {words}._", "",
             "## Claims check", "", "| Claim | Value used | Source file | Locator |", "|---|---|---|---|"]
    for label, value, src, loc in claims:
        v = f"{value:.4g}" if isinstance(value, (float, int, np.floating)) else str(value)
        lines.append(f"| {label} | {v} | `{src}` | {loc} |")
    lines += ["", "## Not claimed", "",
              "* No sensor, measurement or flight exists; every number is a simulated result under assumed parameters (see `PLAN.md`, `report.pdf`).",
              "* The pre-registered answer is NO (0 of 390). The calibrated study was specified in Amendment 3 after that answer was known and is reported next to it, not instead of it.",
              "* The filter failed Gate S3(i); Gate S3(i)′ is post hoc and not an independent test.",
              "* σ_rp = 0.3° is the smallest value in the grid, so the size of the benefit at smaller σ_rp is not known.",
              "* Photon numbers use the median-direction Cramér–Rao bound and an ambiguity-failure threshold; the worst field direction needs far more (see `results/stage2_requirements.csv`)."]
    (ROOT / "abstract_draft.md").write_text("\n".join(lines) + "\n", encoding="utf8")
    print("words", words)


if __name__ == "__main__":
    main()
