"""Generate abstract_submission.md: the five elements the PLANS 2027 call asks for (objectives, results, conclusions,
key innovative steps, significance), with every number read from results/ and listed in a claims check."""
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


def sci(x):
    m, e = f"{x:.1e}".split("e")
    return f"{m}×10^{int(e)}"


def main():
    claims = []

    def C(label, value, text, source, locator):
        claims.append((label, value, source, locator))
        return text

    du, dc = jl("stage3_decision_summary.json"), jl("stage3_cal_decision_summary.json")
    D = cl("stage3_cal_decision_map.csv")
    for r in D:
        r["ok"] = r["criterion_met"] == "True"
        for k in ("sigma_mag", "dip", "sigma_rp_deg", "reduction_600", "reduction_60", "sigma_d", "bias_dph", "arw"):
            r[k] = float(r[k])
    ps = [r for r in D if r["ok"]]
    n_two = sum(1 for r in ps if r["sigma_mag"] == 2.0 or r["dip"] == 15.0)
    noisy = np.mean([r["reduction_600"] for r in ps if r["sigma_mag"] == 2.0 and r["sigma_rp_deg"] == 0.3 and r["arm"] == "C"])
    lowdip = [r["reduction_600"] for r in ps if r["dip"] == 15.0 and r["sigma_rp_deg"] == 0.3]
    base = [r for r in D if r["arm"] == "C" and r["sigma_mag"] == 0.5 and r["dip"] == 45 and r["sigma_d"] == 2.0 and r["bias_dph"] == 10 and r["arw"] == 0.3 and r["sigma_rp_deg"] == 0.3][0]
    imp = [r for r in cl("stage3_cal_B_improvement.csv") if float(r["bias_dph"]) == 10 and float(r["arw"]) == 0.3 and float(r["sigma_d"]) == 2 and float(r["dip"]) == 45 and float(r["sigma_mag"]) == 0.5][0]
    req = req_table()
    n_lo, _ = n_required(0.3, req["1"]); n_hi, _ = n_required(0.3, req["eps_relax_1e+07"])
    F = cl("stage2_failure_rates.csv")
    orth = float(np.mean([float(r["fail_rate_gt_10deg"]) for r in F if r["set"] == "K3orth" and r["eps"] == "1.0"]))
    tor = jl("stage3_tornado_summary.json")

    v = {}
    v["u"] = C("pre-registered passes", du["n_passing"], str(du["n_passing"]), "results/stage3_decision_summary.json", "n_passing")
    v["un"] = C("pre-registered configurations", du["n_configs"], str(du["n_configs"]), "results/stage3_decision_summary.json", "n_configs")
    v["c"] = C("calibrated passes", dc["n_passing"], str(dc["n_passing"]), "results/stage3_cal_decision_summary.json", "n_passing")
    v["cn"] = C("calibrated configurations", dc["n_configs"], str(dc["n_configs"]), "results/stage3_cal_decision_summary.json", "n_configs")
    v["two"] = C("passes with MEMS 2 uT or dip 15 deg", n_two, str(n_two), "results/stage3_cal_decision_map.csv", "criterion met and (sigma_mag == 2 or dip == 15)")
    v["noisy"] = C("arm C reduction at 600 s, MEMS 2 uT, sigma_rp 0.3", noisy, f"{100 * noisy:.0f}%", "results/stage3_cal_decision_map.csv", "arm C, sigma_mag 2.0, sigma_rp 0.3, criterion met")
    v["lo"] = C("reduction at 600 s, dip 15, sigma_rp 0.3 (min)", min(lowdip), f"{100 * min(lowdip):.0f}%", "results/stage3_cal_decision_map.csv", "dip 15, sigma_rp 0.3, criterion met")
    v["hi"] = C("reduction at 600 s, dip 15, sigma_rp 0.3 (max)", max(lowdip), f"{100 * max(lowdip):.0f}%", "results/stage3_cal_decision_map.csv", "dip 15, sigma_rp 0.3, criterion met")
    v["b60"] = C("base setting, arm C, 60 s", base["reduction_60"], f"{100 * base['reduction_60']:.0f}%", "results/stage3_cal_decision_map.csv", "arm C, 0.5 uT, 45 deg, sigma_d 2, sigma_rp 0.3")
    v["b600"] = C("base setting, arm C, 600 s", base["reduction_600"], f"{100 * base['reduction_600']:.0f}%", "results/stage3_cal_decision_map.csv", "same row")
    v["bu"] = C("B yaw RMSE 600 s uncalibrated (deg)", float(imp["rmse600_B_uncal_deg"]), f"{float(imp['rmse600_B_uncal_deg']):.2f}", "results/stage3_cal_B_improvement.csv", "mid gyro, sigma_d 2, 0.5 uT, 45 deg")
    v["bc"] = C("B yaw RMSE 600 s calibrated (deg)", float(imp["rmse600_B_cal_deg"]), f"{float(imp['rmse600_B_cal_deg']):.2f}", "results/stage3_cal_B_improvement.csv", "same row")
    v["nlo"] = C("photons/element for 0.3 deg, ideal contrast", n_lo, sci(n_lo), "results/stage2_requirements.csv + stage2_failure_threshold.json", "N_req(0.3, eps=1)")
    v["nhi"] = C("photons/element for 0.3 deg, relaxation-limited contrast", n_hi, sci(n_hi), "results/stage2_requirements.csv + stage2_failure_threshold.json", "N_req(0.3, eps_relax_1e7)")
    v["orth"] = C("axis-permuted element set ambiguity failure rate", orth, f"{orth:.2f}", "results/stage2_failure_rates.csv", "mean fail_rate_gt_10deg, K3orth, eps=1")
    reach = [p["parameter"] for p in tor["cal"]["ranking"] if p["reaches_20pct"]]
    claims.append(("parameters that lift best-case calibrated benefit to 20%", "; ".join(reach), "results/stage3_tornado_summary.json", "cal.ranking[*].reaches_20pct"))

    title = "Can a Radical-Pair Chemical Compass Aid Inertial Attitude Estimation? A Simulation-Based Requirements Study"
    full = (
        f"**Objectives.** To determine under what conditions a radical-pair (spin-chemical) magnetic compass could improve attitude and heading estimation for GNSS-denied inertial navigation relative to a MEMS magnetometer, "
        f"and what sensor performance that would require. The work is simulation-only: no sensor was built, sensor parameters are assumed and swept, and no new physics is claimed.\n\n"
        f"**Results (simulated).** A spin-dynamics yield model feeds a photon-limited multi-element sensor, a maximum-likelihood axis estimator and a quaternion attitude filter (200 paired runs per configuration). "
        f"Against an uncalibrated magnetometer, none of {v['un']} pre-registered configurations gave a 20% yaw-error reduction at both 60 s and 600 s; that baseline proved inconsistent. "
        f"With online calibration of hard-iron, soft-iron and motor-current fields (magnetometer yaw RMSE at 600 s falls from {v['bu']}° to {v['bc']}°), {v['c']} of {v['cn']} configurations meet the criterion, {v['two']} of them where the magnetometer noise is 2 µT "
        f"(about {v['noisy']} at 600 s) or the field dip is 15° ({v['lo']}–{v['hi']}) and the radical-pair sensor reaches 0.3° per axis. At 0.5 µT and 45° dip the changes are {v['b60']} (60 s) and {v['b600']} (600 s).\n\n"
        f"**Conclusions.** With the magnetic disturbance shared by both sensors, a radical-pair channel is not useful against a good, calibrated magnetometer. It could matter against a poor magnetometer or at low dip, "
        f"and only at sub-degree accuracy, which needs roughly {v['nlo']} to {v['nhi']} detected photons per element depending on the assumed contrast.\n\n"
        f"**Key innovative steps.** An axis-only, sign-invariant measurement model; the identification of a mirror ambiguity caused by molecular symmetry (failure rate {v['orth']} for axis-permuted elements) and its removal by a selected non-orthogonal element set; "
        f"identical disturbance and calibration for both sensors, with no immunity credited; and a pre-registered protocol whose failed checks are reported.\n\n"
        f"**Significance.** The study turns a qualitative bio-inspired idea into measurable targets, namely angular accuracy and contrast, before any hardware is built, and indicates that magnetometer calibration is the stronger lever under these assumptions. "
        f"Attitude only; position is not evaluated."
    )
    short = (
        f"We ask when a radical-pair (spin-chemical) magnetic compass could improve attitude estimation for GNSS-denied inertial navigation relative to a MEMS magnetometer. In simulation only (no sensor built; parameters assumed), "
        f"a spin-yield model, photon-limited estimator and quaternion filter were compared over {v['un']} pre-registered and {v['cn']} calibrated configurations. None of the pre-registered configurations gave a 20% yaw-error reduction at both 60 s and 600 s "
        f"(the baseline was inconsistent). With online magnetometer calibration, {v['c']} of {v['cn']} did, mainly with a noisy 2 µT magnetometer or 15° dip and a radical-pair sensor at 0.3° per axis. "
        f"Under shared disturbance the channel is not useful against a good calibrated magnetometer. Key steps: axis-only sign-invariant fusion, a symmetry-induced mirror ambiguity and its removal, and a pre-registered protocol. "
        f"Significance: measurable sensor targets before hardware."
    )
    wc = lambda s: len(re.findall(r"\S+", re.sub(r"\*\*|\n", " ", s)))
    for s in (full, short):
        assert not [w for w in FORBIDDEN if w in s.lower()]
    lines = ["# abstract_submission.md", "", "> Generated by `experiments/make_submission_abstract.py` from `results/`. Structured around the five elements the IEEE/ION PLANS 2027 call asks for.", "",
             "## Title", "", title, "",
             "## Where to submit", "",
             "* Portal: ion.org/abstracts (create an account, choose PLANS), **due 16 October 2026**; acceptance late November; full manuscript due 1 February 2027 (manuscript must be representative of this abstract).",
             "* Best-fit sessions (check on the portal that the one you pick is not invited-only): **Track C – Multisensor Integrated Systems and Sensor Fusion Technologies** (magnetometer, Bayesian filtering, sensor calibration, performance evaluation) or **Track A – Resilient Inertial Navigation Systems and Alternative Sensors** (GNSS-denied, alternative sensors). Weaker fits: Track A – Quantum Navigation Technologies (that session expects quantum-sensor hardware and field tests) and Track A – Advances in MEMS-Based Inertial Sensors and Integrated Navigation Systems.",
             "* Suggested keywords: radical-pair magnetoreception; chemical compass; magnetometer calibration; attitude estimation; GNSS-denied navigation; Kalman filtering; quantum biology.",
             "* Authors/affiliations: to be filled in (the presenting author must attend in person, pay registration, and attend the speakers' breakfast). The PLANS Student Award and student registration may apply to full-time students.", "",
             f"## Abstract (full, {wc(full)} words)", "", full, "",
             f"## Abstract (short fallback if the portal has a tight limit, {wc(short)} words)", "", short, "",
             "## Claims check", "", "| Claim | Value used | Source file | Locator |", "|---|---|---|---|"]
    for label, value, src, loc in claims:
        vv = f"{value:.4g}" if isinstance(value, (float, int, np.floating)) else str(value)
        lines.append(f"| {label} | {vv} | `{src}` | {loc} |")
    lines += ["", "## Not claimed", "",
              "* No measurement or hardware; all sensor parameters (contrast, error model, noise) are assumed.",
              "* The pre-registered answer is NO; the calibrated study was specified after it and is reported beside it. The filter failed its pre-registered consistency gate; the amended gate passed post hoc.",
              "* 0.3° per axis is the smallest error value simulated; benefit at smaller values is unknown. Attitude only, no position."]
    (ROOT / "abstract_submission.md").write_text("\n".join(lines) + "\n", encoding="utf8")
    print("full words", wc(full), "| short words", wc(short), "| characters full/short", len(re.sub(r"\*\*", "", full)), len(short))


if __name__ == "__main__":
    main()
