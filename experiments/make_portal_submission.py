"""Generate portal_submission.md for the ION Abstract Management Portal form: title, authors, plain-text abstract
(500-2500 words, no symbols), a separate 'what is new' field and up to three session choices.  Every number comes from
results/ and is listed in the claims check.  Text is restricted to ASCII (the form accepts text only)."""
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


def words(x):
    if x >= 1e9:
        return f"{x / 1e9:.1f} billion"
    if x >= 1e6:
        return f"{x / 1e6:.1f} million"
    return f"{round(x, -4):,.0f}"


def pc(x):
    return f"{100 * x:.0f} percent"


def main():
    claims = []

    def C(label, value, text, src, loc):
        claims.append((label, value, src, loc))
        return text

    s1 = jl("stage1_gate.json"); sp = jl("stage2_sprime.json")["gate_S2prime"]; g3 = jl("stage3_gate.json"); g3a = jl("stage3_gate_amended.json")
    du, dc = jl("stage3_decision_summary.json"), jl("stage3_cal_decision_summary.json")
    obs, tor = jl("stage3_observability.json"), jl("stage3_tornado_summary.json")
    margins = cl("stage2_margins.csv")
    F = cl("stage2_failure_rates.csv")
    orth = float(np.mean([float(r["fail_rate_gt_10deg"]) for r in F if r["set"] == "K3orth" and r["eps"] == "1.0"]))
    D = cl("stage3_cal_decision_map.csv")
    for r in D:
        r["ok"] = r["criterion_met"] == "True"
        for k in ("sigma_mag", "dip", "sigma_rp_deg", "reduction_600", "reduction_60", "sigma_d", "bias_dph", "arw"):
            r[k] = float(r[k])
    ps = [r for r in D if r["ok"]]
    n_two = sum(1 for r in ps if r["sigma_mag"] == 2.0 or r["dip"] == 15.0)
    noisy = float(np.mean([r["reduction_600"] for r in ps if r["sigma_mag"] == 2.0 and r["sigma_rp_deg"] == 0.3 and r["arm"] == "C"]))
    lowdip = [r["reduction_600"] for r in ps if r["dip"] == 15.0 and r["sigma_rp_deg"] == 0.3]
    base = [r for r in D if r["arm"] == "C" and r["sigma_mag"] == 0.5 and r["dip"] == 45 and r["sigma_d"] == 2.0 and r["bias_dph"] == 10 and r["arw"] == 0.3 and r["sigma_rp_deg"] == 0.3][0]
    imp = [r for r in cl("stage3_cal_B_improvement.csv") if float(r["bias_dph"]) == 10 and float(r["arw"]) == 0.3 and float(r["sigma_d"]) == 2 and float(r["dip"]) == 45 and float(r["sigma_mag"]) == 0.5][0]
    head = du["headline_mid_gyro_base_config"]
    A600 = head["sigma_d=0.0"]["A_600_deg"]; B0 = head["sigma_d=0.0"]["B_600_deg"]; B10 = head["sigma_d=10.0"]["B_600_deg"]
    S = cl("stage3_summary.csv")
    nb = [r for r in S if r["arm"] == "B" and float(r["bias_dph"]) == 10 and float(r["arw"]) == 0.3 and float(r["sigma_d"]) == 0 and float(r["dip"]) == 45 and float(r["sigma_mag"]) == 0.5][0]
    nees_lo = min(float(nb[k]) for k in ("nees_60", "nees_300", "nees_600")); nees_hi = max(float(nb[k]) for k in ("nees_60", "nees_300", "nees_600"))
    E = cl("stage3_exploratory_tuned.csv")
    e2 = [r for r in E if r["exp"] == "E2_no_soft_iron" and r["arm"] == "C" and float(r["sigma_d"]) == 0 and float(r["sigma_rp_deg"]) == 0.3][0]
    ob = {r["arm"]: r for r in obs if r["window"] == "hover_0_120s" and r["dip"] == 45.0}
    snap = cl("stage3_exploratory_snapshot.csv")
    flip = max(float(r["frac_yaw_gt_90deg_axis_only"]) for r in snap)
    req = req_table(); rrows = cl("stage2_requirements.csv")
    n_lo, _ = n_required(0.3, req["1"]); n_hi, _ = n_required(0.3, req["eps_relax_1e+07"])
    r1 = [x for x in rrows if x["set"] == "K3selected" and x["eps_name"] == "1" and x["realism"] == "baseline"][0]
    worst_ratio = float(r1["N_5deg_worst_dir"]) / float(r1["N_5deg_median_dir"])
    reach = [p["parameter"] for p in tor["cal"]["ranking"] if p["reaches_20pct"]]
    reach_u = [p["parameter"] for p in tor["uncal"]["ranking"] if p["reaches_20pct"]]
    k = cl("stage3_kappa_exploratory_corrected.csv")
    k_c10 = [r for r in k if r["arm"] == "C" and float(r["kappa"]) == 0.5 and float(r["sigma_d"]) == 10.0]
    marg = {r["model"]: float(r["aniso"]) for r in margins}
    r2 = marg["R2"]; relax7 = [v for kk, v in marg.items() if "relax_rate=10000000.0" in kk][0]

    v = {}
    v["an"] = C("minimal-model yield anisotropy", s1["primary"]["aniso"], f"{s1['primary']['aniso']:.2f}", "results/stage1_gate.json", "primary.aniso")
    v["r2"] = C("anisotropy of richer assumed model R2", r2, f"{r2:.2f}", "results/stage2_margins.csv", "model R2")
    v["rx"] = C("relaxation (10 million per second) anisotropy margin over the 0.001 threshold", relax7 / 1e-3, f"{relax7 / 1e-3:.1f}", "results/stage2_margins.csv", "relax_rate=1e7, aniso/1e-3")
    v["orth"] = C("axis-permuted element set ambiguity failure rate", orth, f"{orth:.2f}", "results/stage2_failure_rates.csv", "mean fail_rate_gt_10deg, K3orth, eps=1")
    v["rlo"] = C("Gate S2' ratio min", min(sp["large_N_ratios"]), f"{min(sp['large_N_ratios']):.2f}", "results/stage2_sprime.json", "gate_S2prime.large_N_ratios")
    v["rhi"] = C("Gate S2' ratio max", max(sp["large_N_ratios"]), f"{max(sp['large_N_ratios']):.2f}", "results/stage2_sprime.json", "gate_S2prime.large_N_ratios")
    v["slope"] = C("Gate S2' slope", sp["slope"], f"{sp['slope']:.2f}", "results/stage2_sprime.json", "gate_S2prime.slope")
    v["u"] = C("pre-registered passes", du["n_passing"], str(du["n_passing"]), "results/stage3_decision_summary.json", "n_passing")
    v["un"] = C("pre-registered configurations", du["n_configs"], str(du["n_configs"]), "results/stage3_decision_summary.json", "n_configs")
    v["A"] = C("gravity-only yaw RMSE at 600 s (deg)", A600, f"{A600:.1f}", "results/stage3_decision_summary.json", "headline A_600_deg")
    v["B0"] = C("magnetometer yaw RMSE at 600 s, no disturbance (deg)", B0, f"{B0:.2f}", "results/stage3_decision_summary.json", "headline sigma_d=0 B_600_deg")
    v["B10"] = C("magnetometer yaw RMSE at 600 s, 10 microtesla disturbance (deg)", B10, f"{B10:.2f}", "results/stage3_decision_summary.json", "headline sigma_d=10 B_600_deg")
    v["nlo"] = C("baseline NEES min over checkpoints", nees_lo, f"{nees_lo:.0f}", "results/stage3_summary.csv", "arm B mid gyro sigma_d 0 nees_60/300/600 min")
    v["nhi"] = C("baseline NEES max over checkpoints", nees_hi, f"{nees_hi:.0f}", "results/stage3_summary.csv", "arm B mid gyro sigma_d 0 nees_60/300/600 max")
    v["e60"] = C("no-soft-iron case, arm C 0.3 deg, 60 s", float(e2["reduction_60"]), pc(float(e2["reduction_60"])), "results/stage3_exploratory_tuned.csv", "E2, C, sigma_d 0, 0.3 deg")
    v["e600"] = C("no-soft-iron case, arm C 0.3 deg, 600 s", float(e2["reduction_600"]), pc(float(e2["reduction_600"])), "results/stage3_exploratory_tuned.csv", "E2, C, sigma_d 0, 0.3 deg")
    v["bu"] = C("B yaw RMSE 600 s uncalibrated (deg)", float(imp["rmse600_B_uncal_deg"]), f"{float(imp['rmse600_B_uncal_deg']):.2f}", "results/stage3_cal_B_improvement.csv", "mid gyro, sigma_d 2, 0.5 uT, 45 deg")
    v["bc"] = C("B yaw RMSE 600 s calibrated (deg)", float(imp["rmse600_B_cal_deg"]), f"{float(imp['rmse600_B_cal_deg']):.2f}", "results/stage3_cal_B_improvement.csv", "same row")
    v["c"] = C("calibrated passes", dc["n_passing"], str(dc["n_passing"]), "results/stage3_cal_decision_summary.json", "n_passing")
    v["cn"] = C("calibrated configurations", dc["n_configs"], str(dc["n_configs"]), "results/stage3_cal_decision_summary.json", "n_configs")
    v["two"] = C("passes with magnetometer noise 2 microtesla or dip 15 degrees", n_two, str(n_two), "results/stage3_cal_decision_map.csv", "criterion met and (sigma_mag == 2 or dip == 15)")
    v["noisy"] = C("arm C reduction at 600 s, 2 microtesla magnetometer, 0.3 deg", noisy, pc(noisy), "results/stage3_cal_decision_map.csv", "arm C, sigma_mag 2, 0.3 deg, criterion met")
    v["lo"] = C("reduction at 600 s, dip 15, 0.3 deg (min)", min(lowdip), pc(min(lowdip)), "results/stage3_cal_decision_map.csv", "dip 15, 0.3 deg, criterion met")
    v["hi"] = C("reduction at 600 s, dip 15, 0.3 deg (max)", max(lowdip), pc(max(lowdip)), "results/stage3_cal_decision_map.csv", "dip 15, 0.3 deg, criterion met")
    v["b60"] = C("base setting arm C 60 s", base["reduction_60"], f"{100 * base['reduction_60']:.0f} percent", "results/stage3_cal_decision_map.csv", "arm C, 0.5 uT, 45 deg, sigma_d 2, 0.3 deg")
    v["b600"] = C("base setting arm C 600 s", base["reduction_600"], f"{100 * base['reduction_600']:.0f} percent", "results/stage3_cal_decision_map.csv", "same row")
    v["nl"] = C("photons per element for 0.3 deg, ideal contrast", n_lo, words(n_lo), "results/stage2_requirements.csv + stage2_failure_threshold.json", "N_req(0.3, eps=1)")
    v["nh"] = C("photons per element for 0.3 deg, relaxation-limited contrast", n_hi, words(n_hi), "results/stage2_requirements.csv + stage2_failure_threshold.json", "N_req(0.3, eps_relax_1e7)")
    v["wr"] = C("worst-direction to median-direction photon ratio (5 deg, ideal contrast)", worst_ratio, f"{round(worst_ratio, -3):,.0f}", "results/stage2_requirements.csv", "K3selected eps=1 N_5deg_worst/median")
    v["flip"] = C("maximum snapshot 180-degree flip fraction, axis-only", flip, f"{100 * flip:.1f} percent", "results/stage3_exploratory_snapshot.csv", "max frac_yaw_gt_90deg_axis_only")
    v["ob"] = C("observability rank, gravity only", 5, "5", "results/stage3_observability.json", "arm A rank")
    v["reach"] = C("parameters lifting best-case calibrated benefit to 20 percent", "; ".join(reach), " and ".join(p.replace(" (deg)", "").replace(" (uT)", "").replace("MEMS noise sigma_mag", "magnetometer noise").replace("dip angle", "magnetic dip angle") for p in reach), "results/stage3_tornado_summary.json", "cal.ranking[*].reaches_20pct")
    v["kc"] = C("kappa 0.5, 10 microtesla, arm C criterion met", f"{sum(r['criterion_met'] == 'True' for r in k_c10)} of {len(k_c10)}", f"{sum(r['criterion_met'] == 'True' for r in k_c10)} of {len(k_c10)}", "results/stage3_kappa_exploratory_corrected.csv", "arm C, kappa 0.5, sigma_d 10")
    assert not reach_u, "text below assumes no pre-registered-filter parameter reaches 20 percent"

    title = "Can a Radical-Pair Chemical Compass Aid Inertial Attitude Estimation? A Simulation-Based Requirements Study"
    A = []
    A.append("OBJECTIVES. Magnetometer-based heading is fragile in GNSS-denied environments because motors, batteries and steel distort the local magnetic field. Radical-pair magnetoreception, the spin-chemical mechanism proposed for the avian compass, has inspired the idea of an auxiliary chemical compass whose output is an optical signal. "
             "This study asks under what conditions, if any, such a sensor could improve attitude and heading estimation in an inertial navigation filter relative to a MEMS magnetometer, and what sensor performance that would require. A negative answer was declared a valid outcome before any simulation was run. "
             "The work is simulation-only: no sensor was built, every sensor parameter is assumed and swept, and no new physics is claimed. Radical-pair precision limits, cryptochrome magnetic sensitivity and quantum magnetic navigation with highly sensitive magnetometers are established prior work. A literature search, which was limited, found no earlier evaluation of a radical-pair sensor as an aiding measurement inside an inertial filter.")
    A.append("METHODOLOGY. The study followed a written plan committed before any experiment. Later amendments were dated and committed before the runs they govern, and failed checks are reported rather than removed. "
             "A spin-dynamics model of a flavin-like radical pair, with two nitrogen hyperfine couplings whose tensors were checked against the primary literature, gives the direction dependence of the singlet yield. Two independent implementations of the yield agree to numerical precision. "
             "The yield feeds a photon-limited sensor made of several identically oriented elements. Counts are Poisson distributed, a maximum-likelihood estimator recovers the field axis, and a Cramer-Rao bound gives the corresponding accuracy limit. Because the yield is unchanged when the field is reversed, only the field axis, not its sign, is observable. "
             "The estimator output enters a six hundred second multiplicative quaternion attitude filter with gyroscope bias estimation, accelerometer gravity aiding, and one of four magnetic configurations: none, a MEMS magnetometer, the radical-pair axis alone with a sign-invariant residual, or both. "
             "The scripted trajectory has small tilt oscillations, three slow ninety-degree yaw turns and a figure-eight. Gyroscope bias instability, angle random walk, magnetometer noise, field dip angle, radical-pair accuracy and disturbance level are swept, with 200 paired Monte Carlo runs per configuration. "
             "The magnetic disturbance, made of soft iron, hard iron and a motor-current-correlated field, is applied identically to both magnetic sensors, so no immunity is credited to the radical-pair channel. "
             "The pre-registered decision criterion was a reduction of yaw error by at least twenty percent at both sixty and six hundred seconds relative to the magnetometer, with a paired bootstrap lower confidence bound above zero. "
             "A second study, specified after the first answer was known, gives both sensors the same online calibration of hard iron, soft iron and the motor-current field.")
    A.append(f"RESULTS. Spin model: the minimal model has a yield anisotropy of {v['an']}, an assumed richer hyperfine model gives {v['r2']}, and relaxation at ten million per second leaves {v['rx']} times a pre-registered threshold of 0.001, so the minimal model is an optimistic upper bound. "
             f"Sensor: element sets related by axis permutations have a four-fold mirror ambiguity caused by an exact molecular symmetry, with a failure rate (estimate more than ten degrees off) of {v['orth']}, and this check of the pre-registered plan failed. A non-orthogonal three-element set chosen by a pre-registered rule removes the ambiguity (ratio of simulated error to bound {v['rlo']} to {v['rhi']}, slope {v['slope']} on fresh data). "
             f"Navigation, pre-registered: {v['u']} of {v['un']} configurations met the criterion. Any magnetic aiding helps greatly (yaw error at six hundred seconds falls from {v['A']} degrees with gravity only to {v['B0']} degrees with a magnetometer, and {v['B10']} degrees with a ten microtesla disturbance), but the magnetometer filter was badly inconsistent because soft iron was not in its noise model (normalized estimation error squared between {v['nlo']} and {v['nhi']}, against about 6 for a consistent filter). "
             f"With no soft iron and no disturbance, a sensor twice as accurate as the magnetometer improved yaw by only {v['e60']} at sixty seconds and {v['e600']} at six hundred seconds. "
             f"Navigation, with online calibration: both filters became consistent and magnetometer yaw error at six hundred seconds fell from {v['bu']} to {v['bc']} degrees. {v['c']} of {v['cn']} configurations then met the criterion, {v['two']} of them where magnetometer noise was 2 microtesla (about {v['noisy']} improvement at six hundred seconds) or the field dip was 15 degrees ({v['lo']} to {v['hi']}), with the radical-pair sensor at 0.3 degrees per axis. At 0.5 microtesla and 45 degrees dip the change was {v['b60']} at sixty seconds and {v['b600']} at six hundred seconds. "
             f"Requirements: reaching 0.3 degrees per axis needs roughly {v['nl']} detected photons per element for an ideal contrast and {v['nh']} for a relaxation-limited contrast, and field directions where the response is nearly flat need about {v['wr']} times more. "
             f"Other findings: gravity-only aiding has numerical observability rank {v['ob']} of 6 (yaw unobservable) while all magnetic configurations are fully observable, so the axis-only channel loses nothing locally. A snapshot initialization from the axis alone flips yaw by 180 degrees in at most {v['flip']} of cases. In one-at-a-time sensitivity runs, only {v['reach']} lifted the best-case benefit to twenty percent, and no parameter did so for the pre-registered filters. "
             f"If the radical-pair channel responded to the disturbance only half as strongly as a MEMS magnetometer, the axis-only arm met the criterion in {v['kc']} cases at 10 microtesla; nothing here supports such an assumption.")
    A.append("CONCLUSIONS. When the magnetic disturbance affects both sensors equally, a radical-pair channel is not useful against a good, calibrated magnetometer. It could matter only against a poor magnetometer or at low dip angle, and only if the radical-pair sensor reaches sub-degree accuracy. "
             "Under these assumptions, calibrating the magnetometer is a stronger lever than adding a new magnetic sensing principle. Neither answer rests on measurements, so the results are requirements, not predictions of any real device.")
    A.append("SIGNIFICANCE AND LIMITATIONS. The study turns a qualitative bio-inspired concept into measurable targets, namely angular accuracy and contrast, and gives a reproducible, open framework for deciding whether such a sensor merits hardware work. "
             "The limits are stated plainly. The attitude filter failed its pre-registered consistency check on the original trajectory because of an unmodelled acceleration, and an amended check passed afterwards, which is a confirmation and not an independent test. "
             "The calibrated study followed the first answer and is reported beside it, not instead of it. Only attitude is estimated, not position. Sensor errors are modelled as isotropic Gaussian, the calibration uses idealized priors, 0.3 degrees per axis is the smallest accuracy simulated, and chemistry kinetics, temperature and optics are not modelled.")
    abstract = "\n\n".join(A)
    new = ("This work is not a new sensor or new physics. What is new is an open, reproducible framework that links a spin-chemistry yield model to a photon-limited estimator and then to an inertial attitude filter, so that the benefit of a radical-pair compass can be stated as requirements. "
           "Specific new elements: an axis-only, sign-invariant fusion model; the identification of a mirror ambiguity caused by molecular symmetry in multi-element arrays, with a selection rule that removes it; "
           "an evaluation in which both sensors suffer the same disturbance and receive the same calibration, so no immunity is assumed; and a pre-registered protocol with dated amendments in which failed checks and negative results are kept and reported.")
    wc = lambda s: len(re.findall(r"\S+", s))
    for s in (abstract, new, title):
        assert all(ord(c) < 128 for c in s), [c for c in s if ord(c) >= 128][:5]
        assert not [w for w in FORBIDDEN if w in s.lower()]
    assert 500 <= wc(abstract) <= 2500, wc(abstract)
    L = ["# portal_submission.md", "", "> Generated by `experiments/make_portal_submission.py` from `results/`. Text is plain ASCII (the portal accepts text only: no symbols, equations or images). Fields follow the ION Abstract Management Portal form.", "",
         "## Title of Paper/Presentation", "", title, "",
         "## List all authors (with organizational affiliations) in the order they will appear on the paper", "",
         "Ujesha Sivaramakrishnan [confirm spelling/preferred name], Vellore Institute of Technology, Vellore, India [confirm department]; [Nilesh: full name, affiliation]; [Bibhab: full name, affiliation]; [faculty co-author, if any: name, affiliation]", "",
         f"## Abstract text ({wc(abstract)} words; the form requires 500 to 2,500)", "", abstract, "",
         f"## Describe what is new and/or innovative about your presentation ({wc(new)} words)", "", new, "",
         "## Session choices (up to 3; confirm on the portal that none is invited-only)", "",
         "1. Multisensor Integrated Systems and Sensor Fusion Technologies (Track C)",
         "2. Resilient Inertial Navigation Systems and Alternative Sensors (Track A)",
         "3. Advances in MEMS-Based Inertial Sensors and Integrated Navigation Systems (Track A)", "",
         "## Consent", "", "The form asks you to consent to ION publishing your name, authors, affiliations and the abstract text on its website. Tick this only when every listed author agrees.", "",
         "## Claims check", "", "| Claim | Value used | Source file | Locator |", "|---|---|---|---|"]
    for label, value, src, loc in claims:
        vv = f"{value:.4g}" if isinstance(value, (float, int, np.floating)) else str(value)
        L.append(f"| {label} | {vv} | `{src}` | {loc} |")
    (ROOT / "portal_submission.md").write_text("\n".join(L) + "\n", encoding="utf8")
    print("abstract words", wc(abstract), "| new-field words", wc(new), "| characters", len(abstract))


if __name__ == "__main__":
    main()
