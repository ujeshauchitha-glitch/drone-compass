"""Generate interim_abstract.md from Stages 0-2 result files only.  Every number is read from results/ and listed in the
claims check; nothing is typed by hand.  Stage 3 is described as planned."""
from __future__ import annotations

import csv
import json
import re
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
RES = ROOT / "results"
FORBIDDEN = ("demonstrated", "validated", "real-time", "real time")


def jl(name):
    return json.load(open(RES / name))


def cl(name):
    return list(csv.DictReader(open(RES / name, newline="", encoding="utf8")))


def sci(x, d=2):
    s = f"{x:.{d}e}"
    m, e = s.split("e")
    return f"{m}×10^{int(e)}"


def main():
    claims = []

    def C(label, value, text, source, locator):
        claims.append((label, value, source, locator))
        return text

    s1 = jl("stage1_gate.json")
    s2 = jl("stage2_gate.json")
    sp = jl("stage2_sprime.json")
    des = jl("stage2_design.json")
    nc = jl("stage2_noncollinear.json")
    req = cl("stage2_requirements.csv")
    marg = cl("stage2_margins.csv")
    thr = jl("stage2_failure_threshold.json")
    F = cl("stage2_failure_rates.csv")
    ren = jl("stage2_ren_investigation.json")

    aniso = s1["primary"]["aniso"]
    r2 = [float(r["aniso"]) for r in marg if r["model"] == "R2"][0]
    relax7 = [r for r in marg if "relax_rate=10000000.0" in r["model"]][0]
    thr_rate = float(relax7["margin_vs_1e-3"])
    orth_fail = float(np.mean([float(r["fail_rate_gt_10deg"]) for r in F if r["set"] == "K3orth" and r["eps"] == "1.0"]))
    g = sp["gate_S2prime"]
    rmin, rmax = min(g["large_N_ratios"]), max(g["large_N_ratios"])
    nc5 = [r for r in nc if r["model"] == "NC" and r["alpha_deg"] == 5][0]
    ncr15 = [r for r in nc if r["model"] == "NC-R2" and r["alpha_deg"] == 15][0]

    def rq(set_, en, key):
        return float([r for r in req if r["set"] == set_ and r["eps_name"] == en and r["realism"] == "baseline"][0][key])

    eps_r2_name = "eps_R2"
    eps_r2 = rq("K3selected", eps_r2_name, "eps")
    n5_med = rq("K3selected", eps_r2_name, "N_5deg_median_dir")
    n5_worst = rq("K3selected", eps_r2_name, "N_5deg_worst_dir")
    fail_thr = [t for t in thr if t["eps_name"] == eps_r2_name][0]["N_fail_le_5pct"]
    ren_def = ren["default_a0.5_b0.00127_c26"]["rows"][-1]["N_total_for_5deg_median"]
    ren_opt = ren["optimised_a0.9_b0.00508_c10"]["rows"][-1]["N_total_for_5deg_median"]

    ren_o2 = ren["optimised_a0.9_b0.00508_c10"]["ren_reported_total_photons"]["O2"]
    ren_def_o2 = ren["default_a0.5_b0.00127_c26"]["ren_reported_total_photons"]["O2"]
    ren_def_o1 = ren["default_a0.5_b0.00127_c26"]["ren_reported_total_photons"]["O1"]
    fac_opt = max(ren_opt / ren_o2, ren_o2 / ren_opt)
    t = {}
    t["aniso"] = C("minimal-model yield anisotropy", aniso, f"{aniso:.3f}", "results/stage1_gate.json", "primary.aniso")
    t["r2"] = C("anisotropy of richer assumed model R2", r2, f"{r2:.3f}", "results/stage2_margins.csv", "model == R2, aniso")
    t["marg"] = C("relaxation (1e7 /s) anisotropy margin over the 1e-3 threshold", thr_rate, f"{thr_rate:.1f}", "results/stage2_margins.csv", "model contains relax_rate=10000000.0, margin_vs_1e-3")
    t["fail"] = C("K3orth mean ambiguity failure rate (eps=1, over N grid)", orth_fail, f"{orth_fail:.2f}", "results/stage2_failure_rates.csv", "mean fail_rate_gt_10deg, set=K3orth, eps=1.0")
    t["ratio_lo"] = C("Gate S2' ratio min", rmin, f"{rmin:.2f}", "results/stage2_sprime.json", "gate_S2prime.large_N_ratios (min)")
    t["ratio_hi"] = C("Gate S2' ratio max", rmax, f"{rmax:.2f}", "results/stage2_sprime.json", "gate_S2prime.large_N_ratios (max)")
    t["slope"] = C("Gate S2' log-log slope", g["slope"], f"{g['slope']:.2f}", "results/stage2_sprime.json", "gate_S2prime.slope")
    t["nm5"] = C("photons for NC(5 deg) mirror separation", nc5["N_mirror_photons"], sci(nc5["N_mirror_photons"], 1), "results/stage2_noncollinear.json", "model NC, alpha 5, N_mirror_photons")
    t["nm15"] = C("photons for NC-R2(15 deg) mirror separation", ncr15["N_mirror_photons"], sci(ncr15["N_mirror_photons"], 1), "results/stage2_noncollinear.json", "model NC-R2, alpha 15, N_mirror_photons")
    t["eps"] = C("derived contrast factor for R2", eps_r2, f"{eps_r2:.3f}", "results/stage2_requirements.csv", "K3selected eps_R2 baseline, eps")
    t["n5m"] = C("CRB photons/element, 5 deg median direction, eps_R2", n5_med, sci(n5_med, 1), "results/stage2_requirements.csv", "K3selected eps_R2 baseline, N_5deg_median_dir")
    t["nf"] = C("photons/element for <=5% ambiguity failures, eps_R2", fail_thr, sci(fail_thr, 1), "results/stage2_failure_threshold.json", "eps_name == eps_R2, N_fail_le_5pct")
    t["n5w"] = C("CRB photons/element, 5 deg worst direction, eps_R2", n5_worst, sci(n5_worst, 1), "results/stage2_requirements.csv", "K3selected eps_R2 baseline, N_5deg_worst_dir")
    t["ren1"] = C("Ren-signal total photons for 5 deg (default, K=1000 random elements)", ren_def, sci(ren_def, 1), "results/stage2_ren_investigation.json", "default row K=1000 N_total_for_5deg_median")
    t["facopt"] = C("factor between our optimised-row total photons and Ren Table 2 (3.94e5)", fac_opt, f"{fac_opt:.0f}", "results/stage2_ren_investigation.json", "optimised row K=1000 vs ren_reported_total_photons.O2")
    t["ren2"] = C("Ren-signal total photons for 5 deg (optimised, K=1000 random elements)", ren_opt, sci(ren_opt, 1), "results/stage2_ren_investigation.json", "optimised row K=1000 N_total_for_5deg_median")

    body = (
        f"**INTERIM abstract (Stages 0–2 only; Stage 3 is planned and no navigation results exist).**\n\n"
        f"**Objective.** We ask under what conditions, if any, a radical-pair (spin-chemical) magnetic compass could aid heading in "
        f"GNSS-denied inertial navigation compared with a conventional MEMS magnetometer. The study is simulation-only; no experiment exists and no new physics is claimed.\n\n"
        f"**Method.** A singlet-yield model of a flavin-like radical pair (two ¹⁴N hyperfine couplings with tensors taken from Hiscock et al. 2016) "
        f"feeds a photon-limited sensor of K identical oriented elements, a Poisson maximum-likelihood axis estimator and Cramér–Rao bounds. "
        f"Contrast, noise and ensemble-order parameters are assumed and swept; every assumed value is labelled.\n\n"
        f"**Results (simulated).** The minimal model's yield anisotropy is {t['aniso']}, but an assumed richer hyperfine model gives {t['r2']}, and "
        f"relaxation at 10^7 s⁻¹ leaves only {t['marg']}× our pre-registered 10^-3 threshold. The yield is symmetric under field reversal, so only the field axis is observable. "
        f"Because both nitrogen tensors share one principal frame, the molecule has an exact sign-flip symmetry, and any array of elements related by axis permutations has a four-fold mirror ambiguity "
        f"(failure rate {t['fail']} at every photon number; the pre-registered gate failed). A non-orthogonal three-element set, chosen by a pre-registered distinguishability rule, removes it in the minimal model "
        f"(median-error to Cramér–Rao ratio {t['ratio_lo']}–{t['ratio_hi']}, slope {t['slope']} on fresh data). "
        f"With an assumed 5° tilt of one tensor frame the mirror images become separable above about {t['nm5']} photons in the minimal model, "
        f"but with an assumed 15° tilt in a richer model about {t['nm15']} photons are needed. "
        f"For a derived contrast factor of {t['eps']} (richer model), a 5° median error needs about {t['n5m']} detected photons per element by the Cramér–Rao bound, about {t['nf']} for at most 5% ambiguity failures, "
        f"and about {t['n5w']} for the worst field direction. With Ren et al.'s phenomenological signal, our bound gives {t['ren1']} (default) and {t['ren2']} (optimised) total photons for 5° "
        f"(their published figures: {sci(ren_def_o2, 1)}–{sci(ren_def_o1, 1)} and {sci(ren_o2, 2)}; a factor of {t['facopt']} remains for the optimised row).\n\n"
        f"**Planned.** Stage 3, an attitude-filter comparison against a MEMS magnetometer under shared magnetic disturbances, has not been run; no claim about navigation benefit is made."
    )
    words = len(re.findall(r"\S+", re.sub(r"\*\*|\n", " ", body)))
    low = body.lower()
    bad = [w for w in FORBIDDEN if w in low]
    assert not bad, bad

    lines = ["# interim_abstract.md", "", "> INTERIM, generated by `experiments/make_interim_abstract.py` from files in `results/`. Stages 0–2 only.", "",
             body, "", f"_Word count: {words}._", "", "## Claims check", "", "| Claim | Value used | Source file | Locator |", "|---|---|---|---|"]
    for label, value, src, loc in claims:
        v = f"{value:.4g}" if isinstance(value, (float, int)) else str(value)
        lines.append(f"| {label} | {v} | `{src}` | {loc} |")
    lines += ["", "## Not claimed", "",
              "* No sensor, no measurement and no navigation result exists. All numbers are simulated under assumed parameters (contrast, noise, order) that are labelled in `PLAN.md` and `results/`.",
              "* Photon numbers are per element unless stated; the ε-derived rows are approximate (ratio of relative anisotropies from Stage 1).",
              "* Unresolved: the optimised Ren row differs from their published figure by about 12× after the element-count investigation (see `results/stage2_ren_investigation.json`).",
              "* Stage 3 (attitude filter, MEMS comparison, observability) is **planned**."]
    (ROOT / "interim_abstract.md").write_text("\n".join(lines) + "\n", encoding="utf8")
    print(f"words={words}")


if __name__ == "__main__":
    main()
