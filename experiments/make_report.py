"""Build report.pdf (Stage 4) from files in results/.  Tables and numbers are generated, not typed."""
from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

import numpy as np
from reportlab.lib import colors
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
from reportlab.lib.units import cm
from reportlab.platypus import Image, KeepTogether, PageBreak, Paragraph, SimpleDocTemplate, Spacer, Table, TableStyle

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from experiments.stage3_analysis import n_required, req_table  # noqa: E402

RES = ROOT / "results"
ss = getSampleStyleSheet()
BODY = ParagraphStyle("b", parent=ss["BodyText"], fontSize=8.6, leading=11.0, spaceAfter=3.5)
SMALL = ParagraphStyle("s", parent=BODY, fontSize=7.4, leading=9.2)
H1 = ParagraphStyle("h1", parent=ss["Heading1"], fontSize=13, leading=15, spaceBefore=6, spaceAfter=4)
H2 = ParagraphStyle("h2", parent=ss["Heading2"], fontSize=10, leading=12, spaceBefore=5, spaceAfter=2)
TITLE = ParagraphStyle("t", parent=ss["Title"], fontSize=16, leading=19, spaceAfter=3)
CELL = ParagraphStyle("c", parent=BODY, fontSize=7.0, leading=8.6, spaceAfter=0)
CELLB = ParagraphStyle("cb", parent=CELL, fontName="Helvetica-Bold")


def jl(n):
    return json.load(open(RES / n))


def cl(n):
    return list(csv.DictReader(open(RES / n, newline="", encoding="utf8")))


def P(t, st=BODY):
    return Paragraph(t, st)


def table(data, widths, header=True, font=CELL):
    cells = [[Paragraph(str(c), CELLB if (header and i == 0) else font) for c in row] for i, row in enumerate(data)]
    t = Table(cells, colWidths=widths, repeatRows=1 if header else 0)
    st = [("GRID", (0, 0), (-1, -1), 0.25, colors.HexColor("#999999")), ("VALIGN", (0, 0), (-1, -1), "TOP"),
          ("LEFTPADDING", (0, 0), (-1, -1), 2.5), ("RIGHTPADDING", (0, 0), (-1, -1), 2.5), ("TOPPADDING", (0, 0), (-1, -1), 1.5), ("BOTTOMPADDING", (0, 0), (-1, -1), 1.5)]
    if header:
        st.append(("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#e6eef5")))
    t.setStyle(TableStyle(st))
    return t


def fig(name, width_cm):
    from PIL import Image as PI
    w, h = PI.open(RES / name).size
    return Image(str(RES / name), width=width_cm * cm, height=width_cm * cm * h / w)


def sci(x, d=1):
    m, e = f"{x:.{d}e}".split("e")
    return f"{m}x10<super>{int(e)}</super>"


def main():
    s1, s2, sp = jl("stage1_gate.json"), jl("stage2_gate.json"), jl("stage2_sprime.json")["gate_S2prime"]
    des = jl("stage2_design.json")
    g3, g3a = jl("stage3_gate.json"), jl("stage3_gate_amended.json")
    du, dc = jl("stage3_decision_summary.json"), jl("stage3_cal_decision_summary.json")
    obs = jl("stage3_observability.json")
    tor = jl("stage3_tornado_summary.json")
    gates = cl("gate_table.csv")
    margins = cl("stage2_margins.csv")
    req = req_table()
    reqrows = cl("stage2_requirements.csv")
    nc = cl("stage2_noncollinear.csv")
    snap = cl("stage3_exploratory_snapshot.csv")
    kap = cl("stage3_kappa_exploratory_corrected.csv")
    exp_t = cl("stage3_exploratory_tuned.csv")
    F = cl("stage2_failure_rates.csv")
    orth_fail = np.mean([float(r["fail_rate_gt_10deg"]) for r in F if r["set"] == "K3orth" and r["eps"] == "1.0"])
    dmap = cl("stage3_cal_decision_map.csv")
    calpass = [r for r in dmap if r["criterion_met"] == "True"]
    imp = cl("stage3_cal_B_improvement.csv")

    doc = SimpleDocTemplate(str(ROOT / "report.pdf"), pagesize=A4, leftMargin=1.6 * cm, rightMargin=1.6 * cm, topMargin=1.4 * cm, bottomMargin=1.4 * cm,
                            title="Radical-pair compass: navigation-relevance feasibility study", author="Ujesha Sivaramakrishnan")
    W = 17.8 * cm
    S = []
    # ------------------------------------------------------------------ page 1
    S += [P("Radical-pair compass: a simulation-only navigation-relevance feasibility study", TITLE),
          P("<b>Question.</b> Under what conditions, if any, could a radical-pair (spin-chemical) magnetic compass be useful as an auxiliary heading sensor for GNSS-denied inertial navigation, compared with a conventional MEMS magnetometer? "
            "A \"no, never\" answer was declared a valid outcome in advance.", BODY),
          P("<b>Scope.</b> Everything here is simulation. No experiment exists and no new physics is claimed. Parameters marked ASSUMED are neither measured nor verified; they are swept. "
            "The attitude (not position) filter is the navigation object studied. Pre-registration: <font face='Courier'>PLAN.md</font> (committed before any run) plus dated amendments; every deviation is in <font face='Courier'>CHANGELOG.md</font>.", BODY),
          P("What is and is not new", H2),
          P("<b>Not new:</b> radical-pair compass precision limits and Fisher-information analyses (Ren, Hiscock &amp; Hore 2021; Procopio &amp; Ritz 2013, 2020; Cai 2011; Bezchastnov &amp; Domratcheva 2026), cryptochrome magnetic sensitivity (Xu et al. 2021), "
            "the ChemNav programme, and flight-tested quantum magnetic navigation (arXiv 2504.08167, anomaly matching, not compassing). A web search found no study that models a radical-pair sensor as an aiding measurement in an inertial filter; "
            "that is absence of evidence from a limited search (one patent, US12078486, was read: low overlap). <b>New (the only claim):</b> an open, reproducible requirements-and-feasibility framework that connects a spin-yield model to a photon-limited estimator and an attitude-filter comparison against a MEMS magnetometer. Unverified citations are listed in <font face='Courier'>open_questions.md</font>.", BODY),
          P("Summary of findings", H2)]
    pre_n, pre_t = du["n_passing"], du["n_configs"]
    S += [P(f"1. <b>Spin model.</b> With verified literature tensors (Hiscock 2016), the minimal two-nitrogen model has a yield anisotropy of {s1['primary']['aniso']:.3f}; assumed richer models give about 0.01 and relaxation at 10<super>7</super>/s leaves 1.5x the 10<super>-3</super> threshold. The minimal model is an optimistic upper bound.", BODY),
          P(f"2. <b>Sensor.</b> The yield is symmetric under field reversal, so only the field <i>axis</i> is observable. Because both nitrogen tensors share one principal frame the molecule has exact sign-flip symmetries: element sets related by axis permutations have a four-fold mirror ambiguity (failure rate {orth_fail:.2f}, pre-registered Gate S2 FAILED). A non-orthogonal set chosen by a pre-registered rule removes it (Gate S2' passed on fresh data).", BODY),
          P(f"3. <b>Navigation, pre-registered:</b> {pre_n} of {pre_t} configurations reduce yaw RMSE by at least 20% at both 60 s and 600 s relative to a MEMS magnetometer: <b>answer NO</b>. The baseline filter was badly inconsistent (unmodelled soft iron), and Gate S3(i) FAILED.", BODY),
          P(f"4. <b>Navigation, with online calibration (Amendment 3, specified after the NO):</b> filters become consistent; {dc['n_passing']} of {dc['n_configs']} configurations meet the criterion, mostly in two corners (a noisy 2 uT magnetometer against a 0.3 deg/axis radical-pair sensor, and 15 deg dip; one D row sits at a high-bias, low-ARW gyro). At 0.5 uT and 45 deg dip none do.", BODY),
          P("5. <b>Photons are not the binding constraint</b> for a 0.3 deg/axis sensor (see requirements); the shared magnetic disturbance, the baseline's calibration and the sensor-error model are.", BODY)]
    # ------------------------------------------------------------------ models + stage 1
    S += [P("1. Spin-dynamics model (Stage 1)", H1),
          P("Radical A carries two spin-1 nuclei (N5, N10), radical B none; B<sub>0</sub> = 50 uT, singlet-born, Haberkorn recombination k<sub>S</sub> = k<sub>T</sub> = k = 10<super>6</super>/s. "
            "Yield: &Phi;<sub>S</sub>(n) = (1/M) &Sigma; |&lt;m|P<sub>S</sub>|n&gt;|<super>2</super> k<super>2</super>/(k<super>2</super> + (E<sub>m</sub> - E<sub>n</sub>)<super>2</super>). "
            "<b>Tensors:</b> the values recalled in the study brief did not match the primary source; the verified values (N5 = (-0.087, -0.100, 1.757) mT, N10 = (-0.014, -0.024, 0.605) mT, Hiscock et al. PNAS 2016) are primary, the recalled set is a labelled sensitivity set. "
            "<b>Validation (Gate S1 PASSED):</b> isotropic hyperfine gives a direction-independent yield (spread 1e-15); H = 0 gives yield 1; yields in [0,1]; &Phi;(n) = &Phi;(-n); an independent Liouville-space implementation agrees to 1e-15 on 200 directions; grid and surrogate converge. "
            f"Primary anisotropy &Delta;&Phi; = {s1['primary']['aniso']:.4f} (&Phi; from {s1['primary']['phi_min']:.3f} to {s1['primary']['phi_max']:.3f}), against the pre-registered threshold 10<super>-3</super>.", BODY),
          fig("fig_stage1_yield_map.png", 15.5)]
    mrows = [["Model (50 uT, k = 1e6/s)", "&Delta;&Phi;", "margin vs 1e-3"]]
    pick = ["minimal H16", "R1 a_iso_B=0.0 d_B=0.1", "R1 a_iso_B=1.5 d_B=0.5", "R2", "R3", "EXPLORATORY exchange J_mT=0.5", "EXPLORATORY dipolar D_mT=0.3",
            "EXPLORATORY relaxation (meridian scan) relax_rate=100000.0", "EXPLORATORY relaxation (meridian scan) relax_rate=1000000.0", "EXPLORATORY relaxation (meridian scan) relax_rate=10000000.0"]
    for r in margins:
        if any(r["model"].startswith(p) for p in pick):
            mrows.append([r["model"].replace("EXPLORATORY ", "(exploratory) ").replace(" (meridian scan)", ""), f"{float(r['aniso']):.4f}", f"{float(r['margin_vs_1e-3']):.1f}x"])
    S += [table(mrows, [10.6 * cm, 2.6 * cm, 3.0 * cm]),
          P("R1-R3 and the exchange, dipolar and relaxation terms are ASSUMED parameters (no verified Trp-side tensors were available); relaxation values come from a meridian scan.", SMALL)]
    # ------------------------------------------------------------------ stage 2
    S += [P("2. Photon-limited sensor and estimator (Stage 2)", H1),
          P("K identical oriented elements at rotations R<sub>j</sub>; counts ~ Poisson(N(1 + &epsilon;(&Phi;<sub>j</sub> - &Phi;&#772;)/&Phi;&#772;)); &epsilon; = 1 is an ideal ordered ensemble, &epsilon; &lt; 1 absorbs partial order and incomplete yield transfer (ASSUMED). "
            "Estimator: Poisson maximum likelihood over the field axis (hemisphere grid plus multi-start Fisher scoring); bound: analytic Cramer-Rao in the tangent plane.", BODY),
          P(f"<b>Gate S2 (pre-registered, K3orth) FAILED</b>: ratios of Monte-Carlo to bound {', '.join(f'{r:,.0f}' for r in s2['K3orth']['large_N_ratios'])}, slope {s2['K3orth']['slope']:.3f}. "
            "Cause: exact sign-flip symmetry (the likelihood has equal-height mirror modes), so a local bound is ill-posed. <b>Amendment 1</b> (post hoc, committed before any new run) selected a K=3 set by maximising the minimum Poisson KL divergence between axes more than 10 deg apart (design seeds only), then evaluated it on fresh data: "
            f"<b>Gate S2' PASSED</b> (ratios {', '.join(f'{r:.3f}' for r in sp['large_N_ratios'])}, slope {sp['slope']:.3f}). Amendment 1.1 fixed a blind spot of the first objective (found by a smoke test before any design run); a single-start estimator bug was found after S2' passed and fixed, with all results re-run (CHANGELOG).", BODY),
          fig("fig_stage2_error_vs_N.png", 8.2), ]
    ncr = {(r["model"], r["alpha_deg"]): r for r in nc}
    S += [P(f"<b>Does the ambiguity need a shared principal frame?</b> Yes. Tilting the N10 frame by 5 deg makes mirror images separable (D = 10 nats at the weakest pair) at {sci(float(ncr[('NC','5')]['N_mirror_photons']))} photons in the minimal model; in a richer assumed model with a 15 deg tilt the ambiguity survives at 10<super>8</super> photons (needs {sci(float(ncr[('NC-R2','15')]['N_mirror_photons']))}). Tilts are ASSUMED.", BODY)]
    # requirements
    rr = [["Contrast row", "&epsilon;", "N for 5 deg (median dir)", "N for 5 deg (worst dir)", "N for failures &le; 5%", "N for 0.3 deg/axis (A2.2)"]]
    for en, lab in (("1", "ideal"), ("0.1", "0.1"), ("eps_R2", "richer model (derived)"), ("eps_relax_1e+06", "relaxation 1e6 (derived)"), ("eps_relax_1e+07", "relaxation 1e7 (derived)")):
        r = [x for x in reqrows if x["set"] == "K3selected" and x["eps_name"] == en and x["realism"] == "baseline"][0]
        n03, _ = n_required(0.3, req[en])
        rr.append([lab, f"{float(r['eps']):.3g}", sci(float(r["N_5deg_median_dir"])), sci(float(r["N_5deg_worst_dir"])), sci(req[en]["N_fail"]), sci(n03)])
    S += [table(rr, [3.9 * cm, 1.3 * cm, 3.1 * cm, 3.1 * cm, 3.1 * cm, 3.3 * cm]),
          P("Photons per element per measurement, selected K=3 set. Median-direction values use the Cramer-Rao median error; the failure column is the smallest N with at most 5% of trials more than 10 deg off (300 trials, half-decade grid); "
            "worst-direction values show that near-singular field directions need orders of magnitude more photons. Derived contrasts are ratios of relative anisotropies from Stage 1 (approximate). "
            "At the Cramer-Rao N for a 5 deg median error, 24-26% of Monte-Carlo trials still fail by more than 10 deg (results/stage2_requirements_mc.csv): the bound alone is optimistic at low contrast.", SMALL),
          P("Comparison with Ren et al. (2021): with their phenomenological signal and K=1000 random orientations our bound gives 6.8x10<super>7</super> photons (default) against their 5.3x10<super>7</super>-2.3x10<super>8</super>, and 4.6x10<super>6</super> (optimised) against 3.94x10<super>5</super>: a factor of about 12 remains and is not tuned away. "
            "With only K=3 elements the gap was ~800x; element count explained most of it.", BODY),
          fig("fig_stage2_failure_rates.png", 12.0)]
    # ------------------------------------------------------------------ stage 3 set up + gates
    S += [PageBreak(), P("3. Attitude-filter experiment (Stage 3)", H1),
          P("Multiplicative EKF (quaternion, nav-frame error state, 6 states: attitude error and gyro bias) at 10 Hz on a 100 Hz truth; 600 s scripted trajectory (3-5 deg tilt oscillation, three slow +/-90 deg yaw turns, then a figure-eight, tilt up to 7 deg); "
            "gyro bias instability {1, 10, 100} deg/h, ARW {0.1, 0.3, 1} deg/sqrt(h); accelerometer 0.05 m/s<super>2</super> plus 0.3 m/s<super>2</super> white vibration; Earth field 50 uT at dip {15, 45, 66} deg. "
            "<b>Arms:</b> A gyro + accelerometer; B + MEMS magnetometer (full vector, noise {0.1, 0.5, 2} uT); C + radical-pair axis (random sign, isotropic tangent error sigma_rp {0.3, 1, 2, 5, 10} deg, sign-invariant residual); D = B + C. "
            "<b>Disturbance, identical for both sensors</b> (no immunity is credited): 1% soft iron, hard iron and motor-current-correlated fields with nominal RMS sigma_d {0, 2, 10} uT. 200 paired runs per configuration (common random numbers). "
            "Decision criterion (pre-registered): yaw-RMSE reduction of at least 20% at both 60 s and 600 s with a paired-bootstrap 95% lower bound above 0.", BODY),
          P("Validation gates", H2)]
    gt = [["Gate", "Status", "Detail"]] + [[r["gate"], "<b>" + r["status"] + "</b>", r["detail"]] for r in gates]
    S += [table(gt, [4.6 * cm, 1.9 * cm, 11.3 * cm]),
          P("Gate S2 (K3orth) and S3 stay FAILED. S2' and S3(i)' were defined after the failures with the owner's approval; S3(i)' is post hoc and its outcome was visible beforehand, so it is a confirmation of the diagnosis (unmodelled figure-eight acceleration absorbed by the bias states), not an independent test.", SMALL)]
    # ------------------------------------------------------------------ results
    S += [P("4. Results", H1),
          P(f"<b>Pre-registered answer: NO.</b> {pre_n} of {pre_t} C/D configurations meet the criterion. At 60 s no configuration exceeds 3.5%; at 600 s alone 20 pass, all of them artefacts: in the main grid the MEMS filter is strongly over-confident "
            "(NEES 183-480 against about 6 even with no disturbance) because soft iron is not in its noise model, so a sensor that trusts the same biased measurement less appears to win. "
            "Any magnetic aiding helps hugely (yaw RMSE at 600 s: 10.8 deg gravity-only against 0.18-0.57 deg), and the radical-pair channel adds almost nothing on top of it under shared disturbance.", BODY),
          fig("fig_stage3_yaw_rmse.png", 15.5)]
    e1 = [r for r in exp_t if r["exp"] == "E1_tuned_R"]; e2 = [r for r in exp_t if r["exp"] == "E2_no_soft_iron"]
    S += [P(f"<b>Exploratory E1/E2.</b> Adding the soft-iron variance to R did not restore consistency (NEES still 118-462) because the error is a constant body-frame bias: criterion met in {sum(r['criterion_met']=='True' for r in e1)} of {len(e1)}. "
            f"With no soft iron and no disturbance the filter is consistent, and a sensor twice as accurate as the MEMS one (0.3 deg against about 0.57 deg) improves yaw by only 3% (60 s) and 11% (600 s): criterion met in {sum(r['criterion_met']=='True' for r in e2)} of {len(e2)}.", BODY)]
    # calibrated
    S += [P("Online calibration (Amendment 3)", H2),
          P("Both sensors get the same treatment: 18-state filters that estimate hard iron (3), soft iron (6, symmetric) and the motor-current field per unit throttle (3, throttle profile known), with priors consistent with the truth generator and no inflation of R. "
            f"Filters become consistent (NEES about 6-7 even at sigma_d = 10 uT) and the magnetometer baseline improves (e.g. 600 s yaw RMSE {float([r for r in imp if float(r['bias_dph'])==10 and float(r['arw'])==0.3 and float(r['sigma_d'])==2 and float(r['dip'])==45 and float(r['sigma_mag'])==0.5][0]['rmse600_B_uncal_deg']):.2f} deg to {float([r for r in imp if float(r['bias_dph'])==10 and float(r['arw'])==0.3 and float(r['sigma_d'])==2 and float(r['dip'])==45 and float(r['sigma_mag'])==0.5][0]['rmse600_B_cal_deg']):.2f} deg). "
            f"<b>{dc['n_passing']} of {dc['n_configs']} configurations now meet the criterion</b>:", BODY)]
    pr = [["Arm", "MEMS noise (uT)", "dip", "sigma_d", "sigma_rp", "60 s: reduction (CI low)", "600 s: reduction (CI low)"]]
    for r in calpass:
        pr.append([r["arm"], r["sigma_mag"], f"{float(r['dip']):g}", r["sigma_d"], r["sigma_rp_deg"], f"{float(r['reduction_60'])*100:.0f}% ({float(r['ci_lo_60'])*100:.0f}%)", f"{float(r['reduction_600'])*100:.0f}% ({float(r['ci_lo_600'])*100:.0f}%)"])
    S += [table(pr, [1.0 * cm, 2.4 * cm, 1.2 * cm, 1.8 * cm, 1.8 * cm, 4.8 * cm, 4.8 * cm]),
          P("Mid gyro (10 deg/h, 0.3 deg/sqrt(h)) unless gyro corners are shown in the CSV. All passing rows have sigma_rp = 0.3 deg (the smallest grid value) except one D row at 1 deg, so the benefit at smaller sigma_rp is not known. "
            "At the base setting (0.5 uT, 45 deg) arm C (0.3 deg) changes yaw RMSE by -9% (60 s) and +10% (600 s) at sigma_d = 2 uT (-31% and +10% at sigma_d = 0) relative to the calibrated magnetometer: calibration of an axis-only measurement converges more slowly. The 12 passes cluster in structured regions with lower bounds of at least 11%, so they are unlikely to be multiplicity noise, but 270 configurations were tested.", SMALL),
          fig("fig_stage3_cal_yaw_rmse.png", 15.5)]
    # observability, snapshot, kappa
    ob = [r for r in obs if r["window"] == "hover_0_120s" and r["dip"] == 45.0]
    ow = {r["arm"]: r for r in ob}
    sn = [r for r in snap if float(r["dip"]) == 15.0 and float(r["sigma_rp_deg"]) == 10.0 and float(r["sigma_d"]) == 0.0][0]
    k5 = [r for r in kap if float(r["kappa"]) == 0.5 and float(r["sigma_d"]) == 10.0 and r["arm"] == "C"]
    S += [P(f"<b>Observability.</b> Gravity-only has numerical rank 5 (the single null direction is yaw); arms B, C and D have rank 6 in hover and turn windows at dips 15, 45 and 66 deg, so the axis-only channel loses nothing locally (information-implied yaw std at 45 deg hover: B {ow['B']['yaw_std_deg_from_information']:.3f}, C {ow['C']['yaw_std_deg_from_information']:.3f}, D {ow['D']['yaw_std_deg_from_information']:.3f} deg). "
            f"<b>Global sign ambiguity</b> (exploratory snapshot initialisation): the axis-only sensor flips yaw by 180 deg in at most {float(sn['frac_yaw_gt_90deg_axis_only'])*100:.1f}% of draws (dip 15, sigma_rp 10 deg, no disturbance) and essentially never at 45 deg dip. "
            "<b>Sensor-specific advantage kappa</b> (exploratory, corrected model): the radical-pair channel helps only if it were assumed to respond to the disturbance much less than a MEMS magnetometer; nothing in this study supports that. "
            f"(kappa = 0.5, sigma_d = 10 uT, arm C: criterion met in {sum(r['criterion_met']=='True' for r in k5)} of {len(k5)}.)", BODY)]
    # tornado
    S += [P("5. Sensitivity (Stage 4)", H1),
          P("Every ASSUMED parameter was varied one at a time around the base (mid gyro, dip 45 deg, 0.5 uT, sigma_d = 2 uT) for both filter types; the metric is the best-case benefit m = max over sigma_rp of min(reduction at 60 s, reduction at 600 s) of arm C relative to B. "
            f"Base m: pre-registered filters {tor['uncal']['base_metric']*100:.0f}%, calibrated {tor['cal']['base_metric']*100:.0f}%.", BODY),
          fig("fig_stage4_tornado.png", 17.0)]
    trk = [["Rank", "Parameter (calibrated)", "swing of m (points)", "reaches 20%?", "Parameter (pre-registered)", "swing (points)", "reaches 20%?"]]
    for i in range(min(8, len(tor["cal"]["ranking"]))):
        c, u = tor["cal"]["ranking"][i], tor["uncal"]["ranking"][i]
        trk.append([i + 1, c["parameter"], f"{c['swing']*100:.0f}", "yes" if c["reaches_20pct"] else "no", u["parameter"], f"{u['swing']*100:.0f}", "yes" if u["reaches_20pct"] else "no"])
    S += [table(trk, [0.9 * cm, 4.6 * cm, 2.0 * cm, 1.8 * cm, 4.6 * cm, 1.8 * cm, 1.8 * cm])]
    # limits
    S += [P("6. Limits, what was not modelled, and what was assumed", H1),
          P("<b>Not modelled:</b> position and velocity (no GNSS-denied position drift); real magnetic environments (gradients, anomalies, motor EMI spectrum, time-varying external fields); chemistry kinetics, temperature (the closest experiment, Kerpal 2019, was at 120 K), light level, photobleaching, ageing; "
            "orientational disorder inside the spin model (only through epsilon); hyperfine fluctuations; direction-dependent radical-pair covariance (isotropic Gaussian used) and ambiguity outliers (up to 5% not simulated); correlation of the two sensors' errors in the uncalibrated arm D; "
            "accelerometer bias and scale errors, gyro scale factor and misalignment, coning, timing and latency; hardware of any kind.", BODY),
          P("<b>ASSUMED parameters (all swept or fixed as stated in PLAN.md):</b> hyperfine principal-axis alignment and all Trp-side couplings; exchange, dipolar and relaxation values; contrast epsilon and its derived values; lock-in efficiency 2/pi and 10 Hz updates; background, electronics noise and gain drift; photon caps; "
            "gyro bias/ARW and the 10% bias random walk; accelerometer noise and vibration; initial errors (1, 1, 10 deg); Earth field 50 uT and dip; MEMS noise; soft-iron 1%; disturbance split, throttle profile and nominal RMS; radical-pair per-axis error sigma_rp; kappa = 1; calibration priors; trajectory.", BODY),
          P("<b>Weakest points:</b> (1) no measurement anywhere: epsilon and sigma_rp are assumptions; (2) the calibrated baseline is an idealised EKF with exact priors and known throttle profile, which flatters both sensors; (3) a sequence of post-hoc amendments after failed gates, even though each is dated, committed first and keeps the failure on record. "
            "Reproduction: <font face='Courier'>README.md</font>; seeds from SeedSequence(20261013); Python 3.14, numpy 2.4, scipy 1.17.", BODY)]
    doc.build(S)
    from pypdf import PdfReader
    print("pages", len(PdfReader(str(ROOT / "report.pdf")).pages))


if __name__ == "__main__":
    main()
