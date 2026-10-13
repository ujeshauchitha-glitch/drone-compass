"""Analysis of the calibrated decision study (Amendment 3, A3.3).  Same statistic and criterion as A2.3, baseline B_cal."""
from __future__ import annotations

import csv
import gzip
import json
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from experiments.stage3_analysis import CAPS, EPS_ROWS, M, N_BOOT, n_required, req_table  # noqa: E402

RES = ROOT / "results"


def fnum(x):
    try:
        return float(x)
    except ValueError:
        return float("nan")


def load(prefix):
    rows = list(csv.DictReader(open(RES / f"{prefix}_summary.csv")))
    for r in rows:
        for k in list(r):
            if k not in ("setting", "arm", "group"):
                r[k] = fnum(r[k])
        r["id"] = int(r["id"])
    errs = defaultdict(dict)
    with gzip.open(RES / f"{prefix}_yaw_errors.csv.gz", "rt") as fh:
        rd = csv.reader(fh); next(rd)
        for cid, t, run, v in rd:
            errs[int(cid)].setdefault(int(t), []).append(float(v))
    return rows, {k: {t: np.array(v) for t, v in d.items()} for k, d in errs.items()}


def main():
    rows, errs = load("stage3_cal")
    urows, uerrs = load("stage3")
    req = req_table()
    rng = np.random.default_rng(np.random.SeedSequence(20261013, spawn_key=(11,)))
    idx = rng.integers(0, M, (N_BOOT, M))
    wrap = lambda x: (x + np.pi) % (2 * np.pi) - np.pi

    def red(eb, ex):
        eb, ex = wrap(eb), wrap(ex)
        rb, rx = np.sqrt(np.mean(eb**2)), np.sqrt(np.mean(ex**2))
        bb, bx = np.sqrt(np.mean(eb[idx] ** 2, 1)), np.sqrt(np.mean(ex[idx] ** 2, 1))
        r = (bb - bx) / bb
        return rb, rx, (rb - rx) / rb, float(np.percentile(r, 2.5))

    key = lambda r: (r["bias_dph"], r["arw"], r["sigma_d"], r["dip"], r["sigma_mag"])
    B = {key(r): r for r in rows if r["arm"] == "B"}
    UB = {key(r): r for r in urows if r["arm"] == "B"}
    dec, improve = [], []
    for k_, b in B.items():
        u = UB[k_]
        improve.append(dict(bias_dph=k_[0], arw=k_[1], sigma_d=k_[2], dip=k_[3], sigma_mag=k_[4],
                            rmse600_B_uncal_deg=np.degrees(u["rmse_600"]), rmse600_B_cal_deg=np.degrees(b["rmse_600"]),
                            rmse60_B_uncal_deg=np.degrees(u["rmse_60"]), rmse60_B_cal_deg=np.degrees(b["rmse_60"]),
                            nees600_B_cal=b["nees_600"]))
    for r in rows:
        if r["arm"] not in "CD":
            continue
        b = B[key(r)]
        s60, s600 = red(errs[b["id"]][60], errs[r["id"]][60]), red(errs[b["id"]][600], errs[r["id"]][600])
        row = dict(setting=r["setting"], arm=r["arm"], bias_dph=r["bias_dph"], arw=r["arw"], sigma_d=r["sigma_d"], dip=r["dip"],
                   sigma_mag=r["sigma_mag"], sigma_rp_deg=r["sigma_rp_deg"],
                   rmse_B_60_deg=np.degrees(s60[0]), rmse_X_60_deg=np.degrees(s60[1]), reduction_60=s60[2], ci_lo_60=s60[3],
                   rmse_B_600_deg=np.degrees(s600[0]), rmse_X_600_deg=np.degrees(s600[1]), reduction_600=s600[2], ci_lo_600=s600[3],
                   criterion_met=bool(all(s[2] >= 0.2 and s[3] > 0 for s in (s60, s600))), nees_600=r["nees_600"])
        for en in EPS_ROWS:
            n, _ = n_required(r["sigma_rp_deg"], req[en])
            row[f"N_req_{en}"] = n
            for cap in CAPS:
                row[f"inside_{en}_cap{cap:.0e}"] = bool(n <= cap)
        dec.append(row)
    with open(RES / "stage3_cal_decision_map.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(dec[0].keys())); w.writeheader(); w.writerows(dec)
    with open(RES / "stage3_cal_B_improvement.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(improve[0].keys())); w.writeheader(); w.writerows(improve)

    passing = [d for d in dec if d["criterion_met"]]
    summary = dict(n_configs=len(dec), n_passing=len(passing), n_boot=N_BOOT,
                   passing=[(d["arm"], d["bias_dph"], d["arw"], d["sigma_d"], d["dip"], d["sigma_mag"], d["sigma_rp_deg"]) for d in passing])
    table = {}
    for en in EPS_ROWS:
        for cap in CAPS:
            ok = [d for d in passing if d[f"inside_{en}_cap{cap:.0e}"]]
            table[f"{en}|cap={cap:.0e}"] = dict(any_arm_C_or_D=len(ok) > 0, n_passing_inside=len(ok))
    summary["decision_by_epsilon_row_and_cap"] = table
    summary["literal_overall_answer_yes"] = len(passing) > 0
    summary["max_reduction_60_any"] = float(max(d["reduction_60"] for d in dec))
    summary["max_reduction_600_any"] = float(max(d["reduction_600"] for d in dec))
    summary["n_meeting_600_only"] = int(sum(d["reduction_600"] >= 0.2 and d["ci_lo_600"] > 0 for d in dec))
    summary["n_meeting_60_only"] = int(sum(d["reduction_60"] >= 0.2 and d["ci_lo_60"] > 0 for d in dec))
    # mid-gyro base headline
    mid = lambda r: r["bias_dph"] == 10 and r["arw"] == 0.3 and r["dip"] == 45 and r["sigma_mag"] == 0.5
    head = {}
    for sd in (0.0, 2.0, 10.0):
        b = [r for r in rows if r["arm"] == "B" and mid(r) and r["sigma_d"] == sd][0]
        ub = [r for r in urows if r["arm"] == "B" and mid(r) and r["sigma_d"] == sd][0]
        head[f"sigma_d={sd}"] = dict(
            B_cal_60_deg=float(np.degrees(b["rmse_60"])), B_cal_600_deg=float(np.degrees(b["rmse_600"])),
            B_uncal_600_deg=float(np.degrees(ub["rmse_600"])), nees_B_cal=[b["nees_60"], b["nees_300"], b["nees_600"]],
            C_cal_600_deg={str(r["sigma_rp_deg"]): float(np.degrees(r["rmse_600"])) for r in rows if r["arm"] == "C" and mid(r) and r["sigma_d"] == sd},
            D_cal_600_deg={str(r["sigma_rp_deg"]): float(np.degrees(r["rmse_600"])) for r in rows if r["arm"] == "D" and mid(r) and r["sigma_d"] == sd},
            nees_C_cal_600={str(r["sigma_rp_deg"]): r["nees_600"] for r in rows if r["arm"] == "C" and mid(r) and r["sigma_d"] == sd})
    summary["headline_mid_gyro_base_config"] = head
    # break-even
    be = []
    groups = defaultdict(list)
    for d in dec:
        if d["arm"] == "C":
            groups[(d["bias_dph"], d["arw"], d["sigma_d"], d["dip"], d["sigma_mag"])].append(d)
    for k_, lst in groups.items():
        lst = sorted(lst, key=lambda x: x["sigma_rp_deg"]); s = np.array([x["sigma_rp_deg"] for x in lst])
        out = dict(bias_dph=k_[0], arw=k_[1], sigma_d=k_[2], dip=k_[3], sigma_mag=k_[4])
        for t in (60, 600):
            ratio = np.array([x[f"rmse_X_{t}_deg"] / x[f"rmse_B_{t}_deg"] for x in lst])
            if np.all(ratio > 1):
                out[f"break_even_t{t}"] = "C worse at every grid sigma_rp"
            elif np.all(ratio < 1):
                out[f"break_even_t{t}"] = "C better at every grid sigma_rp"
            else:
                j = int(np.where(np.diff(np.sign(ratio - 1)) != 0)[0][0]); f = (1 - ratio[j]) / (ratio[j + 1] - ratio[j])
                out[f"break_even_t{t}"] = float(np.exp(np.log(s[j]) + f * (np.log(s[j + 1]) - np.log(s[j]))))
            out[f"ratio_C_over_B_t{t}"] = [round(float(x), 4) for x in ratio]
        be.append(out)
    with open(RES / "stage3_cal_break_even.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(be[0].keys())); w.writeheader(); w.writerows(be)
    json.dump(summary, open(RES / "stage3_cal_decision_summary.json", "w"), indent=1)

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    BLUE, ORANGE, GREEN, GREY = "#0072B2", "#D55E00", "#009E73", "#6b6b6b"
    fig, axs = plt.subplots(1, 3, figsize=(11, 3.8), sharey=True)
    srps = (0.3, 1.0, 2.0, 5.0, 10.0)
    for ax, sd in zip(axs, (0.0, 2.0, 10.0)):
        h = head[f"sigma_d={sd}"]
        ax.axhline(h["B_uncal_600_deg"], color=GREY, ls=":", label="B, no calibration (pre-registered)")
        ax.axhline(h["B_cal_600_deg"], color=GREEN, ls="--", label="B, online calibration")
        ax.loglog(srps, [h["C_cal_600_deg"][str(s)] for s in srps], "o-", color=BLUE, label="C, radical pair, calibrated")
        ax.loglog(srps, [h["D_cal_600_deg"][str(s)] for s in srps], "s-", color=ORANGE, label="D, B + C, calibrated")
        ax.set_title(f"disturbance σ_d = {sd:g} µT"); ax.set_xlabel("radical-pair σ_rp per axis (deg)"); ax.grid(alpha=.3, which="both")
    axs[0].set_ylabel("yaw RMSE at t = 600 s (deg)"); axs[0].legend(fontsize=7)
    fig.suptitle("Calibrated filters, mid gyro, dip 45°, MEMS σ = 0.5 µT; 200 paired runs", fontsize=9)
    fig.tight_layout(); fig.savefig(RES / "fig_stage3_cal_yaw_rmse.png", dpi=150); plt.close(fig)
    print(json.dumps(dict(n_passing=summary["n_passing"], n_configs=summary["n_configs"], n600=summary["n_meeting_600_only"], n60=summary["n_meeting_60_only"],
                          max60=summary["max_reduction_60_any"], max600=summary["max_reduction_600_any"]), indent=1))


if __name__ == "__main__":
    main()
