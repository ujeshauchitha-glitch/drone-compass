"""Stage 3 decision analysis (PLAN.md 5.4 + Amendment 2 A2.2-A2.3).  Everything is computed from results/ files."""
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
RES = ROOT / "results"
M = 200
N_BOOT = 2000
CAPS = (1e8, 1e9, 1e10)
EPS_ROWS = ["1", "0.1", "eps_R2", "eps_relax_1e+06", "eps_relax_1e+07"]
K_MED = np.sqrt(2 * np.log(2))


def fnum(x):
    try:
        return float(x)
    except ValueError:
        return float("nan")


def load():
    rows = list(csv.DictReader(open(RES / "stage3_summary.csv")))
    for r in rows:
        for k in ("id", "bias_dph", "arw", "sigma_d", "dip", "sigma_mag", "sigma_rp_deg", "kappa", "rmse_60", "rmse_300", "rmse_600",
                  "nees_60", "nees_300", "nees_600", "tilt_rmse_final"):
            r[k] = fnum(r[k])
        r["id"] = int(r["id"])
    errs = defaultdict(dict)
    with gzip.open(RES / "stage3_yaw_errors.csv.gz", "rt") as fh:
        rd = csv.reader(fh); next(rd)
        for cid, t, run, v in rd:
            errs[int(cid)].setdefault(int(t), []).append(float(v))
    errs = {k: {t: np.array(v) for t, v in d.items()} for k, d in errs.items()}
    return rows, errs


def req_table():
    R = list(csv.DictReader(open(RES / "stage2_requirements.csv")))
    fail = {t["eps_name"]: t["N_fail_le_5pct"] for t in json.load(open(RES / "stage2_failure_threshold.json"))}
    out = {}
    for en in EPS_ROWS:
        r = [x for x in R if x["set"] == "K3selected" and x["eps_name"] == en and x["realism"] == "baseline"][0]
        out[en] = dict(sigma_med_ref=float(r["sigma_med_deg_at_1e8_median_dir"]), N_fail=float(fail[en]), eps=float(r["eps"]))
    return out


def n_required(srp_deg, row):
    n_crb = 1e8 * (row["sigma_med_ref"] / (K_MED * srp_deg)) ** 2
    return max(n_crb, row["N_fail"]), n_crb


def main():
    rows, errs = load()
    req = req_table()
    rng = np.random.default_rng(np.random.SeedSequence(20261013, spawn_key=(11,)))
    idx = rng.integers(0, M, (N_BOOT, M))
    Bidx = {}
    for r in rows:
        if r["arm"] == "B":
            Bidx[(r["bias_dph"], r["arw"], r["sigma_d"], r["dip"], r["sigma_mag"])] = r

    def reduction(rb, rx):
        out = {}
        for t in (60, 600):
            eb, ex = errs[rb["id"]][t], errs[rx["id"]][t]
            eb = (eb + np.pi) % (2 * np.pi) - np.pi; ex = (ex + np.pi) % (2 * np.pi) - np.pi
            rmb, rmx = np.sqrt(np.mean(eb**2)), np.sqrt(np.mean(ex**2))
            bb = np.sqrt(np.mean(eb[idx] ** 2, axis=1)); bx = np.sqrt(np.mean(ex[idx] ** 2, axis=1))
            red = (bb - bx) / bb
            out[t] = dict(rmse_B=float(rmb), rmse_X=float(rmx), reduction=float((rmb - rmx) / rmb),
                          ci_lo=float(np.percentile(red, 2.5)), ci_hi=float(np.percentile(red, 97.5)))
        return out

    dec = []
    for r in rows:
        if r["arm"] not in "CD" or r["group"] not in ("base", "oat"):
            continue
        rb = Bidx[(r["bias_dph"], r["arw"], r["sigma_d"], r["dip"], r["sigma_mag"])]
        red = reduction(rb, r)
        passed = all(red[t]["reduction"] >= 0.20 and red[t]["ci_lo"] > 0 for t in (60, 600))
        row = dict(group=r["group"], arm=r["arm"], bias_dph=r["bias_dph"], arw=r["arw"], sigma_d=r["sigma_d"], dip=r["dip"],
                   sigma_mag=r["sigma_mag"], sigma_rp_deg=r["sigma_rp_deg"],
                   rmse_B_60_deg=np.degrees(red[60]["rmse_B"]), rmse_X_60_deg=np.degrees(red[60]["rmse_X"]),
                   reduction_60=red[60]["reduction"], ci_lo_60=red[60]["ci_lo"],
                   rmse_B_600_deg=np.degrees(red[600]["rmse_B"]), rmse_X_600_deg=np.degrees(red[600]["rmse_X"]),
                   reduction_600=red[600]["reduction"], ci_lo_600=red[600]["ci_lo"], criterion_met=passed)
        for en in EPS_ROWS:
            n, ncrb = n_required(r["sigma_rp_deg"], req[en])
            row[f"N_req_{en}"] = n
            for cap in CAPS:
                row[f"inside_{en}_cap{cap:.0e}"] = bool(n <= cap)
        dec.append(row)
    with open(RES / "stage3_decision_map.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(dec[0].keys())); w.writeheader(); w.writerows(dec)

    # ------------------------------------------------------------------ decision summary
    passing = [d for d in dec if d["criterion_met"]]
    summary = dict(n_configs=len(dec), n_passing=len(passing), n_boot=N_BOOT)
    by = defaultdict(int)
    for d in passing:
        by[(d["arm"], d["sigma_d"])] += 1
    summary["passing_by_arm_and_sigma_d"] = {f"{a}|sigma_d={s}": c for (a, s), c in sorted(by.items())}
    summary["passing_sigma_rp_values"] = sorted({d["sigma_rp_deg"] for d in passing if d["arm"] == "C"})
    summary["passing_arm_C"] = [(d["group"], d["bias_dph"], d["arw"], d["sigma_d"], d["dip"], d["sigma_mag"], d["sigma_rp_deg"]) for d in passing if d["arm"] == "C"]
    table = {}
    for en in EPS_ROWS:
        for cap in CAPS:
            ok = [d for d in passing if d[f"inside_{en}_cap{cap:.0e}"]]
            table[f"{en}|cap={cap:.0e}"] = dict(any_arm_C_or_D=len(ok) > 0, any_arm_C=any(d["arm"] == "C" for d in ok), n_passing_inside=len(ok),
                                                 sigma_d_values=sorted({d["sigma_d"] for d in ok}))
    summary["decision_by_epsilon_row_and_cap"] = table
    summary["literal_overall_answer_yes"] = len(passing) > 0
    summary["literal_overall_answer_yes_arm_C_alone"] = any(d["arm"] == "C" for d in passing)

    # ------------------------------------------------------------------ break-even sigma_rp (arm C vs B)
    be_rows = []
    groups = defaultdict(list)
    for d in dec:
        if d["arm"] == "C":
            groups[(d["group"], d["bias_dph"], d["arw"], d["sigma_d"], d["dip"], d["sigma_mag"])].append(d)
    for key, lst in groups.items():
        lst = sorted(lst, key=lambda x: x["sigma_rp_deg"])
        s = np.array([x["sigma_rp_deg"] for x in lst])
        out = dict(group=key[0], bias_dph=key[1], arw=key[2], sigma_d=key[3], dip=key[4], sigma_mag=key[5])
        for t in (60, 600):
            ratio = np.array([x[f"rmse_X_{t}_deg"] / x[f"rmse_B_{t}_deg"] for x in lst])
            if np.all(ratio > 1):
                be = "C worse than B at every grid sigma_rp (>= 1)"
                val = float("nan")
            elif np.all(ratio < 1):
                be = "C better than B at every grid sigma_rp (up to 10 deg)"
                val = float("nan")
            else:
                j = int(np.where(np.diff(np.sign(ratio - 1)) != 0)[0][0])
                f = (1 - ratio[j]) / (ratio[j + 1] - ratio[j])
                val = float(np.exp(np.log(s[j]) + f * (np.log(s[j + 1]) - np.log(s[j]))))
                be = "interpolated"
            out[f"break_even_sigma_rp_deg_t{t}"] = val
            out[f"break_even_note_t{t}"] = be
            out[f"ratio_C_over_B_t{t}"] = [round(float(x), 4) for x in ratio]
        be_rows.append(out)
    with open(RES / "stage3_break_even.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(be_rows[0].keys())); w.writeheader(); w.writerows(be_rows)

    # ------------------------------------------------------------------ headline numbers for the report
    def lookup(group, arm, gyro, sd, srp=None, dip=45.0, smag=0.5, kappa=1.0):
        for r in rows:
            if (r["group"] == group and r["arm"] == arm and r["bias_dph"] == gyro[0] and r["arw"] == gyro[1] and r["sigma_d"] == sd
                    and r["dip"] == dip and r["sigma_mag"] == smag and r["kappa"] == kappa
                    and (srp is None or abs(r["sigma_rp_deg"] - srp) < 1e-9)):
                return r
    mid = (10.0, 0.3)
    head = {}
    for sd in (0.0, 2.0, 10.0):
        head[f"sigma_d={sd}"] = dict(
            A_600_deg=float(np.degrees(lookup("base", "A", mid, 0.0)["rmse_600"])),
            B_600_deg=float(np.degrees(lookup("base", "B", mid, sd)["rmse_600"])),
            B_60_deg=float(np.degrees(lookup("base", "B", mid, sd)["rmse_60"])),
            C_600_deg={str(s): float(np.degrees(lookup("base", "C", mid, sd, s)["rmse_600"])) for s in (0.3, 1.0, 2.0, 5.0, 10.0)},
            D_600_deg={str(s): float(np.degrees(lookup("base", "D", mid, sd, s)["rmse_600"])) for s in (0.3, 1.0, 2.0, 5.0, 10.0)})
    summary["headline_mid_gyro_base_config"] = head
    json.dump(summary, open(RES / "stage3_decision_summary.json", "w"), indent=1)

    # ------------------------------------------------------------------ EXPLORATORY kappa
    kp = []
    if any(r["group"] == "kappa_EXPLORATORY" for r in rows):
        for r in rows:
            if r["group"] != "kappa_EXPLORATORY":
                continue
            rb = Bidx[(r["bias_dph"], r["arw"], r["sigma_d"], 45.0, 0.5)]
            red = reduction(rb, r)
            kp.append(dict(arm=r["arm"], kappa=r["kappa"], sigma_d=r["sigma_d"], sigma_rp_deg=r["sigma_rp_deg"],
                           reduction_60=red[60]["reduction"], ci_lo_60=red[60]["ci_lo"], reduction_600=red[600]["reduction"], ci_lo_600=red[600]["ci_lo"],
                           criterion_met=all(red[t]["reduction"] >= 0.2 and red[t]["ci_lo"] > 0 for t in (60, 600))))
        with open(RES / "stage3_kappa_exploratory.csv", "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(kp[0].keys())); w.writeheader(); w.writerows(kp)

    # ------------------------------------------------------------------ figures
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    BLUE, ORANGE, GREEN, GREY = "#0072B2", "#D55E00", "#009E73", "#6b6b6b"
    fig, axs = plt.subplots(1, 3, figsize=(11, 3.8), sharey=True)
    srps = (0.3, 1.0, 2.0, 5.0, 10.0)
    for ax, sd in zip(axs, (0.0, 2.0, 10.0)):
        a = np.degrees(lookup("base", "A", mid, 0.0)["rmse_600"]); b = np.degrees(lookup("base", "B", mid, sd)["rmse_600"])
        ax.axhline(a, color=GREY, ls=":", label="A (no magnetic aiding)")
        ax.axhline(b, color=GREEN, ls="--", label="B (MEMS magnetometer)")
        ax.loglog(srps, [np.degrees(lookup("base", "C", mid, sd, s)["rmse_600"]) for s in srps], "o-", color=BLUE, label="C (radical pair, axis only)")
        ax.loglog(srps, [np.degrees(lookup("base", "D", mid, sd, s)["rmse_600"]) for s in srps], "s-", color=ORANGE, label="D (B + C)")
        ax.set_title(f"disturbance σ_d = {sd:g} µT"); ax.set_xlabel("radical-pair σ_rp per axis (deg)"); ax.grid(alpha=.3, which="both")
    axs[0].set_ylabel("yaw RMSE at t = 600 s (deg)"); axs[0].legend(fontsize=7)
    fig.suptitle("Mid gyro (10 °/h, 0.3 °/√h), dip 45°, MEMS σ = 0.5 µT; 200 paired runs", fontsize=9)
    fig.tight_layout(); fig.savefig(RES / "fig_stage3_yaw_rmse.png", dpi=150); plt.close(fig)
    print(json.dumps(dict(n_passing=summary["n_passing"], n_configs=summary["n_configs"], passing_by=summary["passing_by_arm_and_sigma_d"],
                          literal_yes=summary["literal_overall_answer_yes"], literal_yes_C=summary["literal_overall_answer_yes_arm_C_alone"]), indent=1))


if __name__ == "__main__":
    main()
