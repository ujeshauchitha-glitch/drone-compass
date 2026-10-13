"""Build results/gate_table.csv and gate_table.md from the result files (no numbers typed by hand)."""
from __future__ import annotations

import csv
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RES = ROOT / "results"


def main():
    s1 = json.load(open(RES / "stage1_gate.json"))
    s2 = json.load(open(RES / "stage2_gate.json"))
    sp = json.load(open(RES / "stage2_sprime.json"))["gate_S2prime"]
    k3 = s2["K3orth"]
    rows = [
        dict(gate="S1", configuration="minimal two-nitrogen model, verified H16 tensors", status="PASSED" if s1["gate_S1_pass"] else "FAILED",
             detail=f"checks a-f pass; primary anisotropy {s1['primary']['aniso']:.4f} vs threshold {s1['aniso_threshold']:g}", source="results/stage1_gate.json"),
        dict(gate="S2 (pre-registered)", configuration="K3orth, eps=1, Poisson only",
             status="PASSED" if (k3["gate_S2_1_pass"] and k3["gate_S2_2_pass"]) else "FAILED",
             detail=(f"condition 1 ratios {', '.join(f'{r:.4g}' for r in k3['large_N_ratios'])} (need 0.5-2); "
                     f"condition 2 slope {k3['slope']:.3f} (need -0.5 +/- 0.1)"), source="results/stage2_gate.json"),
        dict(gate="S2' (amended, Amendment 1)", configuration="K3selected (selected by the pre-registered rule), eps=1, Poisson only, fresh seeds/directions",
             status="PASSED" if (sp["cond1_pass"] and sp["cond2_pass"]) else "FAILED",
             detail=(f"condition 1 ratios {', '.join(f'{r:.3f}' for r in sp['large_N_ratios'])} (need 0.5-2); "
                     f"condition 2 slope {sp['slope']:.3f} (need -0.5 +/- 0.1); condition 3: Ren comparison documented in results/stage2_sprime.json"),
             source="results/stage2_sprime.json"),
        dict(gate="S2' diagnostic (not a gate)", configuration="K3generic_EXPLORATORY",
             status="DIAGNOSTIC", detail=(f"ratios {', '.join(f'{r:.3f}' for r in s2['K3generic_EXPLORATORY']['large_N_ratios'])}; "
                                          f"slope {s2['K3generic_EXPLORATORY']['slope']:.3f}"), source="results/stage2_gate.json"),
    ]
    s3 = json.load(open(RES / "stage3_gate.json"))
    s3a = json.load(open(RES / "stage3_gate_amended.json"))
    pre = [r for r in s3["nees_rows"] if r["condition"].startswith("PRE")]
    diag = [r for r in s3["nees_rows"] if r["condition"].startswith("DIAG")]
    fails = ", ".join(f"arm {a} at {t:g} s: {v:.2f}" for a, t, v in s3["gate_S3_i_failures"])
    iii_max = max(x["max_abs_diff_rad"] for x in s3["gate_S3_iii_rows"])
    ii_lo = min(x["ratio"] for x in s3["gate_S3_ii_rows"]); ii_hi = max(x["ratio"] for x in s3["gate_S3_ii_rows"])
    lo, hi = s3["nees_intervals"]["bonferroni_family_95"]
    rows += [
        dict(gate="S3(i) NEES (pre-registered)", configuration="nominal, mid gyro, arms A, B, C and D at sigma_rp=2 deg, t=60/300/600 s",
             status="PASSED" if s3["gate_S3_i_pass_bonferroni"] else "FAILED",
             detail=(f"outside the family-wise 95% interval [{lo:.2f}, {hi:.2f}]: {fails if fails else 'none'}; "
                     f"{sum(r['in_bonferroni_95'] for r in pre)} of {len(pre)} tests inside"), source="results/stage3_gate.json"),
        dict(gate="S3(i) diagnostic (not a gate)", configuration="same, but with the linear acceleration of the figure-eight removed",
             status="DIAGNOSTIC", detail=f"{sum(r['in_bonferroni_95'] for r in diag)} of {len(diag)} tests inside the same interval", source="results/stage3_gate.json"),
        dict(gate="S3(i)' (amended, POST HOC)", configuration="same 12 tests on trajectory T0 (linear acceleration = 0); outcome was visible in the diagnostic before the amendment",
             status="PASSED" if s3a["gate_S3_i_prime_pass_bonferroni"] else "FAILED",
             detail=f"{sum(r['in_bonferroni_95'] for r in s3a['rows'])} of {len(s3a['rows'])} tests inside [{lo:.2f}, {hi:.2f}]; not an independent test", source="results/stage3_gate_amended.json"),
        dict(gate="S3(ii) arm A vs analytic drift", configuration="nine gyro settings, t=60 and 600 s", status="PASSED" if s3["gate_S3_ii_pass"] else "FAILED",
             detail=f"ratio simulated/analytic {ii_lo:.3f} to {ii_hi:.3f} (need within a factor 2)", source="results/stage3_gate.json"),
        dict(gate="S3(iii) independent filter", configuration="rotation matrices, numerical Jacobians, 5 runs x arms A-D", status="PASSED" if s3["gate_S3_iii_pass"] else "FAILED",
             detail=f"max |yaw difference| {iii_max:.2e} rad (need < 1e-6)", source="results/stage3_gate.json"),
        dict(gate="S3 overall (pre-registered)", configuration="conditions (i)-(iii)", status="PASSED" if s3["gate_S3_pass"] else "FAILED",
             detail="(i) failed", source="results/stage3_gate.json"),
        dict(gate="S3' overall (amended, post hoc)", configuration="(i)' and (ii) and (iii)", status="PASSED" if s3a["S3_prime_overall_pass"] else "FAILED",
             detail="reported next to, never instead of, S3", source="results/stage3_gate_amended.json"),
    ]
    with open(RES / "gate_table.csv", "w", newline="", encoding="utf8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
    with open(RES / "gate_table.md", "w", encoding="utf8") as fh:
        fh.write("| Gate | Configuration | Status | Detail |\n|---|---|---|---|\n")
        for r in rows:
            fh.write(f"| {r['gate']} | {r['configuration']} | **{r['status']}** | {r['detail']} |\n")
    print(open(RES / "gate_table.md", encoding="utf8").read())


if __name__ == "__main__":
    main()
