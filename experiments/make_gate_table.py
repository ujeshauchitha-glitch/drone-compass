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
    with open(RES / "gate_table.csv", "w", newline="", encoding="utf8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
    with open(RES / "gate_table.md", "w", encoding="utf8") as fh:
        fh.write("| Gate | Configuration | Status | Detail |\n|---|---|---|---|\n")
        for r in rows:
            fh.write(f"| {r['gate']} | {r['configuration']} | **{r['status']}** | {r['detail']} |\n")
    print(open(RES / "gate_table.md", encoding="utf8").read())


if __name__ == "__main__":
    main()
