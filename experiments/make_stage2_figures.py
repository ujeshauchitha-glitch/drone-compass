"""Stage 2 figures from results/ files."""
from __future__ import annotations

import csv
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
RES = ROOT / "results"
BLUE, ORANGE, GREEN, GREY = "#0072B2", "#D55E00", "#009E73", "#6b6b6b"  # Okabe-Ito


def main():
    s2 = json.load(open(RES / "stage2_gate.json"))
    sp = json.load(open(RES / "stage2_sprime.json"))["gate_S2prime"]
    F = list(csv.DictReader(open(RES / "stage2_failure_rates.csv")))
    R = list(csv.DictReader(open(RES / "stage2_requirements.csv")))

    # ---- error vs N
    fig, ax = plt.subplots(figsize=(6.2, 4.2))
    N = [r["N"] for r in sp["rows"]]
    ax.loglog(N, [r["crb_median_of_sigmed"] for r in sp["rows"]], "-", color=GREY, label="Cramér–Rao prediction (selected set)")
    ax.loglog(N, [r["mc_median"] for r in sp["rows"]], "o", color=BLUE, label="Monte-Carlo median, selected set (Gate S2′)")
    ko = s2["K3orth"]["rows"]
    ax.loglog([r["N"] for r in ko], [r["mc_median"] for r in ko], "s", color=ORANGE, label="Monte-Carlo median, K3orth (Gate S2 FAILED)")
    ax.set_xlabel("detected photons per element N"); ax.set_ylabel("axis error (deg)")
    ax.grid(alpha=.3, which="both"); ax.legend(fontsize=8); ax.set_title("Median axis error vs photons (ε = 1, Poisson only)")
    fig.tight_layout(); fig.savefig(RES / "fig_stage2_error_vs_N.png", dpi=150); plt.close(fig)

    # ---- failure rate vs N
    fig, axs = plt.subplots(1, 2, figsize=(9.5, 3.8), sharey=True)
    cols = {"K3orth": ORANGE, "K3generic_EXPLORATORY": GREEN, "K3selected": BLUE}
    for ax, eps in zip(axs, ("1.0", "0.1")):
        for name, c in cols.items():
            rows = [r for r in F if r["set"] == name and r["eps"] == eps]
            ax.semilogx([float(r["N"]) for r in rows], [float(r["fail_rate_gt_10deg"]) for r in rows], "-o", ms=3, color=c, label=name)
        ax.set_title(f"ε = {eps}"); ax.set_xlabel("detected photons per element N"); ax.grid(alpha=.3)
    axs[0].set_ylabel("fraction of trials with error > 10°"); axs[0].legend(fontsize=7)
    fig.tight_layout(); fig.savefig(RES / "fig_stage2_failure_rates.png", dpi=150); plt.close(fig)

    # ---- required photons vs contrast (selected set, baseline realism)
    fig, ax = plt.subplots(figsize=(6.2, 4.2))
    rows = [r for r in R if r["set"] == "K3selected" and r["realism"] == "baseline"]
    eps = np.array([float(r["eps"]) for r in rows]); o = np.argsort(eps)
    ax.loglog(eps[o], np.array([float(r["N_5deg_median_dir"]) for r in rows])[o], "o-", color=BLUE, label="5°, median direction")
    ax.loglog(eps[o], np.array([float(r["N_5deg_worst_dir"]) for r in rows])[o], "s-", color=ORANGE, label="5°, worst direction")
    ax.loglog(eps[o], np.array([float(r["N_1deg_median_dir"]) for r in rows])[o], "o--", color=BLUE, alpha=.6, label="1°, median direction")
    ax.loglog(eps[o], np.array([float(r["N_1deg_worst_dir"]) for r in rows])[o], "s--", color=ORANGE, alpha=.6, label="1°, worst direction")
    ax.set_xlabel("relative contrast factor ε (1 = minimal model; DERIVED values for richer/relaxed models)")
    ax.set_ylabel("photons per element (CRB-derived)"); ax.grid(alpha=.3, which="both"); ax.legend(fontsize=8)
    ax.set_title("Photon requirement, selected K=3 set")
    fig.tight_layout(); fig.savefig(RES / "fig_stage2_requirements.png", dpi=150); plt.close(fig)


if __name__ == "__main__":
    main()
