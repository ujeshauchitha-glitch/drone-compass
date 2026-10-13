"""Tornado analysis (Amendment 3 A3.4): rank ASSUMED parameters by their swing in the best-case radical-pair benefit."""
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
BLUE, ORANGE, GREY = "#0072B2", "#D55E00", "#6b6b6b"


def main():
    rows = list(csv.DictReader(open(RES / "stage3_tornado.csv")))
    out = {}
    fig, axs = plt.subplots(1, 2, figsize=(11, 4.6), sharex=True)
    for ax, filt, title in zip(axs, ("uncal", "cal"), ("pre-registered (uncalibrated) filters", "online-calibrated filters")):
        base = [r for r in rows if r["filter"] == filt and r["parameter"] == "base"][0]
        b = float(base["metric_best_min60_600"])
        params = []
        for name in dict.fromkeys(r["parameter"] for r in rows):
            if name == "base":
                continue
            vals = [(r["value"], float(r["metric_best_min60_600"]), float(r["metric_best_reduction_600"])) for r in rows if r["filter"] == filt and r["parameter"] == name]
            lo, hi = min(v[1] for v in vals), max(v[1] for v in vals)
            params.append(dict(parameter=name, values=[(v[0], round(v[1], 4), round(v[2], 4)) for v in vals], swing=max(hi, b) - min(lo, b), min_metric=min(lo, b), max_metric=max(hi, b),
                               reaches_20pct=bool(hi >= 0.2), best_reduction_600=max(v[2] for v in vals)))
        params.sort(key=lambda p: -p["swing"])
        out[filt] = dict(base_metric=b, base_best_reduction_600=float(base["metric_best_reduction_600"]), ranking=params)
        y = np.arange(len(params))[::-1]
        for yy, p in zip(y, params):
            ax.barh(yy, p["max_metric"] - p["min_metric"], left=p["min_metric"], color=BLUE if filt == "cal" else ORANGE, alpha=.85)
            ax.text(p["max_metric"] + 0.01, yy, f"{p['swing']*100:.0f} pts", va="center", fontsize=7)
        ax.axvline(b, color=GREY, lw=1); ax.axvline(0.2, color="crimson", ls="--", lw=1, label="20 % criterion")
        ax.set_yticks(y); ax.set_yticklabels([p["parameter"] for p in params], fontsize=7)
        ax.set_title(title, fontsize=9); ax.set_xlabel("best-case benefit m = max over σ_rp of min(reduction₆₀, reduction₆₀₀)")
        ax.grid(alpha=.3, axis="x")
    axs[0].legend(fontsize=7, loc="lower right")
    fig.tight_layout(); fig.savefig(RES / "fig_stage4_tornado.png", dpi=150); plt.close(fig)
    json.dump(out, open(RES / "stage3_tornado_summary.json", "w"), indent=1)
    for filt in ("uncal", "cal"):
        print(filt, "base m = %.3f" % out[filt]["base_metric"], [(p["parameter"], round(p["swing"], 3), p["reaches_20pct"]) for p in out[filt]["ranking"][:6]])


if __name__ == "__main__":
    main()
