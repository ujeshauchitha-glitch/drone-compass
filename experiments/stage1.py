"""Stage 1: spin-dynamics model, Gate S1 checks, anisotropy tables, robustness.  Writes results/stage1_*."""
from __future__ import annotations

import csv
import json
import sys
import time
from itertools import product
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from model import checks  # noqa: E402
from model.sphere import angles_to_vec, fibonacci_sphere, sphere_extrema  # noqa: E402
from model.spin import GAMMA_E, H16, RECALLED, RPModel, primary_model  # noqa: E402
from scipy.optimize import minimize_scalar  # noqa: E402

RES = ROOT / "results"
RES.mkdir(exist_ok=True)


def aniso_row(model, B0, k, label, n_grid=2000, **extra):
    f = lambda v: model.yield_eig(v, B0, k)
    r = sphere_extrema(f, n_grid)
    row = dict(label=label, B0_uT=B0 * 1000, k=k, phi_min=r["min"], phi_max=r["max"], phi_mean=r["mean"],
               aniso=r["aniso"], aniso_rel=r["aniso"] / r["mean"],
               max_at=np.round(r["max_at"], 3).tolist(), min_at=np.round(r["min_at"], 3).tolist())
    row.update(extra)
    return row


def axial(a_iso, d):
    return np.diag([a_iso - d, a_iso - d, a_iso + 2 * d])


def write_csv(path, rows):
    keys = list(rows[0].keys())
    for r in rows:
        for kk in r:
            if kk not in keys:
                keys.append(kk)
    with open(path, "w", newline="", encoding="utf8") as fh:
        w = csv.DictWriter(fh, fieldnames=keys)
        w.writeheader()
        w.writerows(rows)


def meridian_aniso(model, B0, k, relax):
    """Scan-based anisotropy for the Liouville (relaxation) solver: max/min along two meridians + 1-D polish."""
    th = np.linspace(0, np.pi / 2, 13)
    best = {"max": -1, "min": 2}
    for ph in (0.0, np.pi / 2):
        v = np.array([angles_to_vec(t, ph) for t in th])
        y = model.yield_liouville(v, B0, k, relax=relax)
        for key, sgn in (("max", -1), ("min", 1)):
            i = int(np.argmax(y) if key == "max" else np.argmin(y))
            lo, hi = th[max(i - 1, 0)], th[min(i + 1, len(th) - 1)]
            g = lambda t: sgn * float(model.yield_liouville(angles_to_vec(t, ph)[None], B0, k, relax=relax)[0])
            r = minimize_scalar(g, bounds=(lo, hi), method="bounded", options={"xatol": 1e-4, "maxiter": 25})
            val = sgn * r.fun
            best[key] = max(best["max"], val) if key == "max" else min(best["min"], val)
    return best["max"], best["min"]


def main():
    t0 = time.time()
    out = {}
    B0, k = 0.05, 1e6

    # ------------------------------------------------------------------ Gate S1 checks (primary condition)
    gate = {}
    for nm, ten in (("H16", H16), ("R_recalled", RECALLED)):
        g = {}
        g["a"] = checks.check_a(B0, k, ten)
        g["b"] = checks.check_b(ten, k)
        g["c"] = checks.check_c(B0, k, ten)
        g["d"] = checks.check_d(B0, k, ten)
        g["e"] = checks.check_e(B0, k, ten)
        g["f"] = checks.check_f(B0, k, ten)
        gate[nm] = g
        print(nm, {kk: v["pass"] for kk, v in g.items()}, f"{time.time() - t0:.0f}s", flush=True)

    # sign-flip of B0 (reported, not a gate)
    m = primary_model(H16)
    v = checks.rand_dirs(200, 7)
    gate["B0_sign_flip_max_diff"] = float(np.max(np.abs(m.yield_eig(v, B0, k) - m.yield_eig(v, -B0, k))))

    # ------------------------------------------------------------------ anisotropy over conditions
    rows = []
    for (nm, ten), kk, b in product((("H16", H16), ("R_recalled", RECALLED)), (1e5, 1e6, 1e7), (0.025, 0.05, 0.065)):
        rows.append(aniso_row(primary_model(ten), b, kk, nm, tensors=nm))
    write_csv(RES / "stage1_conditions.csv", rows)
    prim = [r for r in rows if r["tensors"] == "H16" and r["k"] == 1e6 and abs(r["B0_uT"] - 50) < 1e-9][0]
    out["primary"] = prim
    print("primary", prim["phi_min"], prim["phi_max"], prim["aniso"], flush=True)

    # ------------------------------------------------------------------ robustness R1 / R2 / R3 (ASSUMED)
    rob = []
    for a_iso, d in product((0.0, 0.5, 1.5), (0.1, 0.5)):
        mdl = RPModel([(1, H16["N5"]), (1, H16["N10"])], [(0.5, axial(a_iso, d))], name="R1")
        rob.append(aniso_row(mdl, B0, k, "R1", n_grid=600, a_iso_B=a_iso, d_B=d))
        print("R1", a_iso, d, rob[-1]["aniso"], f"{time.time() - t0:.0f}s", flush=True)
    r2 = RPModel([(1, H16["N5"]), (1, H16["N10"]), (0.5, axial(0.8, 0.1))], [(0.5, axial(0.5, 0.5))], name="R2")
    rob.append(aniso_row(r2, B0, k, "R2", n_grid=600))
    print("R2", rob[-1]["aniso"], flush=True)
    r3 = RPModel([(1, H16["N5"]), (1, H16["N10"]), (0.5, axial(0.8, 0.1)), (0.5, axial(0.4, 0.05))],
                 [(0.5, axial(0.5, 0.5)), (0.5, axial(0.3, 0.05))], name="R3")
    rob.append(aniso_row(r3, B0, k, "R3 (EXPLORATORY)", n_grid=400))
    print("R3", rob[-1]["aniso"], f"{time.time() - t0:.0f}s", flush=True)

    # ------------------------------------------------------------------ EXPLORATORY: exchange, dipolar, relaxation
    for J in (0.0, 0.05, -0.05, 0.5, -0.5):
        mdl = RPModel([(1, H16["N5"]), (1, H16["N10"])], J=J)
        rob.append(aniso_row(mdl, B0, k, "EXPLORATORY exchange", J_mT=J))
    for D in (0.1, 0.3):
        mdl = RPModel([(1, H16["N5"]), (1, H16["N10"])], D=D)
        rob.append(aniso_row(mdl, B0, k, "EXPLORATORY dipolar", D_mT=D))
    for J, D in ((0.05, 0.3), (-0.5, 0.3)):
        mdl = RPModel([(1, H16["N5"]), (1, H16["N10"])], J=J, D=D)
        rob.append(aniso_row(mdl, B0, k, "EXPLORATORY exchange+dipolar", J_mT=J, D_mT=D))
    print("exchange/dipolar done", f"{time.time() - t0:.0f}s", flush=True)
    base = primary_model(H16)
    for r_rate in (0.0, 1e5, 1e6, 1e7):
        mx, mn = meridian_aniso(base, B0, k, r_rate)
        rob.append(dict(label="EXPLORATORY relaxation (meridian scan)", B0_uT=50, k=k, relax_rate=r_rate,
                        phi_max=mx, phi_min=mn, aniso=mx - mn, aniso_rel=(mx - mn) / (0.5 * (mx + mn))))
        print("relax", r_rate, mx - mn, f"{time.time() - t0:.0f}s", flush=True)
    write_csv(RES / "stage1_robustness.csv", rob)

    # ------------------------------------------------------------------ yield map of the primary condition (5 deg grid)
    th = np.deg2rad(np.arange(0, 181, 5.0))
    ph = np.deg2rad(np.arange(0, 360, 5.0))
    TT, PP = np.meshgrid(th, ph, indexing="ij")
    vv = np.stack([np.sin(TT) * np.cos(PP), np.sin(TT) * np.sin(PP), np.cos(TT)], -1).reshape(-1, 3)
    Y = base.yield_eig(vv, B0, k).reshape(TT.shape)
    with open(RES / "stage1_yield_map_primary.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["theta_deg", "phi_deg", "phi_S"])
        for i, j in product(range(len(th)), range(len(ph))):
            w.writerow([np.rad2deg(th[i]), np.rad2deg(ph[j]), Y[i, j]])

    # ------------------------------------------------------------------ gate verdict
    h = gate["H16"]
    g_pass = all(h[c]["pass"] for c in "abcdef")
    aniso_pass = prim["aniso"] >= checks.ANISO_MIN
    out["gate"] = gate
    out["gate_S1_pass"] = bool(g_pass and aniso_pass)
    out["aniso_threshold"] = checks.ANISO_MIN
    out["runtime_s"] = time.time() - t0
    json.dump(out, open(RES / "stage1_gate.json", "w"), indent=1, default=lambda o: o.tolist() if hasattr(o, "tolist") else str(o))
    print("GATE S1:", out["gate_S1_pass"], f"runtime {out['runtime_s']:.0f}s")

    # ------------------------------------------------------------------ figures
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(1, 2, figsize=(10, 4))
    im = ax[0].pcolormesh(np.rad2deg(PP), np.rad2deg(TT), Y, shading="auto", cmap="viridis")
    ax[0].invert_yaxis()
    ax[0].set_xlabel("azimuth φ (deg)"); ax[0].set_ylabel("polar angle θ from flavin z (deg)")
    ax[0].set_title("Singlet yield Φ_S(n), 50 µT, k = 1e6/s")
    fig.colorbar(im, ax=ax[0], label="Φ_S")
    ax[1].plot(np.rad2deg(th), Y[:, 0], label="φ = 0 (x–z plane)")
    ax[1].plot(np.rad2deg(th), Y[:, 18], label="φ = 90° (y–z plane)")
    ax[1].set_xlabel("θ (deg)"); ax[1].set_ylabel("Φ_S"); ax[1].legend(); ax[1].grid(alpha=.3)
    ax[1].set_title(f"ΔΦ = {prim['aniso']:.3f} (verified H16 tensors)")
    fig.tight_layout(); fig.savefig(RES / "fig_stage1_yield_map.png", dpi=150); plt.close(fig)

    fig, ax = plt.subplots(figsize=(9, 4))
    sub = [r for r in rows if r["tensors"] == "H16"]
    labels = [f"k={r['k']:.0e}\n{r['B0_uT']:.0f}µT" for r in sub]
    ax.bar(range(len(sub)), [r["aniso"] for r in sub], color="#2a6f97")
    ax.set_xticks(range(len(sub))); ax.set_xticklabels(labels, fontsize=7)
    ax.axhline(checks.ANISO_MIN, color="crimson", ls="--", label="pre-registered threshold 1e-3")
    ax.set_yscale("log"); ax.set_ylabel("ΔΦ = Φmax − Φmin"); ax.legend(); ax.grid(alpha=.3, axis="y")
    ax.set_title("Anisotropy of the minimal two-nitrogen model (H16)")
    fig.tight_layout(); fig.savefig(RES / "fig_stage1_aniso_conditions.png", dpi=150); plt.close(fig)


if __name__ == "__main__":
    main()
