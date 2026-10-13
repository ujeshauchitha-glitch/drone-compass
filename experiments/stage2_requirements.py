"""PLAN.md 4.8 + Amendment 1 A1.6: requirements table with ambiguity failure rate and relaxation margin."""
from __future__ import annotations

import csv
import json
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from experiments.stage2_sprime import FAIL_DEG, mc, sets  # noqa: E402
from model.sphere import fibonacci_hemisphere  # noqa: E402
from model.spin import primary_model  # noqa: E402
from model.surrogate import PolySurrogate  # noqa: E402
from sensor.poisson import Sensor, hemisphere_directions  # noqa: E402

RES = ROOT / "results"
SEED = 20261013
B0, K_RATE = 0.05, 1e6
ETA, T_UPDATE = 2 / np.pi, 0.1  # ASSUMED lock-in efficiency; 10 Hz update
N_REF = 1e8
TARGETS = (5.0, 1.0)
ANISO_THRESHOLD = 1e-3


def read_csv(p):
    with open(p, newline="", encoding="utf8") as fh:
        return list(csv.DictReader(fh))


def derived_eps():
    cond = read_csv(RES / "stage1_conditions.csv")
    prim = [r for r in cond if r["tensors"] == "H16" and float(r["k"]) == 1e6 and abs(float(r["B0_uT"]) - 50) < 1e-6][0]
    rel0 = float(prim["aniso"]) / float(prim["phi_mean"])
    rob = read_csv(RES / "stage1_robustness.csv")
    out = {}
    for r in rob:
        lab = r["label"]
        if lab == "R2":
            out["eps_R2"] = float(r["aniso_rel"]) / rel0
        if lab.startswith("EXPLORATORY relaxation") and float(r["relax_rate"]) in (1e6, 1e7):
            out[f"eps_relax_{float(r['relax_rate']):.0e}"] = float(r["aniso_rel"]) / rel0
    return out, rel0


def stat_sigma(sen, grid, N):
    _, smed, _ = sen.crb(grid, N)
    deg = np.degrees(smed)
    return float(np.median(deg)), float(np.max(deg))


def n_required(sen, grid, target, which):
    """Photons per element for the median-over-directions ('med') or worst-direction ('worst') sigma_med to reach target."""
    s_ref = stat_sigma(sen, grid, N_REF)[0 if which == "med" else 1]
    if not np.isfinite(s_ref):
        return np.inf
    if sen.sigma_g == 0:
        return N_REF * (s_ref / target) ** 2  # J proportional to N exactly when sigma_g = 0
    lo, hi = 1e2, 1e20
    for _ in range(80):
        mid = np.sqrt(lo * hi)
        s = stat_sigma(sen, grid, mid)[0 if which == "med" else 1]
        lo, hi = (mid, hi) if s > target else (lo, mid)
    return float(np.sqrt(lo * hi))


def main():
    t0 = time.time()
    design = json.load(open(RES / "stage2_design.json"))
    m = primary_model()
    phi = PolySurrogate.fit(lambda v: m.yield_eig(v, B0, K_RATE), 8)
    S = sets(design)
    grid = fibonacci_hemisphere(2000)
    eps_d, rel0 = derived_eps()
    realism = [("baseline", dict()), ("background beta=1", dict(beta=1.0)), ("background beta=10", dict(beta=10.0)),
               ("electronics rho=1", dict(rho=1.0)), ("electronics rho=10", dict(rho=10.0)),
               ("gain drift sigma_g=1e-3", dict(sigma_g=1e-3)), ("gain drift sigma_g=1e-2", dict(sigma_g=1e-2)),
               ("gain drift sigma_g=1e-1", dict(sigma_g=1e-1))]
    eps_list = [(1.0, "ASSUMED"), (0.1, "ASSUMED")] + [(v, "DERIVED from Stage 1 anisotropy ratio (approximate)") for v in eps_d.values()]
    eps_names = ["1", "0.1"] + list(eps_d.keys())
    rows = []
    for name in ("K3selected", "K3orth"):
        for (eps, tag), en in zip(eps_list, eps_names):
            for rname, kw in realism:
                if en not in ("1", "0.1") and rname != "baseline":
                    continue
                sen = Sensor(phi, S[name], eps=eps, **kw)
                row = dict(set=name, eps_name=en, eps=eps, eps_tag=tag, realism=rname)
                med1, worst1 = stat_sigma(sen, grid, N_REF)
                row["sigma_med_deg_at_1e8_median_dir"], row["sigma_med_deg_at_1e8_worst_dir"] = med1, worst1
                for tgt in TARGETS:
                    nm, nw = n_required(sen, grid, tgt, "med"), n_required(sen, grid, tgt, "worst")
                    row[f"N_{tgt:g}deg_median_dir"], row[f"N_{tgt:g}deg_worst_dir"] = nm, nw
                    row[f"rate_{tgt:g}deg_median_dir_per_s"] = nm / (ETA * T_UPDATE)
                    rho = kw.get("rho", 0.0)
                    row[f"nu_e_{tgt:g}deg_median_dir"] = rho * np.sqrt(nm / T_UPDATE) if rho else 0.0
                rows.append(row)
        print("requirements:", name, f"{time.time() - t0:.0f}s", flush=True)
    keys = list(rows[0].keys())
    with open(RES / "stage2_requirements.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=keys); w.writeheader(); w.writerows(rows)

    # ---------------- Monte-Carlo verification at the CRB-derived N (baseline realism), fresh seeds (prefix 7)
    dirs = hemisphere_directions(300, np.random.default_rng(np.random.SeedSequence(SEED, spawn_key=(7, 0))))
    ver = []
    idx = 0
    for name in ("K3selected", "K3orth"):
        for eps, en in ((1.0, "1"), (0.1, "0.1")):
            sen = Sensor(phi, S[name], eps=eps)
            sen.make_grid(4000)
            for tgt in TARGETS:
                r0 = [r for r in rows if r["set"] == name and r["eps_name"] == en and r["realism"] == "baseline"][0]
                N = r0[f"N_{tgt:g}deg_median_dir"]
                rng = np.random.default_rng(np.random.SeedSequence(SEED, spawn_key=(7, 1, idx)))
                idx += 1
                err = mc(sen, dirs, N, rng)
                ver.append(dict(set=name, eps=eps, target_deg=tgt, N=N, mc_median_deg=float(np.median(err)),
                                mc_p90_deg=float(np.quantile(err, .9)), ambiguity_failure_rate_gt_10deg=float((err > FAIL_DEG).mean())))
        print("verification:", name, f"{time.time() - t0:.0f}s", flush=True)
    with open(RES / "stage2_requirements_mc.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(ver[0].keys())); w.writeheader(); w.writerows(ver)

    # ---------------- relaxation / richness margin against the pre-registered 1e-3 anisotropy threshold
    marg = []
    cond = read_csv(RES / "stage1_conditions.csv")
    prim = [r for r in cond if r["tensors"] == "H16" and float(r["k"]) == 1e6 and abs(float(r["B0_uT"]) - 50) < 1e-6][0]
    marg.append(dict(model="minimal H16 (primary)", aniso=float(prim["aniso"]), source="stage1_conditions.csv"))
    for r in read_csv(RES / "stage1_robustness.csv"):
        lab = r["label"]
        if lab in ("R1", "R2") or lab.startswith("R3") or lab.startswith("EXPLORATORY"):
            tag = lab
            for k_ in ("a_iso_B", "d_B", "J_mT", "D_mT", "relax_rate"):
                if r.get(k_, "") != "":
                    tag += f" {k_}={r[k_]}"
            marg.append(dict(model=tag, aniso=float(r["aniso"]), source="stage1_robustness.csv"))
    for r in marg:
        r["margin_vs_1e-3"] = r["aniso"] / ANISO_THRESHOLD
        r["above_threshold"] = r["aniso"] >= ANISO_THRESHOLD
    with open(RES / "stage2_margins.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(marg[0].keys())); w.writeheader(); w.writerows(marg)
    json.dump(dict(eps_derived=eps_d, rel_aniso_primary=rel0, eta=ETA, T_update_s=T_UPDATE, runtime_s=time.time() - t0),
              open(RES / "stage2_requirements_meta.json", "w"), indent=1)
    print("done", f"{time.time() - t0:.0f}s")


if __name__ == "__main__":
    main()
