"""Amendment 1, A1.2: select the element set by the pre-registered rule (design seeds / design grid only)."""
from __future__ import annotations

import json
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
from scipy.optimize import differential_evolution, minimize

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from model.sphere import fibonacci_hemisphere  # noqa: E402
from model.spin import primary_model  # noqa: E402
from model.surrogate import PolySurrogate  # noqa: E402
from sensor.design import N_DESIGN, X_DEG, build_rots, min_kl, objective, refined_jref, s_matrix, separation_mask  # noqa: E402
from sensor.poisson import ELEMENT_SETS, rotvec  # noqa: E402

RES = ROOT / "results"
SEED = 20261013
B0, K_RATE = 0.05, 1e6
DESIGN_GRID_N, CONFIRM_GRID_N = 600, 3000


def _phi():
    m = primary_model()
    return PolySurrogate.fit(lambda v: m.yield_eig(v, B0, K_RATE), 8)


def run_one(task):
    K, r = task
    phi = _phi()
    phibar = phi.mean
    grid = fibonacci_hemisphere(DESIGN_GRID_N)
    mask = separation_mask(grid)
    n_par = 3 * (K - 1)
    rng = np.random.default_rng(np.random.SeedSequence(SEED, spawn_key=(3, r)))
    f = lambda p: -objective(p, K, phi, phibar, grid, mask)
    t0 = time.time()
    de = differential_evolution(f, [(-np.pi, np.pi)] * n_par, popsize=15, maxiter=60, polish=False, rng=rng, tol=0)
    nm = minimize(f, de.x, method="Nelder-Mead", options={"xatol": 1e-6, "fatol": 1e-9, "maxiter": 600})
    best = nm.x if nm.fun < de.fun else de.x
    Jd = -min(nm.fun, de.fun)
    # confirmation on an independent, finer grid
    cg = fibonacci_hemisphere(CONFIRM_GRID_N)
    rots = build_rots(best, K)
    Jc = min_kl(s_matrix(phi, phibar, rots, cg), separation_mask(cg))
    Jr_d = refined_jref(phi, phibar, rots, grid, mask)
    Jr_c = refined_jref(phi, phibar, rots, cg)
    return dict(K=K, restart=r, params=best.tolist(), J_design=Jd, J_confirm=Jc, Jref_design=Jr_d, Jref_confirm=Jr_c,
                de_J=-de.fun, seconds=time.time() - t0)


def main():
    t0 = time.time()
    tasks = [(K, r) for K in (3, 2, 4) for r in range(3)]
    with ProcessPoolExecutor(max_workers=2) as ex:
        results = list(ex.map(run_one, tasks))
    out = {"X_deg": X_DEG, "N_design": N_DESIGN, "design_grid": DESIGN_GRID_N, "confirm_grid": CONFIRM_GRID_N, "runs": results}
    sel = {}
    for K in (3, 2, 4):
        rs = [r for r in results if r["K"] == K]
        b = max(rs, key=lambda r: r["Jref_design"])  # Amendment 1.1: selection by the refined objective
        sel[str(K)] = dict(restart=b["restart"], params=b["params"], J_design=b["J_design"], J_confirm=b["J_confirm"],
                           Jref_design=b["Jref_design"], Jref_confirm=b["Jref_confirm"],
                           J_all_restarts=[r["J_design"] for r in rs], Jref_all_restarts=[r["Jref_design"] for r in rs])
    out["selected"] = sel
    # comparison scores for the labelled reference sets (design grid and confirmation grid)
    phi = _phi()
    phibar = phi.mean
    grid, cg = fibonacci_hemisphere(DESIGN_GRID_N), fibonacci_hemisphere(CONFIRM_GRID_N)
    gen = [np.eye(3), rotvec(np.array([0.3, 0.5, 0.2])), rotvec(np.array([-0.4, 0.2, 0.6]))]
    refs = {}
    for name, rots in (("K3orth", ELEMENT_SETS["K3orth"]), ("K3generic_EXPLORATORY", gen),
                       ("K3selected", build_rots(sel["3"]["params"], 3))):
        refs[name] = dict(J_design=min_kl(s_matrix(phi, phibar, rots, grid), separation_mask(grid)),
                          J_confirm=min_kl(s_matrix(phi, phibar, rots, cg), separation_mask(cg)),
                          Jref_design=refined_jref(phi, phibar, rots, grid),
                          Jref_confirm=refined_jref(phi, phibar, rots, cg))
    out["reference_scores"] = refs
    out["selected_rotvecs_K3"] = sel["3"]["params"]
    out["runtime_s"] = time.time() - t0
    json.dump(out, open(RES / "stage2_design.json", "w"), indent=1)
    for K in (3, 2, 4):
        print("K", K, "grid J per restart", np.round(sel[str(K)]["J_all_restarts"], 3), "Jref per restart",
              np.round(sel[str(K)]["Jref_all_restarts"], 3), "selected", sel[str(K)]["restart"],
              "Jref_confirm", sel[str(K)]["Jref_confirm"])
    print("reference (grid J, Jref design, Jref confirm)",
          {k: (round(v["J_design"], 3), round(v["Jref_design"], 4), round(v["Jref_confirm"], 4)) for k, v in refs.items()})
    print("runtime", out["runtime_s"])


if __name__ == "__main__":
    main()
