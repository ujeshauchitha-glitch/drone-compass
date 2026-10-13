"""Investigate the disagreement (> 10x) with Ren et al. Table 2 (A1.4): is it the number of independently oriented elements?
CRB only, Ren's phenomenological signal, K random element orientations (seeded).  Nothing is tuned."""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
from scipy.spatial.transform import Rotation

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from model.sphere import fibonacci_hemisphere  # noqa: E402
from sensor.poisson import RenPhi, Sensor  # noqa: E402

RES = ROOT / "results"
SEED = 20261013
N_REF = 1e8


def main():
    grid = fibonacci_hemisphere(1000)
    out = {}
    cases = {"default_a0.5_b0.00127_c26": (RenPhi(0.5, 0.00127, 26.0), 2.29e8, 5.26e7),
             "optimised_a0.9_b0.00508_c10": (RenPhi(0.9, 0.00508, 10.0), None, 3.94e5)}
    for nm, (rp, ren_o1, ren_o2) in cases.items():
        rows = []
        for K in (3, 12, 50, 200, 1000):
            rots = Rotation.random(K, random_state=np.random.default_rng(np.random.SeedSequence(SEED, spawn_key=(9, K)))).as_matrix()
            sen = Sensor(rp, list(rots), eps=1.0)
            _, smed, sing = sen.crb(grid, N_REF)
            sd = np.degrees(smed[np.isfinite(smed)])
            med = float(np.median(sd))
            rows.append(dict(K=K, sigma_med_deg_at_1e8_per_element=med, N_total_for_5deg_median=K * N_REF * (med / 5.0) ** 2,
                             worst_direction_deg=float(sd.max()), singular=int(sing.sum())))
        out[nm] = dict(rows=rows, ren_reported_total_photons={"O1": ren_o1, "O2": ren_o2})
        print(nm, [(r["K"], f"{r['N_total_for_5deg_median']:.3g}") for r in rows], "Ren O1/O2:", ren_o1, ren_o2, flush=True)
    json.dump(out, open(RES / "stage2_ren_investigation.json", "w"), indent=1)


if __name__ == "__main__":
    main()
