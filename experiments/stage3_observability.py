"""Amendment 2 A2.4: numerical local observability of the 6-state error model along the true attitude."""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from nav import attitude as A  # noqa: E402

RES = ROOT / "results"
SIG_ACC = float(np.hypot(0.05, 0.3))
SIG_MAG, SIG_RP_DEG = 0.5, 1.0
S_B = 100.0 / 3600.0 * A.DEG  # state scaling for the bias block (100 deg/h) so singular values are comparable


def gramian(arm, dip, k0, k1, scen):
    be = A.earth_field(dip); bn = np.linalg.norm(be); un = be / bn
    gz = np.array([0, 0, A.G])
    sig_rp = SIG_RP_DEG * A.DEG
    W = np.zeros((6, 6))
    Phi = np.eye(6)
    for k in range(k0, k1):
        R = scen.R[k + 1]
        F = np.eye(6)
        F[:3, 3:] = -scen.R[k] * 0 - scen.R[k] * A.DT   # error propagation with R_hat ~ R_true
        Phi = F @ Phi
        rows = [(np.hstack([R.T @ A.skew(gz), np.zeros((3, 3))]), SIG_ACC**2)]
        if arm in ("B", "D"):
            rows.append((np.hstack([R.T @ A.skew(be), np.zeros((3, 3))]), SIG_MAG**2))
        if arm in ("C", "D"):
            rows.append((np.hstack([R.T @ A.skew(un), np.zeros((3, 3))]), sig_rp**2))
        for H, rv in rows:
            HP = H @ Phi
            W += HP.T @ HP / rv
    return W


def analyse(W):
    D = np.diag([1, 1, 1, S_B, S_B, S_B])
    Ws = D @ W @ D
    ev, V = np.linalg.eigh(Ws)
    ev = np.clip(ev, 0, None)
    rank = int(np.sum(ev > 1e-12 * ev.max()))
    null = [np.round(V[:, i], 3).tolist() for i in range(6) if ev[i] <= 1e-12 * ev.max()]
    # initial-state yaw uncertainty implied by the information (unbounded if rank deficient in the yaw direction)
    yaw_dir = np.zeros(6); yaw_dir[2] = 1.0
    unbounded = any(abs(v[2]) > 1e-3 for v in null) if null else False
    if unbounded:
        yaw_std = float("inf")
    else:
        yaw_std = float(np.degrees(np.sqrt(yaw_dir @ np.linalg.solve(Ws + 1e-300 * np.eye(6), yaw_dir))))
    return dict(singular_values=[float(x) for x in ev[::-1]], rank=rank, null_directions_scaled=null, yaw_std_deg_from_information=yaw_std)


def main():
    sc = A.Scenario()
    out = []
    windows = {"hover_0_120s": (0, 1200), "yaw_turn_130_148s": (1300, 1480)}
    for dip in (15.0, 45.0, 66.0):
        for wname, (k0, k1) in windows.items():
            for arm in "ABCD":
                res = analyse(gramian(arm, dip, k0, k1, sc))
                res.update(dip=dip, window=wname, arm=arm)
                out.append(res)
                print(f"dip {dip:>4} {wname:18s} arm {arm}: rank {res['rank']}  yaw std from information {res['yaw_std_deg_from_information']:.4g} deg")
    json.dump(out, open(RES / "stage3_observability.json", "w"), indent=1)


if __name__ == "__main__":
    main()
