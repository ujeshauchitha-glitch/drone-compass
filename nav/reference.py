"""Independent second code path for Gate S3(iii): single run, plain Python loop, rotation MATRICES (no quaternions),
scipy expm/Rotation for the rotation maps, and NUMERICAL (central-difference) Jacobians for both the error propagation
and the measurement models.  Same models as nav/attitude.py (PLAN.md Amendment 2), written independently."""
from __future__ import annotations

import numpy as np
from scipy.linalg import expm
from scipy.spatial.transform import Rotation as Rot

from nav.attitude import DEG, DT, G, N_STEP, earth_field, tangent_basis_rows

H_FD = 1e-6


def sk(v):
    return np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0.0]])


def E(v):
    return expm(sk(v))


def logv(R):
    return Rot.from_matrix(R).as_rotvec()


def numjac(f, x0, h=H_FD):
    cols = []
    for i in range(len(x0)):
        d = np.zeros(len(x0)); d[i] = h
        cols.append((f(x0 + d) - f(x0 - d)) / (2 * h))
    return np.stack(cols, axis=1)


def ref_run(arm, cfg, scen, rnd, r, n_steps=N_STEP):
    use_mag, use_rp = arm in ("B", "D"), arm in ("C", "D")
    be = earth_field(cfg["dip"]); bn = np.linalg.norm(be); un = be / bn
    sd = cfg["sigma_d"]
    bias_rps = cfg["bias_dph"] / 3600.0 * DEG
    arw = cfg["arw"] / 60.0 * DEG
    sig_bw = 0.1 * bias_rps / np.sqrt(600.0)
    sig_a = cfg["acc_sigma"]
    sig_m2 = cfg["sigma_mag"] ** 2 + sd**2 / 3.0
    sig_rp = cfg["sigma_rp_deg"] * DEG
    sig_r2 = sig_rp**2 + sd**2 / (3.0 * bn**2)
    Smat = rnd.soft_iron()[r] if cfg["soft_iron"] else np.eye(3)
    d_hi, d_mc = sd / np.sqrt(6.0) * rnd.dhi[r], sd / np.sqrt(6.0) * rnd.dmc[r]
    dth0 = rnd.init[r] * np.array([1.0, 1.0, cfg["init_yaw_deg"]]) * DEG
    Rh = E(-dth0) @ scen.R[0]
    bh = np.zeros(3)
    b_true = bias_rps * rnd.b0[r]
    P = np.zeros((6, 6))
    P[:3, :3] = np.diag(np.array([1.0, 1.0, cfg["init_yaw_deg"]]) * DEG) ** 2
    P[3:, 3:] = np.eye(3) * bias_rps**2
    gz = np.array([0.0, 0.0, G])
    yaw = np.empty(n_steps + 1)
    yaw[0] = logv(scen.R[0] @ Rh.T)[2]

    def update(res, hfun, Rm, Rh, bh, P):
        Hth = numjac(hfun, np.zeros(3))
        H = np.hstack([Hth, np.zeros((3, 3))])
        S = H @ P @ H.T + Rm * np.eye(3)
        K = P @ H.T @ np.linalg.inv(S)
        dx = K @ res
        IKH = np.eye(6) - K @ H
        P = IKH @ P @ IKH.T + Rm * K @ K.T
        Rh = E(dx[:3]) @ Rh
        return Rh, bh + dx[3:], P

    for k in range(n_steps):
        dtm = scen.dtheta[k] + b_true * DT + arw * np.sqrt(DT) * rnd.gyro[r, k]
        u = dtm - bh * DT
        Rh_old = Rh

        def g(x):
            Rerr = E(x[:3]) @ Rh_old @ E(u - x[3:] * DT) @ E(-u) @ Rh_old.T
            return logv(Rerr)
        Jx = numjac(g, np.zeros(6))
        F = np.eye(6)
        F[:3, :] = Jx
        P = F @ P @ F.T
        P[:3, :3] += arw**2 * DT * np.eye(3)
        P[3:, 3:] += sig_bw**2 * DT * np.eye(3)
        Rh = Rh_old @ E(u)
        b_true = b_true + sig_bw * np.sqrt(DT) * rnd.bw[r, k]
        Rt = scen.R[k + 1]
        # accelerometer
        za = scen.fb[k + 1] + sig_a * rnd.acc[r, k]
        Rh, bh, P = update(za - Rh.T @ gz, lambda d, Rh=Rh: (E(d) @ Rh).T @ gz, sig_a**2, Rh, bh, P)
        sv = Smat @ (Rt.T @ be)
        dist = d_hi + scen.c[k + 1] * d_mc
        b_mems = sv + dist
        b_body = sv + cfg["kappa"] * dist
        if use_mag:
            zm = b_mems + cfg["sigma_mag"] * rnd.mag[r, k]
            Rh, bh, P = update(zm - Rh.T @ be, lambda d, Rh=Rh: (E(d) @ Rh).T @ be, sig_m2, Rh, bh, P)
        if use_rp:
            ut = b_body / np.linalg.norm(b_body)
            e1, e2 = tangent_basis_rows(ut[None])
            um = ut + sig_rp * (rnd.rp[r, k, 0] * e1[0] + rnd.rp[r, k, 1] * e2[0])
            um = rnd.sgn[r, k] * um / np.linalg.norm(um)
            up = Rh.T @ un
            s = 1.0 if um @ up >= 0 else -1.0
            Rh, bh, P = update(um - s * up, lambda d, Rh=Rh, s=s: s * ((E(d) @ Rh).T @ un), sig_r2, Rh, bh, P)
        yaw[k + 1] = logv(Rt @ Rh.T)[2]
    return yaw
