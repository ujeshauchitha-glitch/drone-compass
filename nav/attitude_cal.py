"""Online-calibration MEKF (PLAN.md Amendment 3, A3.2).  18 states: dtheta(3), db(3), dh(3, uT), de(6), dd(3, uT).
Same truth generator and random inputs as nav/attitude.py.  Applies equally to the vector sensor (B), the axis-only
sensor (C) and D (shared calibration states)."""
from __future__ import annotations

import numpy as np

from nav.attitude import (DEG, DT, G, R_to_q, R_to_rotvec, earth_field, left_jacobian, q_from_rotvec, q_mul, q_to_R,
                          rotvec_to_R, skew, tangent_basis_rows)

NS = 18
IT, IB, IH, IE, ID = slice(0, 3), slice(3, 6), slice(6, 9), slice(9, 15), slice(15, 18)


def e_matrix(e):
    """(...,6) -> symmetric (...,3,3) with order (xx, yy, zz, xy, xz, yz)."""
    out = np.zeros(e.shape[:-1] + (3, 3))
    out[..., 0, 0], out[..., 1, 1], out[..., 2, 2] = e[..., 0], e[..., 1], e[..., 2]
    out[..., 0, 1] = out[..., 1, 0] = e[..., 3]
    out[..., 0, 2] = out[..., 2, 0] = e[..., 4]
    out[..., 1, 2] = out[..., 2, 1] = e[..., 5]
    return out


def e_jacobian(v):
    """d(E v)/d e : (...,3) -> (...,3,6)."""
    z = np.zeros_like(v[..., 0])
    v1, v2, v3 = v[..., 0], v[..., 1], v[..., 2]
    cols = [np.stack([v1, z, z], -1), np.stack([z, v2, z], -1), np.stack([z, z, v3], -1),
            np.stack([v2, v1, z], -1), np.stack([v3, z, v1], -1), np.stack([z, v3, v2], -1)]
    return np.stack(cols, -1)


def field_model(Rh, be, ehat, hhat, dhat, c):
    """m_hat = (I+E) R^T b_e + h + c d ; returns (m_hat, v = R^T b_e)."""
    v = np.einsum("mji,j->mi", Rh, be)
    m = v + np.einsum("mij,mj->mi", e_matrix(ehat), v) + hhat + c * dhat
    return m, v


def h_field(Rh, be, ehat, v, c):
    """3x18 Jacobian of m_hat w.r.t. the error state (nav-frame left attitude error)."""
    M = Rh.shape[0]
    H = np.zeros((M, 3, NS))
    IpE = np.eye(3) + e_matrix(ehat)
    H[:, :, IT] = IpE @ np.swapaxes(Rh, 1, 2) @ skew(be)
    H[:, :, IH] = np.eye(3)
    H[:, :, IE] = e_jacobian(v)
    H[:, :, ID] = c * np.eye(3)
    return H


def run_filter_cal(arm, cfg, scen, rnd, nees_steps=(), prior_scale=1.0, force_floor=None):
    """arm in {'B','C','D'} with online calibration.  Returns dict(yaw_err (M,T+1), nees {step: (M,)} for the 6-dof block)."""
    M, T = rnd.M, rnd.T
    use_mag, use_rp = arm in ("B", "D"), arm in ("C", "D")
    be = earth_field(cfg["dip"]); bnorm = np.linalg.norm(be); un = be / bnorm
    sd = cfg["sigma_d"]
    bias_rps = cfg["bias_dph"] / 3600.0 * DEG
    arw = cfg["arw"] / 60.0 * DEG
    sig_bw = 0.1 * bias_rps / np.sqrt(600.0)
    q_tt, q_bb = arw**2 * DT, sig_bw**2 * DT
    sig_a2 = cfg["acc_sigma"] ** 2
    sig_m2 = cfg["sigma_mag"] ** 2
    sig_rp = cfg["sigma_rp_deg"] * DEG
    sig_r2 = sig_rp**2
    S = rnd.soft_iron(cfg.get("soft_iron_scale", 0.01)) if cfg["soft_iron"] else np.broadcast_to(np.eye(3), (M, 3, 3)).copy()
    d_hi = sd / np.sqrt(6.0) * rnd.dhi
    d_mc = sd / np.sqrt(6.0) * rnd.dmc
    kappa = cfg["kappa"]
    dth0 = rnd.init * np.array([1.0, 1.0, cfg["init_yaw_deg"]]) * DEG
    R0 = scen.R[0]
    q = R_to_q(rotvec_to_R(-dth0) @ R0)
    bhat = np.zeros((M, 3)); hhat = np.zeros((M, 3)); dhat = np.zeros((M, 3)); ehat = np.zeros((M, 6))
    b_true = bias_rps * rnd.b0
    floor = 0.05 if force_floor is None else force_floor
    cal_std = max(sd / np.sqrt(6.0), floor) * prior_scale
    P = np.zeros((M, NS, NS))
    P[:, :3, :3] = np.diag(np.array([1.0, 1.0, cfg["init_yaw_deg"]]) * DEG) ** 2
    P[:, IB, IB] = np.eye(3) * bias_rps**2
    P[:, IH, IH] = np.eye(3) * cal_std**2
    P[:, ID, ID] = np.eye(3) * cal_std**2
    P[:, IE, IE] = np.eye(6) * (cfg.get("soft_iron_scale", 0.01) * prior_scale) ** 2
    INS = np.eye(NS)
    gz = np.array([0.0, 0.0, G])
    skew_g = skew(gz)
    yaw_err = np.empty((M, T + 1))
    nees = {}

    def err_state(Rt, Rh):
        return R_to_rotvec(Rt @ np.swapaxes(Rh, 1, 2))

    Rh = q_to_R(q)
    yaw_err[:, 0] = err_state(R0, Rh)[:, 2]

    def update(r, H, Rm, q, bhat, hhat, ehat, dhat, P):
        PHt = P @ np.swapaxes(H, 1, 2)
        Sm = H @ PHt + Rm * np.eye(3)
        K = np.swapaxes(np.linalg.solve(Sm, np.swapaxes(PHt, 1, 2)), 1, 2)
        dx = np.einsum("mij,mj->mi", K, r)
        IKH = INS - K @ H
        P = IKH @ P @ np.swapaxes(IKH, 1, 2) + Rm * (K @ np.swapaxes(K, 1, 2))
        q = q_mul(q_from_rotvec(dx[:, IT]), q)
        q = q / np.linalg.norm(q, axis=1, keepdims=True)
        return q, bhat + dx[:, IB], hhat + dx[:, IH], ehat + dx[:, IE], dhat + dx[:, ID], P

    for k in range(T):
        dtm = scen.dtheta[k] + b_true * DT + arw * np.sqrt(DT) * rnd.gyro[:, k]
        u = dtm - bhat * DT
        F = np.broadcast_to(INS, (M, NS, NS)).copy()
        F[:, IT, IB] = -(Rh @ left_jacobian(u)) * DT
        P = F @ P @ np.swapaxes(F, 1, 2)
        P[:, IT, IT] += q_tt * np.eye(3)
        P[:, IB, IB] += q_bb * np.eye(3)
        q = q_mul(q, q_from_rotvec(u))
        q = q / np.linalg.norm(q, axis=1, keepdims=True)
        b_true = b_true + sig_bw * np.sqrt(DT) * rnd.bw[:, k]
        Rh = q_to_R(q)
        Rt = scen.R[k + 1]
        c = scen.c[k + 1]
        # accelerometer
        za = scen.fb[k + 1] + np.sqrt(sig_a2) * rnd.acc[:, k]
        H = np.zeros((M, 3, NS)); H[:, :, IT] = np.swapaxes(Rh, 1, 2) @ skew_g
        q, bhat, hhat, ehat, dhat, P = update(za - np.einsum("mji,j->mi", Rh, gz), H, sig_a2, q, bhat, hhat, ehat, dhat, P)
        Rh = q_to_R(q)
        # truth field
        v_body = Rt.T @ be
        b_body = np.einsum("mij,j->mi", S, v_body) + kappa * (d_hi + c * d_mc)
        if use_mag:
            zm = b_body + cfg["sigma_mag"] * rnd.mag[:, k]
            mh, v = field_model(Rh, be, ehat, hhat, dhat, c)
            H = h_field(Rh, be, ehat, v, c)
            q, bhat, hhat, ehat, dhat, P = update(zm - mh, H, sig_m2, q, bhat, hhat, ehat, dhat, P)
            Rh = q_to_R(q)
        if use_rp:
            ut = b_body / np.linalg.norm(b_body, axis=1, keepdims=True)
            e1, e2 = tangent_basis_rows(ut)
            um = ut + sig_rp * (rnd.rp[:, k, 0, None] * e1 + rnd.rp[:, k, 1, None] * e2)
            um = rnd.sgn[:, k, None] * um / np.linalg.norm(um, axis=1, keepdims=True)
            mh, v = field_model(Rh, be, ehat, hhat, dhat, c)
            mn = np.linalg.norm(mh, axis=1, keepdims=True)
            up = mh / mn
            s = np.where((um * up).sum(1) >= 0, 1.0, -1.0)
            Pm = (np.eye(3) - up[:, :, None] * up[:, None, :]) / mn[:, :, None]
            H = s[:, None, None] * (Pm @ h_field(Rh, be, ehat, v, c))
            q, bhat, hhat, ehat, dhat, P = update(um - s[:, None] * up, H, sig_r2, q, bhat, hhat, ehat, dhat, P)
            Rh = q_to_R(q)
        e = err_state(Rt, Rh)
        yaw_err[:, k + 1] = e[:, 2]
        if (k + 1) in nees_steps:
            x = np.concatenate([e, b_true - bhat], axis=1)
            nees[k + 1] = np.einsum("mi,mi->m", x, np.linalg.solve(P[:, :6, :6], x[:, :, None])[:, :, 0])
    return dict(yaw_err=yaw_err, nees=nees)
