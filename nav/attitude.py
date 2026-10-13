"""Stage 3: attitude MEKF with magnetic aiding (PLAN.md section 5 + Amendment 2).

Batched over Monte-Carlo runs (numpy).  Error state: R_true = exp([dtheta]x) R_hat (nav-frame, left-multiplicative),
db = b_true - b_hat.  Storage: unit quaternion.  All sensor parameters are ASSUMED (see PLAN.md).
"""
from __future__ import annotations

import numpy as np
from scipy.spatial.transform import Rotation as Rot

G = 9.81
DT = 0.1
T_TOTAL = 600.0
N_STEP = 6000
DEG = np.pi / 180.0
SEED = 20261013
B_EARTH_UT = 50.0


# ----------------------------------------------------------------------------------------------- rotation helpers
def skew(v):
    v = np.asarray(v, float)
    z = np.zeros(v.shape[:-1])
    return np.stack([np.stack([z, -v[..., 2], v[..., 1]], -1),
                     np.stack([v[..., 2], z, -v[..., 0]], -1),
                     np.stack([-v[..., 1], v[..., 0], z], -1)], -2)


def rotvec_to_R(v):
    v = np.asarray(v, float)
    th = np.linalg.norm(v, axis=-1)[..., None, None]
    K = skew(v)
    small = th < 1e-8
    ths = np.where(small, 1.0, th)
    A = np.where(small, 1.0 - th**2 / 6, np.sin(ths) / ths)
    B = np.where(small, 0.5 - th**2 / 24, (1 - np.cos(ths)) / ths**2)
    return np.eye(3) + A * K + B * (K @ K)


def R_to_rotvec(R):
    R = np.asarray(R, float)
    skewpart = 0.5 * (R - np.swapaxes(R, -1, -2))
    v = np.stack([skewpart[..., 2, 1], skewpart[..., 0, 2], skewpart[..., 1, 0]], -1)
    s = np.linalg.norm(v, axis=-1)
    c = 0.5 * (np.trace(R, axis1=-2, axis2=-1) - 1.0)
    th = np.arctan2(s, c)
    f = np.where(s > 1e-12, th / np.where(s > 1e-12, s, 1.0), 1.0)
    return v * f[..., None]


def q_from_rotvec(v):
    v = np.asarray(v, float)
    th = np.linalg.norm(v, axis=-1)
    half = 0.5 * th
    k = np.where(th < 1e-8, 0.5 - th**2 / 48, np.sin(half) / np.where(th < 1e-8, 1.0, th))
    return np.concatenate([np.cos(half)[..., None], v * k[..., None]], -1)


def q_mul(a, b):
    w1, x1, y1, z1 = a[..., 0], a[..., 1], a[..., 2], a[..., 3]
    w2, x2, y2, z2 = b[..., 0], b[..., 1], b[..., 2], b[..., 3]
    return np.stack([w1 * w2 - x1 * x2 - y1 * y2 - z1 * z2,
                     w1 * x2 + x1 * w2 + y1 * z2 - z1 * y2,
                     w1 * y2 - x1 * z2 + y1 * w2 + z1 * x2,
                     w1 * z2 + x1 * y2 - y1 * x2 + z1 * w2], -1)


def q_to_R(q):
    w, x, y, z = q[..., 0], q[..., 1], q[..., 2], q[..., 3]
    return np.stack([np.stack([1 - 2 * (y * y + z * z), 2 * (x * y - z * w), 2 * (x * z + y * w)], -1),
                     np.stack([2 * (x * y + z * w), 1 - 2 * (x * x + z * z), 2 * (y * z - x * w)], -1),
                     np.stack([2 * (x * z - y * w), 2 * (y * z + x * w), 1 - 2 * (x * x + y * y)], -1)], -2)


def R_to_q(R):
    xyzw = Rot.from_matrix(R).as_quat()
    return np.concatenate([xyzw[..., 3:4], xyzw[..., :3]], -1)


def left_jacobian(u):
    """J_l(u) = I + (1-cos th)/th^2 [u]x + (th - sin th)/th^3 [u]x^2  (exact linearisation of the bias coupling)."""
    th = np.linalg.norm(u, axis=-1)[..., None, None]
    K = skew(u)
    small = th < 1e-6
    ths = np.where(small, 1.0, th)
    a = np.where(small, 0.5 - th**2 / 24, (1 - np.cos(ths)) / ths**2)
    b = np.where(small, 1 / 6 - th**2 / 120, (ths - np.sin(ths)) / ths**3)
    return np.eye(3) + a * K + b * (K @ K)


def Rz(a):
    c, s = np.cos(a), np.sin(a)
    o, z = np.ones_like(a), np.zeros_like(a)
    return np.stack([np.stack([c, -s, z], -1), np.stack([s, c, z], -1), np.stack([z, z, o], -1)], -2)


def Ry(a):
    c, s = np.cos(a), np.sin(a)
    o, z = np.ones_like(a), np.zeros_like(a)
    return np.stack([np.stack([c, z, s], -1), np.stack([z, o, z], -1), np.stack([-s, z, c], -1)], -2)


def Rx(a):
    c, s = np.cos(a), np.sin(a)
    o, z = np.ones_like(a), np.zeros_like(a)
    return np.stack([np.stack([o, z, z], -1), np.stack([z, c, -s], -1), np.stack([z, s, c], -1)], -2)


def smoothstep(u):
    u = np.clip(u, 0.0, 1.0)
    return u * u * (3 - 2 * u)


# ----------------------------------------------------------------------------------------------- truth scenario
class Scenario:
    """Deterministic trajectory of Amendment 2 A2.1, 100 Hz truth, sampled at 10 Hz."""

    def __init__(self, a_scale=1.0):
        t = np.arange(0.0, T_TOTAL + 1e-9, 0.01)
        n = len(t)
        # roll / pitch
        ramp = smoothstep((t - 300.0) / 10.0)
        amp = (3.0 + 2.0 * ramp) * DEG
        fr = 0.3 + (0.2 - 0.3) * ramp
        fp = 0.3 + (0.13 - 0.3) * ramp
        ph_r = np.concatenate([[0.0], np.cumsum(2 * np.pi * 0.5 * (fr[1:] + fr[:-1]) * 0.01)])
        ph_p = np.concatenate([[0.0], np.cumsum(2 * np.pi * 0.5 * (fp[1:] + fp[:-1]) * 0.01)]) + np.pi / 2
        roll, pitch = amp * np.sin(ph_r), amp * np.sin(ph_p)
        # yaw
        yaw = np.zeros(n)
        for tk, sg in ((130.0, 1.0), (190.0, -1.0), (250.0, 1.0)):
            u = np.clip((t - tk) / 18.0, 0.0, 1.0)
            yaw += sg * (np.pi / 2) * (u - np.sin(2 * np.pi * u) / (2 * np.pi))
        A, om = 5.0, 2 * np.pi / 90.0
        tau = np.clip(t - 300.0, 0.0, None)
        fig = t >= 300.0
        vx, vy = A * om * np.cos(om * tau), A * om * np.cos(2 * om * tau)
        head = np.unwrap(np.arctan2(vy, vx))
        yaw = np.where(fig, yaw[np.searchsorted(t, 300.0)] + (head - head[np.searchsorted(t, 300.0)]), yaw)
        ax = np.where(fig, -A * om**2 * np.sin(om * tau), 0.0)
        ay = np.where(fig, -2 * A * om**2 * np.sin(2 * om * tau), 0.0)
        a_nav = a_scale * np.stack([ax, ay, np.zeros(n)], -1)
        R100 = Rz(yaw) @ Ry(pitch) @ Rx(roll)
        yawrate = np.gradient(yaw, 0.01)
        c = np.clip(0.5 + 0.3 * np.sin(2 * np.pi * 0.05 * t) + 0.2 * np.abs(yawrate) / (10 * DEG), 0.0, 1.0)
        self.t100 = t
        self.R = R100[::10]                                    # (6001, 3, 3)
        dR = np.swapaxes(self.R[:-1], 1, 2) @ self.R[1:]
        self.dtheta = R_to_rotvec(dR)                          # (6000, 3) body delta-angle
        self.fb = np.einsum("kji,kj->ki", self.R, a_nav[::10] + np.array([0, 0, G]))  # specific force, body frame
        self.c = c[::10]
        self.max_yaw_rate_dps = float(np.max(np.abs(yawrate)) / DEG)
        self.max_acc = float(np.max(np.linalg.norm(a_nav, axis=1)))
        self.max_tilt_deg = float(np.max(np.degrees(np.arccos(np.clip(R100[:, 2, 2], -1, 1)))))
        self.yaw_true = yaw[::10]


# ----------------------------------------------------------------------------------------------- common random numbers
class Randoms:
    """Standard-normal inputs indexed by run only (paired across arms and parameter settings)."""

    def __init__(self, M=200, T=N_STEP, seed=SEED):
        self.M, self.T = M, T
        self.gyro = np.empty((M, T, 3)); self.bw = np.empty((M, T, 3)); self.acc = np.empty((M, T, 3))
        self.mag = np.empty((M, T, 3)); self.rp = np.empty((M, T, 2)); self.sgn = np.empty((M, T))
        self.b0 = np.empty((M, 3)); self.dhi = np.empty((M, 3)); self.dmc = np.empty((M, 3))
        self.E = np.empty((M, 6)); self.init = np.empty((M, 3))
        for r in range(M):
            g = np.random.default_rng(np.random.SeedSequence(seed, spawn_key=(10, r)))
            self.gyro[r] = g.standard_normal((T, 3)); self.bw[r] = g.standard_normal((T, 3))
            self.acc[r] = g.standard_normal((T, 3)); self.mag[r] = g.standard_normal((T, 3))
            self.rp[r] = g.standard_normal((T, 2)); self.sgn[r] = np.where(g.random(T) < 0.5, -1.0, 1.0)
            self.b0[r] = g.standard_normal(3); self.dhi[r] = g.standard_normal(3); self.dmc[r] = g.standard_normal(3)
            self.E[r] = g.standard_normal(6); self.init[r] = g.standard_normal(3)

    def soft_iron(self, scale=0.01):
        e = self.E * scale
        S = np.empty((self.M, 3, 3))
        S[:] = np.eye(3)
        S[:, 0, 0] += e[:, 0]; S[:, 1, 1] += e[:, 1]; S[:, 2, 2] += e[:, 2]
        S[:, 0, 1] += e[:, 3]; S[:, 1, 0] += e[:, 3]
        S[:, 0, 2] += e[:, 4]; S[:, 2, 0] += e[:, 4]
        S[:, 1, 2] += e[:, 5]; S[:, 2, 1] += e[:, 5]
        return S


def default_cfg(**kw):
    c = dict(bias_dph=10.0, arw=0.3, sigma_d=0.0, dip=45.0, sigma_mag=0.5, sigma_rp_deg=2.0, kappa=1.0,
             soft_iron=True, soft_iron_scale=0.01, init_yaw_deg=10.0, vib_sigma=0.3)
    c.update(kw)
    if "acc_sigma" not in c:
        c["acc_sigma"] = float(np.hypot(0.05, c["vib_sigma"]))
    return c


def earth_field(dip_deg):
    d = dip_deg * DEG
    return B_EARTH_UT * np.array([0.0, np.cos(d), -np.sin(d)])


def tangent_basis_rows(u):
    """(M,3) unit vectors -> two orthonormal tangent vectors (deterministic, same rule as the sensor code)."""
    a = np.where((np.abs(u[:, 0]) < 0.9)[:, None], np.array([1.0, 0, 0]), np.array([0, 1.0, 0]))
    e1 = np.cross(u, a)
    e1 /= np.linalg.norm(e1, axis=1, keepdims=True)
    e2 = np.cross(u, e1)
    return e1, e2


# ----------------------------------------------------------------------------------------------- batched MEKF
def run_filter(arm, cfg, scen, rnd, nees_steps=(), progress=None):
    """arm in {'A','B','C','D'}.  Returns dict(yaw_err (M, T+1) [rad], nees {step: (M,)}, tilt_err_final)."""
    M, T = rnd.M, rnd.T
    use_mag = arm in ("B", "D")
    use_rp = arm in ("C", "D")
    be = earth_field(cfg["dip"])
    bnorm = np.linalg.norm(be)
    un = be / bnorm
    sd = cfg["sigma_d"]
    bias_rps = cfg["bias_dph"] / 3600.0 * DEG
    arw = cfg["arw"] / 60.0 * DEG                      # rad / sqrt(s)
    sig_bw = 0.1 * bias_rps / np.sqrt(600.0)
    q_tt = arw**2 * DT
    q_bb = sig_bw**2 * DT
    sig_a2 = cfg["acc_sigma"] ** 2
    si_R = cfg.get("si_R_um", 0.0)                      # EXPLORATORY: per-axis soft-iron std (uT) added to the filter's R (0 = pre-registered)
    sig_m2 = cfg["sigma_mag"] ** 2 + sd**2 / 3.0 + si_R**2
    sig_rp = cfg["sigma_rp_deg"] * DEG
    sig_r2 = sig_rp**2 + sd**2 / (3.0 * bnorm**2) + (si_R / bnorm) ** 2
    S = rnd.soft_iron(cfg.get("soft_iron_scale", 0.01)) if cfg["soft_iron"] else np.broadcast_to(np.eye(3), (M, 3, 3)).copy()
    d_hi = sd / np.sqrt(6.0) * rnd.dhi
    d_mc = sd / np.sqrt(6.0) * rnd.dmc
    kappa = cfg["kappa"]
    # truth initial / filter initial
    dth0 = rnd.init * np.array([1.0, 1.0, cfg["init_yaw_deg"]]) * DEG
    R0 = scen.R[0]
    Rhat0 = rotvec_to_R(-dth0) @ R0
    q = R_to_q(Rhat0)
    bhat = np.zeros((M, 3))
    b_true = bias_rps * rnd.b0
    P = np.zeros((M, 6, 6))
    P[:, :3, :3] = np.diag(np.array([1.0, 1.0, cfg["init_yaw_deg"]]) * DEG) ** 2
    P[:, 3:, 3:] = np.eye(3) * bias_rps**2
    I6 = np.eye(6)
    gz = np.array([0.0, 0.0, G])
    skew_g, skew_be, skew_un = skew(gz), skew(be), skew(un)
    yaw_err = np.empty((M, T + 1))
    nees = {}

    def err_state(Rt, Rh):
        return R_to_rotvec(Rt @ np.swapaxes(Rh, 1, 2))

    Rh = q_to_R(q)
    yaw_err[:, 0] = err_state(R0, Rh)[:, 2]

    def update(r, Hth, Rm, q, bhat, P):
        H = np.concatenate([Hth, np.zeros_like(Hth)], axis=2)
        PHt = P @ np.swapaxes(H, 1, 2)
        Sm = H @ PHt + Rm * np.eye(3)
        K = np.swapaxes(np.linalg.solve(Sm, np.swapaxes(PHt, 1, 2)), 1, 2)
        dx = np.einsum("mij,mj->mi", K, r)
        IKH = I6 - K @ H
        P = IKH @ P @ np.swapaxes(IKH, 1, 2) + Rm * (K @ np.swapaxes(K, 1, 2))
        q = q_mul(q_from_rotvec(dx[:, :3]), q)
        q = q / np.linalg.norm(q, axis=1, keepdims=True)
        return q, bhat + dx[:, 3:], P

    for k in range(T):
        # ---- truth gyro
        dtm = scen.dtheta[k] + b_true * DT + arw * np.sqrt(DT) * rnd.gyro[:, k]
        u = dtm - bhat * DT
        # ---- covariance propagation (uses R_hat and the increment before propagation)
        F = np.broadcast_to(I6, (M, 6, 6)).copy()
        F[:, :3, 3:] = -(Rh @ left_jacobian(u)) * DT
        P = F @ P @ np.swapaxes(F, 1, 2)
        P[:, :3, :3] += q_tt * np.eye(3)
        P[:, 3:, 3:] += q_bb * np.eye(3)
        q = q_mul(q, q_from_rotvec(u))
        q = q / np.linalg.norm(q, axis=1, keepdims=True)
        # ---- truth bias random walk to t_{k+1}
        b_true = b_true + sig_bw * np.sqrt(DT) * rnd.bw[:, k]
        Rh = q_to_R(q)
        Rt = scen.R[k + 1]
        # ---- accelerometer (gravity direction)
        za = scen.fb[k + 1] + np.sqrt(sig_a2) * rnd.acc[:, k]
        pred = np.einsum("mji,j->mi", Rh, gz)
        Hth = np.swapaxes(Rh, 1, 2) @ skew_g
        q, bhat, P = update(za - pred, Hth, sig_a2, q, bhat, P)
        Rh = q_to_R(q)
        # ---- total body-frame field (shared by both magnetic sensors)
        v_body = Rt.T @ be                                   # (3,)
        sv = np.einsum("mij,j->mi", S, v_body)
        dist = d_hi + scen.c[k + 1] * d_mc
        b_mems = sv + dist                                   # the MEMS magnetometer always sees the full disturbance
        b_body = sv + kappa * dist                           # the radical-pair sensor sees kappa x disturbance (kappa = 1 by default)
        if use_mag:
            zm = b_mems + cfg["sigma_mag"] * rnd.mag[:, k]
            predm = np.einsum("mji,j->mi", Rh, be)
            Hth = np.swapaxes(Rh, 1, 2) @ skew_be
            q, bhat, P = update(zm - predm, Hth, sig_m2, q, bhat, P)
            Rh = q_to_R(q)
        if use_rp:
            ut = b_body / np.linalg.norm(b_body, axis=1, keepdims=True)
            e1, e2 = tangent_basis_rows(ut)
            um = ut + sig_rp * (rnd.rp[:, k, 0, None] * e1 + rnd.rp[:, k, 1, None] * e2)
            um = rnd.sgn[:, k, None] * um / np.linalg.norm(um, axis=1, keepdims=True)
            upred = np.einsum("mji,j->mi", Rh, un)
            s = np.where((um * upred).sum(1) >= 0, 1.0, -1.0)
            Hth = s[:, None, None] * (np.swapaxes(Rh, 1, 2) @ skew_un)
            q, bhat, P = update(um - s[:, None] * upred, Hth, sig_r2, q, bhat, P)
            Rh = q_to_R(q)
        e = err_state(Rt, Rh)
        yaw_err[:, k + 1] = e[:, 2]
        if (k + 1) in nees_steps:
            x = np.concatenate([e, b_true - bhat], axis=1)
            nees[k + 1] = np.einsum("mi,mi->m", x, np.linalg.solve(P, x[:, :, None])[:, :, 0])
        if progress and (k + 1) % 1000 == 0:
            progress(k + 1)
    return dict(yaw_err=yaw_err, nees=nees, tilt_err=np.linalg.norm(e[:, :2], axis=1))


def yaw_rmse(yaw_err, t_s):
    k = int(round(t_s / DT))
    e = (yaw_err[:, k] + np.pi) % (2 * np.pi) - np.pi
    return float(np.sqrt(np.mean(e**2)))
