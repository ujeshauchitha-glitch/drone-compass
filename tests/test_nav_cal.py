import numpy as np
from scipy.linalg import expm

from nav import attitude as A
from nav import attitude_cal as C


def test_field_jacobian_matches_central_difference():
    rng = np.random.default_rng(3)
    be = A.earth_field(45.0)
    Rh = A.rotvec_to_R(rng.normal(size=(1, 3)) * 0.5)
    ehat = rng.normal(size=(1, 6)) * 0.02
    hhat = rng.normal(size=(1, 3)); dhat = rng.normal(size=(1, 3)); c = 0.7
    _, v = C.field_model(Rh, be, ehat, hhat, dhat, c)
    H = C.h_field(Rh, be, ehat, v, c)[0]

    def m_of(x):
        Rt = expm(A.skew(x[:3])) @ Rh[0]            # true attitude = exp([dtheta]) R_hat
        e = ehat[0] + x[9:15]
        vv = Rt.T @ be
        return vv + C.e_matrix(e) @ vv + hhat[0] + x[6:9] + c * (dhat[0] + x[15:18])

    num = np.zeros((3, 18))
    h = 1e-6
    for i in range(18):
        d = np.zeros(18); d[i] = h
        num[:, i] = (m_of(d) - m_of(-d)) / (2 * h)
    assert np.max(np.abs(num - H)) < 1e-6


def test_axis_jacobian_matches_central_difference():
    rng = np.random.default_rng(4)
    be = A.earth_field(30.0)
    Rh = A.rotvec_to_R(rng.normal(size=(1, 3)) * 0.4)
    ehat = rng.normal(size=(1, 6)) * 0.02; hhat = rng.normal(size=(1, 3)) * 2; dhat = rng.normal(size=(1, 3)); c = 0.4
    mh, v = C.field_model(Rh, be, ehat, hhat, dhat, c)
    up = mh / np.linalg.norm(mh)
    Pm = (np.eye(3) - up[0][:, None] * up[0][None, :]) / np.linalg.norm(mh)
    H = Pm @ C.h_field(Rh, be, ehat, v, c)[0]

    def u_of(x):
        Rt = expm(A.skew(x[:3])) @ Rh[0]
        vv = Rt.T @ be
        m = vv + C.e_matrix(ehat[0] + x[9:15]) @ vv + hhat[0] + x[6:9] + c * (dhat[0] + x[15:18])
        return m / np.linalg.norm(m)

    num = np.zeros((3, 18)); h = 1e-6
    for i in range(18):
        d = np.zeros(18); d[i] = h
        num[:, i] = (u_of(d) - u_of(-d)) / (2 * h)
    assert np.max(np.abs(num - H)) < 1e-6


def test_cal_filter_reduces_to_uncalibrated_when_priors_vanish():
    sc = A.Scenario()
    rnd = A.Randoms(M=6, T=400)
    cfg = A.default_cfg(soft_iron=False, sigma_d=0.0, sigma_rp_deg=1.0)
    for arm in "BCD":
        e0 = A.run_filter(arm, cfg, sc, rnd)["yaw_err"]
        e1 = C.run_filter_cal(arm, cfg, sc, rnd, prior_scale=1e-9)["yaw_err"]
        assert np.max(np.abs(e0 - e1)) < 1e-6, arm
