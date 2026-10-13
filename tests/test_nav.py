import numpy as np
from scipy.spatial.transform import Rotation as Rot

from nav import attitude as A


def test_rotation_helpers_match_scipy():
    rng = np.random.default_rng(0)
    v = rng.normal(size=(50, 3)) * 0.8
    assert np.allclose(A.rotvec_to_R(v), Rot.from_rotvec(v).as_matrix(), atol=1e-12)
    assert np.allclose(A.R_to_rotvec(A.rotvec_to_R(v)), v, atol=1e-10)
    q = A.q_from_rotvec(v)
    assert np.allclose(A.q_to_R(q), Rot.from_rotvec(v).as_matrix(), atol=1e-12)
    assert np.allclose(A.q_to_R(A.q_mul(q, q)), A.rotvec_to_R(2 * v), atol=1e-10)


def test_left_jacobian_is_the_derivative_of_exp():
    rng = np.random.default_rng(1)
    for _ in range(5):
        u = rng.normal(size=3) * 0.7
        d = rng.normal(size=3) * 1e-6
        lhs = Rot.from_matrix(A.rotvec_to_R(u + d) @ A.rotvec_to_R(u).T).as_rotvec()
        assert np.allclose(lhs, A.left_jacobian(u) @ d, atol=1e-11)


def test_scenario_respects_pre_registered_bounds():
    s = A.Scenario()
    assert s.max_tilt_deg <= 10.0
    assert s.max_acc <= 1.0
    # three turns of +90, -90, +90 degrees end the turn phase at +90 degrees
    assert abs(np.degrees(s.yaw_true[3000]) - 90.0) < 1e-6


def test_magnetic_aiding_beats_gravity_only_on_yaw():
    sc = A.Scenario()
    rnd = A.Randoms(M=12, T=600)
    cfg = A.default_cfg(soft_iron=False)
    a = A.run_filter("A", cfg, sc, rnd)["yaw_err"][:, 600]
    b = A.run_filter("B", cfg, sc, rnd)["yaw_err"][:, 600]
    assert np.sqrt(np.mean(b**2)) < 0.2 * np.sqrt(np.mean(a**2))


def test_axis_sign_is_irrelevant_to_the_filter():
    """Flipping all random measurement signs must not change the arm-C result (sign-invariant residual)."""
    sc = A.Scenario()
    rnd = A.Randoms(M=6, T=300)
    cfg = A.default_cfg(soft_iron=False, sigma_rp_deg=1.0)
    e1 = A.run_filter("C", cfg, sc, rnd)["yaw_err"]
    rnd.sgn = -rnd.sgn
    e2 = A.run_filter("C", cfg, sc, rnd)["yaw_err"]
    assert np.max(np.abs(e1 - e2)) < 1e-12
