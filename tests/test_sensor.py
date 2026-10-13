import numpy as np
import pytest

from model.sphere import fibonacci_sphere
from model.spin import primary_model
from model.surrogate import PolySurrogate
from sensor.poisson import ELEMENT_SETS, Sensor, hemisphere_directions, tangent_basis, axis_error_deg, RenPhi

B0, K_RATE = 0.05, 1e6


@pytest.fixture(scope="module")
def setup():
    m = primary_model()
    f = lambda v: m.yield_eig(v, B0, K_RATE)
    return m, f, PolySurrogate.fit(f, 8)


def test_poly_matches_eig_holdout(setup):
    m, f, s = setup
    rng = np.random.default_rng(5)
    v = rng.normal(size=(1000, 3)); v /= np.linalg.norm(v, axis=1, keepdims=True)
    assert np.max(np.abs(s(v) - f(v))) < 1e-7


def test_poly_mean_exact(setup):
    m, f, s = setup
    g = fibonacci_sphere(20000)
    assert abs(s.mean - f(g).mean()) < 1e-6


def test_tangent_gradient_matches_finite_difference_of_eigensolver(setup):
    """Independent check: analytic polynomial gradient vs central differences of the eigenbasis yield."""
    m, f, s = setup
    rng = np.random.default_rng(6)
    n = rng.normal(size=(15, 3)); n /= np.linalg.norm(n, axis=1, keepdims=True)
    e1, e2 = tangent_basis(n)
    h = 1e-5
    from model.surrogate import tangent_project
    _, g = s.value_grad(n)
    g = tangent_project(g, n)
    for e in (e1, e2):
        up = n + h * e; up /= np.linalg.norm(up, axis=1, keepdims=True)
        dn = n - h * e; dn /= np.linalg.norm(dn, axis=1, keepdims=True)
        fd = (f(up) - f(dn)) / (2 * h)
        assert np.max(np.abs(fd - (g * e).sum(1))) < 1e-6


def test_fisher_matches_numeric_derivative(setup):
    m, f, s = setup
    sensor = Sensor(s, ELEMENT_SETS["K3orth"], eps=1.0)
    rng = np.random.default_rng(7)
    n = hemisphere_directions(5, rng)
    N = 1e6
    J = sensor.fisher(n, N)
    e1, e2 = tangent_basis(n)
    h = 1e-5
    for b in range(5):
        cols = []
        for e in (e1[b], e2[b]):
            up = n[b] + h * e; up /= np.linalg.norm(up)
            dn = n[b] - h * e; dn /= np.linalg.norm(dn)
            cols.append((sensor.lam(up[None], N)[0] - sensor.lam(dn[None], N)[0]) / (2 * h))
        D = np.stack(cols, axis=1)  # (K,2)
        lam = sensor.lam(n[b:b + 1], N)[0]
        Jnum = (D.T / lam) @ D
        assert np.allclose(J[b], Jnum, rtol=1e-4, atol=1e-6 * np.abs(Jnum).max())


def test_crb_scales_inverse_sqrt_N(setup):
    m, f, s = setup
    sensor = Sensor(s, ELEMENT_SETS["K3orth"], eps=1.0)
    n = hemisphere_directions(20, np.random.default_rng(8))
    r1 = sensor.crb(n, 1e8)[0]
    r4 = sensor.crb(n, 4e8)[0]
    ok = np.isfinite(r1)
    assert np.allclose(r1[ok] / r4[ok], 2.0, rtol=1e-9)


def test_mle_recovers_truth_near_crb(setup):
    m, f, s = setup
    sensor = Sensor(s, ELEMENT_SETS["K3orth"], eps=1.0)
    sensor.make_grid(3000)
    rng = np.random.default_rng(9)
    N = 5e8
    errs, preds = [], []
    for n in hemisphere_directions(40, rng):
        c = sensor.simulate(n, N, rng)
        errs.append(axis_error_deg(sensor.mle(c, N), n))
        preds.append(np.degrees(sensor.crb(n[None], N)[1][0]))
    errs, preds = np.array(errs), np.array(preds)
    ok = np.isfinite(preds) & (preds < 3)
    assert ok.sum() > 10
    ratio = np.median(errs[ok]) / np.median(preds[ok])
    assert 0.5 < ratio < 2.0


def test_ren_phi_symmetry_and_value():
    r = RenPhi()
    v = np.array([[0, 0, 1.0], [0, 0, -1.0], [1.0, 0, 0]])
    y = r(v)
    assert abs(y[0] - y[1]) < 1e-12
    assert abs(y[2] - (0.5 - 0.00127)) < 1e-12
