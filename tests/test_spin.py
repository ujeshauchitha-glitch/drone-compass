import numpy as np
import pytest

from model import checks
from model.spin import H16, RECALLED


@pytest.mark.parametrize("tensors", [H16, RECALLED])
def test_a_isotropic(tensors):
    assert checks.check_a(tensors=tensors)["pass"]


def test_b_zero_hamiltonian():
    assert checks.check_b()["pass"]


@pytest.mark.parametrize("tensors", [H16, RECALLED])
def test_c_range(tensors):
    assert checks.check_c(tensors=tensors)["pass"]


@pytest.mark.parametrize("tensors", [H16, RECALLED])
def test_d_inversion(tensors):
    assert checks.check_d(tensors=tensors)["pass"]


def test_e_two_implementations_agree():
    assert checks.check_e(n=30)["pass"]


def test_f_convergence():
    assert checks.check_f()["pass"]
