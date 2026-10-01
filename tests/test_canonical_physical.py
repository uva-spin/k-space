"""Independent implementation checks for the final canonical physical audit."""
import sys
from pathlib import Path

import numpy as np
import pytest
from scipy.linalg import eigh

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'src'))
sys.path.insert(0, str(ROOT/'scripts'))
from canonical_physical import bin_chebyshev_moments
from canonical_independent_reference import TransversePair
from css2_coefficients import qcd_coefficients
from rational_radial import RadialGrid
from resummed_operator import CanonicalRG, apply_chebyshev_resolvent


def test_infinite_canonical_endpoint_is_exact_constant():
    rg = CanonicalRG(qcd_coefficients(5))
    assert rg.tc(np.inf) == 0.
    assert rg.v(np.inf, 20., 400.) == rg.v(2*rg.s0, 20., 400.)
    with pytest.raises(ValueError):
        rg.v(np.nan, 20., 400.)


def test_shared_bin_moments_against_spectral_and_clenshaw_oracles():
    grid = RadialGrid(8., 160)
    edges = np.array([0., .5, 1., 2., 4.])
    pair = TransversePair(.9, 1.3)
    inputs = pair.native(grid.k)
    degree, tau = 70, .08
    coefficients = np.random.default_rng(403).normal(size=(degree+1, 3))/(1+np.arange(degree+1)[:, None])**2
    moments = bin_chebyshev_moments(grid, inputs, edges, degree, tau)
    actual = np.einsum('nbi,nc->bci', moments, coefficients)
    eigenvalues, vectors = eigh(grid.matrix().toarray())
    weights = np.sqrt(grid.weights)
    for c in range(3):
        multiplier = np.polynomial.chebyshev.chebval(2*tau/(tau+eigenvalues)-1, coefficients[:, c])
        evolved = (vectors@(multiplier[:, None]*(vectors.T@(weights[:, None]*inputs))))/weights[:, None]
        expected = grid.bins(evolved, edges)
        np.testing.assert_allclose(actual[:, c], expected, rtol=1e-9, atol=2e-10)
        clenshaw = apply_chebyshev_resolvent(grid.matrix(), weights[:, None]*inputs,
                                           tau, coefficients[:, c])/weights[:, None]
        np.testing.assert_allclose(actual[:, c], grid.bins(clenshaw, edges), rtol=1e-9, atol=2e-10)


@pytest.mark.parametrize('gaussian', [False, True])
def test_two_leg_native_convolution_normalization_and_gradient(gaussian):
    pair = TransversePair(.97, 1.32, gaussian=gaussian)
    edges = np.array([0., .5, 1., 4., 40., 1e8])
    bins = pair.bins(edges)
    assert abs(bins[:, 0].sum()-1) < 2e-11
    assert abs(bins[:, 1].sum()) < 2e-11
    from scipy.integrate import quad_vec
    numerical = np.array([quad_vec(lambda k: 2*np.pi*k*pair.native(np.array([k]))[0],
                                   lo, hi, epsabs=2e-10)[0] for lo, hi in zip(edges[:4], edges[1:5])])
    np.testing.assert_allclose(bins[:4], numerical, atol=3e-9, rtol=1e-8)
    k = np.linspace(0., 8., 51)
    h = 1e-5
    plus = TransversePair(pair.a*np.exp(h), pair.b, gaussian=gaussian).native(k)[:, 0]
    minus = TransversePair(pair.a*np.exp(-h), pair.b, gaussian=gaussian).native(k)[:, 0]
    np.testing.assert_allclose(pair.native(k)[:, 1], (plus-minus)/(2*h), atol=2e-10, rtol=1e-7)


def test_moment_integrator_rejects_fractional_bin_faces():
    grid = RadialGrid(8., 80)
    with pytest.raises(ValueError, match='Bin edges'):
        bin_chebyshev_moments(grid, np.ones((80, 1)), [0., .333, 1.], 10, .2)
