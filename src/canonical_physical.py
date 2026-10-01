"""Uncapped CSS2 canonical two-leg strict-NLO physical symbol.

Matching is canonical at the completed coordinate b_eval=b0/mu_c(s).
There is no Q-dependent scale cap. The numerical PDF provider is never called
outside its tabulated domain. The analytic zero endpoint follows under the
declared polynomial-logarithmic growth assumption for longitudinal boundaries.
"""
from __future__ import annotations

import numpy as np
from scipy.fft import dct
from scipy.integrate import solve_ivp
from scipy.interpolate import CubicSpline
from scipy.sparse import eye
from scipy.sparse.linalg import splu

from css2_coefficients import qcd_coefficients
from css2_matching import match_nlo, hard_nlo
from resummed_operator import CanonicalRG


COMPONENTS = ('born', 'hard_1', 'ope_qq_leg_a', 'ope_qg_leg_a',
              'ope_qq_leg_b', 'ope_qg_leg_b')
FLAVORS = (-5, -4, -3, -2, -1, 1, 2, 3, 4, 5)


def physical_rg(mu0=5., alpha_mz=.1180024, mz=91.1876):
    c = qcd_coefficients(5)
    sol = solve_ivp(lambda t, y: [-2*c.series(y[0], c.beta, start=2)],
                    (np.log(mz), np.log(mu0)), [alpha_mz/(4*np.pi)],
                    method='DOP853', rtol=2e-13, atol=2e-15)
    if not sol.success:
        raise RuntimeError(sol.message)
    return CanonicalRG(c, mu0, float(sol.y[0, -1]))


class CanonicalPhysicalSymbol:
    """Shared-flavor transverse boundary, with physical longitudinal matching.

    The symbol acts on the native convolution of the two transverse inputs.
    The six returned components are strictly linear in hard/OPE insertions.
    Fixed nf=5 is a massless diagnostic, including above the top mass; it does
    not claim a massive Standard Model threshold treatment.
    """
    def __init__(self, pdf, rg, Q, xA, xB, nodes=513, order=256, pdf_max=1e5):
        if not (Q >= rg.mu0 and 0 < xA < 1 and 0 < xB < 1 and pdf_max > Q):
            raise ValueError('Invalid physical canonical benchmark kinematics')
        self.pdf, self.rg, self.Q = pdf, rg, Q
        self.xA, self.xB, self.order, self.pdf_max = xA, xB, order, pdf_max
        self.ts = np.linspace(0., np.log(pdf_max/rg.mu0), nodes)
        blocks = np.array([self.direct_blocks(min(pdf_max, rg.mu0*np.exp(t)))
                           for t in self.ts])
        self.cache = CubicSpline(self.ts, blocks, axis=0, extrapolate=False)
        aq = float(rg.values(np.log(Q/rg.mu0))[0])
        self.hard = hard_nlo(aq, Q, Q)

    def direct_blocks(self, mu, order=None):
        a = float(self.rg.values(np.log(mu/self.rg.mu0))[0])
        return np.array([[[m.born, m.qq, m.qg] for pid in FLAVORS
                          for m in [match_nlo(x, lambda y: self.pdf.f(pid, y, mu),
                                              lambda y: self.pdf.f(21, y, mu), a,
                                              order=order or self.order)]]
                         for x in (self.xA, self.xB)])

    def components(self, s):
        s = np.asarray(s, float)
        if np.isnan(s).any() or (s < 0).any():
            raise ValueError('Require nonnegative spectral coordinate')
        shape = s.shape
        flat = s.ravel()
        out = np.zeros((flat.size, 6))
        positive = flat > 0
        if np.any(positive):
            tc = self.rg.tc(flat[positive])
            if np.any(tc > self.ts[-1]+1e-12):
                raise ValueError('Canonical scale exceeds PDF table; no hidden cap or extrapolation')
            blocks = self.cache(np.minimum(tc, self.ts[-1]))
            A, B = blocks[:, 0], blocks[:, 1, ::-1]
            charges = np.array([4/9 if abs(pid) in (2, 4) else 1/9 for pid in FLAVORS])
            product = lambda i, j: np.sum(charges*A[:, :, i]*B[:, :, j], axis=-1)
            born = product(0, 0)
            terms = np.stack([born, self.hard*born, product(1, 0), product(2, 0),
                              product(0, 1), product(0, 2)], axis=-1)
            evolution = np.exp(2*self.rg.log_v_tc(tc, self.Q, self.Q**2))
            out[positive] = terms*evolution[:, None]
        return out.reshape(shape+(6,))

    def polynomials(self, degree, tau):
        x = np.cos(np.pi*np.arange(degree+1)/degree)
        r = (x+1)/2
        s = np.empty_like(r)
        s[0], s[-1] = 0., np.inf
        s[1:-1] = tau*(1-r[1:-1])/r[1:-1]
        values = self.components(s)
        c = dct(values, type=1, axis=0)/degree
        c[[0, -1]] *= .5
        return c


def bin_chebyshev_moments(grid, inputs, edges, degree, tau):
    """Native finite-bin moments of Chebyshev polynomials of the resolvent.

    Shared recurrence avoids repeating elliptic solves for each physical
    component. Only three input-shaped arrays and the small bin-moment table
    are retained. No transform of the input is performed.
    """
    values = np.asarray(inputs, float)
    if values.ndim != 2 or len(values) != grid.cells or degree < 1 or tau <= 0:
        raise ValueError('Invalid input shape or polynomial parameters')
    edges = np.asarray(edges, float)
    # Bin boundaries must coincide with grid faces for this integrator.
    if np.max(np.abs(edges/(grid.radius/grid.cells)-np.round(edges/(grid.radius/grid.cells)))) > 1e-7:
        raise ValueError('Bin edges must coincide with radial grid faces')
    weights = np.sqrt(grid.weights)
    masks = [(grid.k >= lo) & (grid.k < hi) for lo, hi in zip(edges[:-1], edges[1:])]
    def bins(v):
        return np.array([np.sum(weights[m, None]*v[m], axis=0) for m in masks])
    solve = splu(grid.matrix()+tau*eye(grid.cells, format='csc'))
    def X(v):
        return 2*tau*solve.solve(v)-v
    b0 = weights[:, None]*values
    b1 = X(b0)
    result = np.empty((degree+1, len(edges)-1, values.shape[1]))
    result[0], result[1] = bins(b0), bins(b1)
    for n in range(2, degree+1):
        b0, b1 = b1, 2*X(b1)-b0
        result[n] = bins(b1)
    return result
