"""Independent CSS2 physical reference and analytic transverse fixtures.

No production RG, matching, flavor-assembly or physical-symbol functions are
imported. Physical coefficient numbers and the official PDF provider are shared
inputs. RG is integrated from the final endpoint with RK45; matching is integrated
in log(xi), distinct from production integration in z.
"""
from __future__ import annotations

import json
from math import gamma, comb, factorial
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp, quad_vec
from scipy.interpolate import CubicSpline
from scipy.special import roots_legendre, kv, j1, expit, gammaincc

BZERO = 2*np.exp(-np.euler_gamma)


class IndependentCanonical:
    def __init__(self, pdf, coefficient_path, Q, xA, xB, mu0=5.,
                 pdf_max=1e5, nodes=769, order=384, rtol=2e-11):
        self.pdf, self.Q, self.xs, self.mu0, self.pdf_max = pdf, Q, (xA, xB), mu0, pdf_max
        self.c = json.loads(Path(coefficient_path).read_text())
        self.order = order
        self.pids = list(range(1, 6)) + list(range(-1, -6, -1))
        def series(a, coefficients, start=1):
            return sum(v*a**(n+start) for n, v in enumerate(coefficients))
        self.series = series
        self.beta = lambda a: -2*series(a, self.c['beta'], 2)
        anchor = solve_ivp(lambda t, y: [self.beta(y[0])],
                           (np.log(91.1876), np.log(Q)), [.1180024/(4*np.pi)],
                           method='RK45', rtol=rtol, atol=rtol*1e-2)
        if not anchor.success:
            raise RuntimeError(anchor.message)
        self.aQ = float(anchor.y[0, -1])
        # I integrates the full virtuality anomalous dimension directly. It
        # does not construct the production G,N,J primitive integrals.
        def rhs(t, y):
            a = y[0]
            return [self.beta(a), series(a, self.c['noncusp'])
                    - 2*np.log(Q/np.exp(t))*series(a, self.c['cusp'])]
        self.down = solve_ivp(rhs, (np.log(Q), np.log(mu0)), [self.aQ, 0.],
                              method='RK45', rtol=rtol, atol=rtol*1e-2, dense_output=True)
        self.up = solve_ivp(rhs, (np.log(Q), np.log(pdf_max)), [self.aQ, 0.],
                            method='RK45', rtol=rtol, atol=rtol*1e-2, dense_output=True)
        if not (self.down.success and self.up.success):
            raise RuntimeError('Reference RG failed')
        # Different knot count and independently evaluated longitudinal integrals.
        self.ts = np.linspace(np.log(mu0), np.log(pdf_max), nodes)
        self.cache = CubicSpline(self.ts, np.array([
            self.direct_blocks(min(pdf_max, np.exp(t))) for t in self.ts]), axis=0)

    def running(self, t):
        t = np.atleast_1d(t)
        out = np.empty((2, len(t)))
        low = t < np.log(self.Q)
        if np.any(low):
            out[:, low] = self.down.sol(t[low])
        if np.any(~low):
            out[:, ~low] = self.up.sol(t[~low])
        return out

    def scale(self, b):
        y = (np.asarray(b)*self.mu0/BZERO)**2
        logmu = np.full_like(y, np.log(self.mu0))
        low = y <= .5
        logmu[low] = np.log(BZERO/np.asarray(b)[low])
        middle = (y > .5) & (y < 1)
        u = 2*y[middle]-1
        step = expit(1/(1-u)-1/u)
        logmu[middle] -= .5*(1-step)*np.log(y[middle])
        return logmu

    def direct_blocks(self, mu, order=None, adaptive=False):
        order = order or self.order
        a = float(self.running(np.log(mu))[0, 0])
        result = []
        for x in self.xs:
            def integrand(u):
                u = np.asarray(u)
                xi = np.exp(u)
                ratio = x/xi
                q = np.array([self.pdf.f(pid, xi, mu) for pid in self.pids])
                g = self.pdf.f(21, xi, mu)
                return np.concatenate((2*(1-ratio)*q,
                    np.broadcast_to(2*ratio*(1-ratio)*g, q.shape)), axis=0)
            if adaptive:
                values, _ = quad_vec(integrand, np.log(x), 0., epsabs=2e-10,
                                     epsrel=2e-10, limit=300)
            else:
                z, w = roots_legendre(order)
                lo = np.log(x)
                nodes = lo+(z+1)*(-lo)/2
                values = integrand(nodes)@(w*(-lo)/2)
            q0 = np.array([self.pdf.f(pid, x, mu) for pid in self.pids])
            qq = a*(4/3)*(values[:10]-np.pi**2*q0/6)
            qg = a*values[10:]
            result.append(np.stack((q0, qq, qg), axis=-1))
        return np.array(result)

    def components(self, b):
        b = np.atleast_1d(b)
        t = self.scale(b)
        if np.max(t) > np.log(self.pdf_max)+1e-10:
            raise ValueError('Reference requests untabulated PDF')
        a, integral = self.running(t)
        d = self.series(a, self.c['d'])
        evolution = np.exp(-2*integral-2*d*(np.log(self.Q)-t))
        blocks = self.cache(t)
        terms = np.zeros((len(b), 6))
        hard = self.aQ*(4/3)*(7*np.pi**2/3-16)
        for i, pid in enumerate(self.pids):
            j = self.pids.index(-pid)
            charge = (2/3)**2 if abs(pid) % 2 == 0 else (-1/3)**2
            A, B = blocks[:, 0, i], blocks[:, 1, j]
            terms[:, 0] += charge*A[:, 0]*B[:, 0]
            terms[:, 2] += charge*A[:, 1]*B[:, 0]
            terms[:, 3] += charge*A[:, 2]*B[:, 0]
            terms[:, 4] += charge*A[:, 0]*B[:, 1]
            terms[:, 5] += charge*A[:, 0]*B[:, 2]
        terms[:, 1] = hard*terms[:, 0]
        return terms*evolution[:, None]


class TransversePair:
    """Two positive normalized Student legs; optional Gaussian pair.

    Native convolution is a one-dimensional Feynman-parameter integral derived
    directly from the two-dimensional momentum convolution. The reference uses
    the product of analytic Bessel-K transforms instead. Derivatives are with
    respect to log(lambda_A), including both normalization and denominator.
    """
    def __init__(self, a, b, p=3, q=4, gaussian=False):
        self.a, self.b, self.p, self.q, self.gaussian = a, b, p, q, gaussian

    def native(self, k, order=160):
        k = np.asarray(k)
        if self.gaussian:
            width = self.a**2+self.b**2
            value = np.exp(-k*k/width)/(np.pi*width)
            derivative = value*2*self.a**2*(k*k/width-1)/width
            return np.stack((value, derivative), axis=-1)
        nodes, weights = roots_legendre(order)
        u, w = (nodes+1)/2, weights/2
        prefactor = ((self.p-1)*(self.q-1)*self.a**(2*self.p-2)
                     * self.b**(2*self.q-2)/np.pi
                     * gamma(self.p+self.q-1)/(gamma(self.p)*gamma(self.q)))
        m = self.p+self.q-1
        result = np.empty((len(k), 2))
        # Bound temporary memory for high-resolution radial grids.
        for start in range(0, len(k), 2048):
            ks = k[start:start+2048, None]
            D = u*self.a**2+(1-u)*self.b**2+u*(1-u)*ks**2
            value = prefactor*u**(self.p-1)*(1-u)**(self.q-1)/D**m
            result[start:start+2048, 0] = value@w
            result[start:start+2048, 1] = (value*((2*self.p-2)-2*m*u*self.a**2/D))@w
        return result

    def bins(self, edges, order=256):
        if self.gaussian:
            width = self.a**2+self.b**2
            lo, hi = edges[:-1], edges[1:]
            value = np.exp(-lo**2/width)-np.exp(-hi**2/width)
            derivative = 2*self.a**2/width**2*(lo**2*np.exp(-lo**2/width)-hi**2*np.exp(-hi**2/width))
            return np.stack((value, derivative), axis=-1)
        z, w = roots_legendre(order)
        u, w = (z+1)/2, w/2
        m = self.p+self.q-2
        factor = ((self.p-1)*(self.q-1)*self.a**(2*self.p-2)*self.b**(2*self.q-2)
                  * gamma(self.p+self.q-1)/(gamma(self.p)*gamma(self.q))/m
                  *u**(self.p-2)*(1-u)**(self.q-2))
        D = u*self.a**2+(1-u)*self.b**2+edges[:, None]**2*u*(1-u)
        primitive = factor/D**m
        derivative = primitive*((2*self.p-2)-2*m*u*self.a**2/D)
        return np.stack(((primitive[:-1]-primitive[1:])@w,
                         (derivative[:-1]-derivative[1:])@w), axis=-1)

    def transforms(self, b):
        if self.gaussian:
            value = np.exp(-(self.a**2+self.b**2)*b*b/4)
            return np.array([value, -self.a**2*b*b/2*value])
        def one(scale, power):
            z = scale*b
            if z == 0:
                return 1., 0.
            coefficient = 2**(2-power)/gamma(power-1)
            return (coefficient*z**(power-1)*kv(power-1, z),
                    -coefficient*z**power*kv(power-2, z))
        va, da = one(self.a, self.p)
        vb, _ = one(self.b, self.q)
        return np.array([va*vb, da*vb])


def continuum_components(reference, pair, edges, eps=2e-10):
    contact = reference.components([2*BZERO/reference.mu0])[0]
    bcut, bmax = BZERO/reference.pdf_max, BZERO/reference.mu0
    lo, hi = edges[:-1], edges[1:]
    def kernel(b):
        return hi*j1(hi*b)-lo*j1(lo*b)
    def transforms(b):
        return np.r_[pair.transforms(b), 1.]
    def integrand(b):
        return ((reference.components([b])[0]-contact)[None, :, None]
                *kernel(b)[:, None, None]*transforms(b)[None, None, :])
    remainder, error = quad_vec(integrand, bcut, bmax,
                                points=[bmax/np.sqrt(2)], epsabs=eps, epsrel=eps, limit=500)
    # Below the last tabulated PDF scale, remove the known contact subtraction
    # exactly. The omitted full M contribution is bounded separately, not hidden
    # in the quadrature error.
    low_contact, low_error = quad_vec(lambda b: kernel(b)[:, None]*transforms(b),
                                     0., bcut, epsabs=eps*1e-3, epsrel=eps)
    input_bins = np.column_stack((pair.bins(edges), (edges[:-1] == 0).astype(float)))
    bins = contact[None, :, None]*(input_bins-low_contact)[:, None, :]+remainder
    return bins, float(error+np.max(np.abs(contact))*low_error)


def ultraviolet_envelope(reference, edges, log_power=4):
    """Conditional PDF-table tail bound for a declared longitudinal class.

    Beyond pdf_max, assume each |f(x,mu)| is bounded by
    |f(x,pdf_max)|*(1+ln(mu/pdf_max))**log_power, including convolution arguments.
    |v| is bounded above by its value at pdf_max once the polynomial sign tests
    in the audit establish monotonic decay. This conservative bound does not
    assume a physical PDF extrapolation or change the RG prescription.
    """
    mu = reference.pdf_max
    blocks = reference.direct_blocks(mu, adaptive=True)
    a, integral = reference.running(np.log(mu))[:, 0]
    ev = np.exp(-2*integral-2*reference.series(a, reference.c['d'])*np.log(reference.Q/mu))
    # Absolute convolution terms are required: cancellation in cqq cannot be
    # used to bound a different ultraviolet continuation.
    bound_blocks = blocks.copy()
    for leg, x in enumerate(reference.xs):
        for i, pid in enumerate(reference.pids):
            q0 = abs(reference.pdf.f(pid, x, mu))
            def integrand(u):
                xi = np.exp(u); z = x/xi
                return np.array([2*(1-z)*abs(reference.pdf.f(pid, xi, mu)),
                                 2*z*(1-z)*abs(reference.pdf.f(21, xi, mu))])
            c, _ = quad_vec(integrand, np.log(x), 0., epsabs=2e-10, epsrel=2e-10)
            bound_blocks[leg, i] = [q0, a*(4/3)*(np.pi**2*q0/6+c[0]), a*c[1]]
    A = np.zeros(6)
    for i, pid in enumerate(reference.pids):
        j = reference.pids.index(-pid)
        q = (2/3)**2 if abs(pid) % 2 == 0 else (1/3)**2
        u, v = bound_blocks[0, i], bound_blocks[1, j]
        A += q*np.array([u[0]*v[0], 0., u[1]*v[0], u[2]*v[0], u[0]*v[1], u[0]*v[2]])
    A[1] = abs(reference.aQ*(4/3)*(7*np.pi**2/3-16))*A[0]
    power = 2*log_power
    # Use the proven monotonicity on small intervals, rather than the very
    # loose constant v(Qmax) over the entire unbounded interval. The RG is
    # well-defined there; no PDF is queried beyond its table.
    t0 = np.log(mu)
    def rhs(u, y):
        aa = y[0]
        return [reference.beta(aa), reference.series(aa, reference.c['noncusp'])
                -2*(np.log(reference.Q)-t0-u)*reference.series(aa, reference.c['cusp'])]
    tail = solve_ivp(rhs, (0., 24.), [a, integral], method='DOP853',
                     rtol=2e-13, atol=2e-15, dense_output=True)
    if not tail.success:
        raise RuntimeError('UV RG bound integration failed')
    u = np.linspace(0., 24., 241)
    aa, ii = tail.sol(u)
    vv = np.exp(-2*ii-2*reference.series(aa, reference.c['d'])*(np.log(reference.Q)-t0-u))
    def remaining_moment(lower):
        return sum(comb(power, j)*factorial(j)/2**(j+1)*gammaincc(j+1, 2*lower)
                   for j in range(power+1))
    # On each interval v^2 <= v^2(left); the last constant also bounds u>24.
    remaining = remaining_moment(u)
    weighted_moment = np.sum(vv[:-1]*(remaining[:-1]-remaining[1:]))+vv[-1]*remaining[-1]
    # Margin is far larger than requested ODE relative accuracy. The report
    # still distinguishes this numerical evaluation from interval arithmetic.
    weighted_moment *= 1+1e-8
    geometry = (edges[1:]**2-edges[:-1]**2)/2*(BZERO/mu)**2
    return geometry[:, None]*A[None, :]*weighted_moment


def resolved_ultraviolet(reference, pair, edges, log_power=0, upper=24.):
    """Resolve the PDF-table remainder for an explicitly specified input.

    For this validation fixture only, each PDF above Qmax is its endpoint value
    times (1+ln(mu/Qmax))**log_power. This is NOT claimed to be DGLAP or QCD
    knowledge outside the source grid. Coupling and canonical RG remain uncapped.
    The production polynomial nodes stay inside the real PDF table; testing
    power=0 and power=4 quantifies sensitivity to this input specification.
    """
    mu = reference.pdf_max
    t0 = np.log(mu)
    a0, i0 = reference.running(t0)[:, 0]
    def rhs(u, y):
        aa = y[0]
        return [reference.beta(aa), reference.series(aa, reference.c['noncusp'])
                -2*(np.log(reference.Q)-t0-u)*reference.series(aa, reference.c['cusp'])]
    ode = solve_ivp(rhs, (0., upper), [a0, i0], method='DOP853',
                    rtol=2e-13, atol=2e-15, dense_output=True)
    if not ode.success:
        raise RuntimeError(ode.message)
    v0 = np.exp(-2*i0-2*reference.series(a0, reference.c['d'])*np.log(reference.Q/mu))
    base = reference.components([BZERO/mu])[0]/v0
    lo, hi = edges[:-1], edges[1:]
    def integrand(u):
        aa, ii = ode.sol(u)
        evolution = np.exp(-2*ii-2*reference.series(aa, reference.c['d'])*(np.log(reference.Q)-t0-u))
        terms = base*evolution*(1+u)**(2*log_power)
        terms[2:] *= aa/a0
        b = BZERO/mu*np.exp(-u)
        kernel = b*(hi*j1(hi*b)-lo*j1(lo*b))
        transforms = np.r_[pair.transforms(b), 1.]
        return kernel[:, None, None]*terms[None, :, None]*transforms[None, None, :]
    value, error = quad_vec(integrand, 0., upper, epsabs=1e-15, epsrel=2e-11)
    # Above upper, all component constants keep their signs in these fixtures.
    # a and v^2 decrease. Their values at upper bound the remaining magnitude.
    aa, ii = ode.sol(upper)
    vv = np.exp(-2*ii-2*reference.series(aa, reference.c['d'])*(np.log(reference.Q)-t0-upper))
    constants = abs(base)*vv
    constants[2:] *= aa/a0
    power = 2*log_power
    moment = sum(comb(power, j)*factorial(j)/2**(j+1)*gammaincc(j+1, 2*upper)
                 for j in range(power+1))
    bound = (hi**2-lo**2)[:, None]/2*(BZERO/mu)**2*constants[None, :]*moment
    return value, float(error), bound
