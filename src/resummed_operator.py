"""CSS2 canonical RG symbols and momentum-space resolvent evaluation.

No Fourier/Hankel transform is used here. Spectral scalar nodes are used to
construct polynomials of the momentum Laplacian; input arrays stay in momentum
space. Both a lower-order demonstration input and a sourced physical N3LL
evolution input are shipped. Neither is a physical N3LL-prime boundary payload.
"""
from __future__ import annotations
from dataclasses import dataclass
import json
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
from scipy.fft import dct
from scipy.sparse import csc_matrix, eye
from scipy.sparse.linalg import splu

B0 = 2.0 * np.exp(-np.euler_gamma)

@dataclass(frozen=True)
class Coefficients:
    beta: tuple[float, ...]
    cusp: tuple[float, ...]
    noncusp: tuple[float, ...]
    d: tuple[float, ...]  # d1, d2, d3; d1 must be zero in this implementation.
    description: str
    physical_n3llprime_payload: bool = False
    scheme: str | None = None
    evolution_accuracy: str | None = None
    hard_and_ope_order: str | None = None
    nf: int | None = None
    source_manifest: str | None = None

    @classmethod
    def load(cls, path: str | Path) -> "Coefficients":
        data = json.loads(Path(path).read_text())
        obj = cls(*(tuple(float(v) for v in data[k]) for k in
                    ("beta", "cusp", "noncusp", "d")),
                  description=data["description"],
                  physical_n3llprime_payload=bool(data.get("physical_n3llprime_payload", False)),
                  scheme=data.get("scheme"),
                  evolution_accuracy=data.get("evolution_accuracy"),
                  hard_and_ope_order=data.get("hard_and_ope_order"),
                  nf=data.get("nf"),
                  source_manifest=data.get("source_manifest"))
        for name, vals, limit in (("beta", obj.beta, 4), ("cusp", obj.cusp, 4),
                                 ("noncusp", obj.noncusp, 3), ("d", obj.d, 3)):
            if not 1 <= len(vals) <= limit or not np.isfinite(vals).all():
                raise ValueError(f"Invalid {name} coefficients")
        if obj.beta[0] <= 0 or obj.cusp[0] <= 0 or obj.d[0] != 0:
            raise ValueError("Requires beta0>0, Gamma0>0, and canonical d1=0")
        if obj.evolution_accuracy == "N3LL" and tuple(map(len,
                (obj.beta, obj.cusp, obj.noncusp, obj.d))) != (4, 4, 3, 3):
            raise ValueError("N3LL evolution metadata requires 4/4/3/3 coefficients")
        if obj.nf is not None and (isinstance(obj.nf, bool) or obj.nf not in (3, 4, 5, 6)):
            raise ValueError("Unsupported fixed-flavor metadata")
        return obj

    @staticmethod
    def series(a: np.ndarray | float, coeff: tuple[float, ...], start: int = 1):
        return np.polynomial.polynomial.polyval(a, (0.0,) * start + coeff)

class CanonicalRG:
    """Unexpanded ODE solution; t=ln(nu/mu0) avoids ultraviolet overflow."""
    def __init__(self, coeff: Coefficients, mu0: float = 2.0,
                 a0: float = 0.02, tmax: float = 1200.0,
                 rtol: float = 2e-12, atol: float = 2e-14):
        if mu0 <= 0 or a0 <= 0 or tmax <= 0:
            raise ValueError("Positive mu0, a0 and tmax are required")
        self.c, self.mu0, self.a0, self.tmax = coeff, mu0, a0, tmax
        self.s0 = B0**2 / mu0**2
        # Verify positivity of the beta polynomial throughout [0,a0], by
        # checking endpoints and all real stationary points (not just a grid).
        pol = np.polynomial.Polynomial(coeff.beta)
        pts = [0., a0]
        pts += [float(z.real) for z in pol.deriv().roots()
                if abs(z.imag) < 1e-10 and 0 < z.real < a0]
        if min(pol(pts)) <= 0:
            raise ValueError("Selected a0 is not on the required AF branch")
        self.beta_min = float(min(pol(pts)))
        def rhs(t, y):
            a = y[0]
            ga = coeff.series(a, coeff.cusp)
            return [-2*coeff.series(a, coeff.beta, start=2),
                    ga, coeff.series(a, coeff.noncusp), t*ga]
        sol = solve_ivp(rhs, (0., tmax), (a0, 0., 0., 0.),
                        method="DOP853", rtol=rtol, atol=atol, dense_output=True)
        if not sol.success or np.min(sol.y[0]) <= 0:
            raise RuntimeError(f"RG integration failed: {sol.message}")
        self.sol = sol.sol

    def values(self, t):
        t = np.asarray(t, dtype=float)
        if np.any(t < -1e-13) or np.any(t > self.tmax):
            raise ValueError("Requested scale outside the solved RG interval")
        return self.sol(np.maximum(t, 0.))

    @staticmethod
    def chi(u):
        u = np.asarray(u, dtype=float)
        out = np.zeros_like(u)
        out[u >= 1] = 1
        sel = (u > 0) & (u < 1)
        x = u[sel]
        # stable logistic form of exp(-1/u)/(exp(-1/u)+exp(-1/(1-u)))
        from scipy.special import expit
        out[sel] = expit(-1/x + 1/(1-x))
        return out

    def tc(self, s):
        s = np.asarray(s, dtype=float)
        if np.isnan(s).any() or np.any(s <= 0):
            raise ValueError("Use analytic endpoint limits at s=0")
        # The exact constant infrared endpoint includes positive infinity.
        # Avoid 0*log(infinity), which otherwise produces NaN at that endpoint.
        out = np.zeros_like(s)
        sel = s < self.s0
        y = s[sel]/self.s0
        out[sel] = -0.5*(1-self.chi(2*y-1))*np.log(y)
        return out

    def endpoint(self, mu: float, zeta: float):
        if mu < self.mu0 or zeta <= 0:
            raise ValueError("mu>=mu0 and zeta>0 are required")
        return np.log(mu/self.mu0), np.log(np.sqrt(zeta)/self.mu0)

    def D_tc(self, tc, mu):
        tf = np.log(mu/self.mu0)
        ac, gc, _, _ = self.values(tc)
        _, gf, _, _ = self.values(tf)
        return self.c.series(ac, self.c.d) + 2*(gf-gc)

    def D(self, s, mu):
        return self.D_tc(self.tc(s), mu)

    def log_v_tc(self, tc, mu, zeta):
        tf, w = self.endpoint(mu, zeta)
        ac, gc, nc, jc = self.values(tc)
        _, gf, nf, jf = self.values(tf)
        return (-self.c.series(ac, self.c.d)*(w-tc)
                + nf-nc - 2*w*(gf-gc) + 2*(jf-jc))

    def v(self, s, mu, zeta):
        s = np.asarray(s, dtype=float)
        if np.isnan(s).any() or np.any(s < 0):
            raise ValueError("s must be nonnegative")
        out = np.zeros_like(s)
        sel = s > 0
        if not np.any(sel):
            return out
        out[sel] = np.exp(self.log_v_tc(self.tc(s[sel]), mu, zeta))
        return out

    def log_u(self, s, target: tuple[float,float], initial: tuple[float,float]):
        tf, wf = self.endpoint(*target)
        ti, wi = self.endpoint(*initial)
        _, gf, nf, jf = self.values(tf)
        _, gi, ni, ji = self.values(ti)
        R = nf-ni - 2*wi*(gf-gi)+2*(jf-ji)
        return R-(wf-wi)*self.D(s, target[0])

    def h_nodes(self, degree: int, tau: float, mu: float, zeta: float):
        if degree < 1 or tau <= 0:
            raise ValueError("degree>=1 and tau>0 are required")
        r = np.linspace(0., 1., degree+1)
        minf = float(self.v(self.s0*2, mu, zeta))
        h = np.zeros_like(r)
        h[-1] = -minf
        sel = (r > 0) & (r < 1)
        s = tau*(1-r[sel])/r[sel]
        h[sel] = (self.v(s, mu, zeta)-minf)/r[sel]
        return minf, h

    def chebyshev_coefficients(self, degree: int, tau: float,
                               mu: float, zeta: float):
        """Lobatto interpolant of v(tau*(1-r)/r), r in [0,1].

        The endpoints r=0 (s=infinity) and r=1 (s=0) are inserted
        analytically, so the constant/contact and zero endpoints are retained.
        """
        if degree < 1 or tau <= 0:
            raise ValueError("degree>=1 and tau>0 are required")
        x = np.cos(np.pi*np.arange(degree+1)/degree)
        r = (x+1)/2
        values = np.empty_like(r)
        values[0] = 0.0
        values[-1] = float(self.v(np.array([2*self.s0]), mu, zeta)[0])
        interior = np.arange(1, degree)
        s = tau*(1-r[interior])/r[interior]
        values[interior] = self.v(s, mu, zeta)
        coefficients = dct(values, type=1)/degree
        coefficients[[0, -1]] *= 0.5
        return coefficients


def bernstein_scalar(nodes: np.ndarray, r: np.ndarray) -> np.ndarray:
    """Stable de Casteljau evaluation, used for scalar error diagnostics."""
    r = np.atleast_1d(r)
    v = np.broadcast_to(np.asarray(nodes)[:,None], (len(nodes),len(r))).copy()
    for size in range(len(nodes)-1, 0, -1):
        v = (1-r)*v[:size] + r*v[1:size+1]
    return v[0]


def apply_resolvent_polynomial(T, f, tau: float, minf: float, h_nodes):
    """Production recurrence: sparse elliptic solves only; no eigensolver/FFT.

    Computes (minf*I + R*B_N[h](R)) f, R=tau*(tau+T)^(-1).
    The caller owns the spatial-discretization error; this routine does not
    silently certify a finite grid as the continuum.
    """
    T = csc_matrix(T)
    f = np.asarray(f, dtype=float)
    if T.shape != (len(f),len(f)) or tau <= 0:
        raise ValueError("Inconsistent grid operator/input, or invalid tau")
    solve = splu(T+tau*eye(len(f), format="csc"))
    R = lambda x: tau*solve.solve(x)
    v = f[:,None]*np.asarray(h_nodes)[None,:]
    for size in range(len(h_nodes)-1, 0, -1):
        v = v[:,:size] + R(v[:,1:size+1]-v[:,:size])
    return minf*f + R(v[:,0])


def apply_chebyshev_resolvent(T, f, tau: float, coefficients):
    """Apply a Chebyshev polynomial of R=tau*(tau+T)^(-1).

    Clenshaw recurrence keeps only three input-shaped arrays and uses sparse
    elliptic solves. ``coefficients`` follow NumPy's Chebyshev convention for
    x=2R-I. No eigensolver or transform is used in this application path.
    """
    T = csc_matrix(T)
    f = np.asarray(f, dtype=float)
    coefficients = np.asarray(coefficients, dtype=float)
    if T.shape[0] != T.shape[1] or T.shape[0] != f.shape[0] or tau <= 0:
        raise ValueError("Inconsistent grid operator/input, or invalid tau")
    if coefficients.ndim != 1 or len(coefficients) < 2 or not np.isfinite(coefficients).all():
        raise ValueError("At least two finite Chebyshev coefficients are required")
    solve = splu(T+tau*eye(T.shape[0], format="csc"))
    def X(value):
        return 2*tau*solve.solve(value)-value
    b1 = np.zeros_like(f)
    b2 = np.zeros_like(f)
    for coefficient in coefficients[:0:-1]:
        b0 = 2*X(b1)-b2+coefficient*f
        b2, b1 = b1, b0
    return X(b1)-b2+coefficients[0]*f
