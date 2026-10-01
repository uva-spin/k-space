# Canonical CSS2 derivation and implementation correspondence

This document records the analytic reasoning behind the canonical method and
its separately coded checks. It does not constitute independent human peer
review. The imported source is
Collins--Rogers, arXiv:1705.07167v3, especially Eqs. (7)--(11), (58)--(62),
(74)--(76); the archived TeX and hashes are in the source manifest.

## 1. Signs, two scales, and path independence

With `a=alpha_s/(4 pi)`, `D=-K`, `Gamma=gamma_K/2`, and `gamma=gamma_j`,

    d a / d ln(mu) = -2 B(a), B(a)=sum beta_n a^(n+2)
    d ln(F) / d ln(mu) = gamma(a)-Gamma(a) ln(zeta/mu^2)
    d ln(F) / d ln(sqrt(zeta)) = -D
    d D / d ln(mu) = 2 Gamma(a).

Both mixed derivatives equal `-2 Gamma`. At fixed spectral coordinate and fixed
completion, evolving first in mu and then in rapidity gives

    ln U(f,i) = int_mu_i^mu_f [gamma-Gamma ln(zeta_i/nu^2)] dln(nu)
                - ln(sqrt(zeta_f/zeta_i)) D(s;mu_f).

The alternative order agrees after inserting the integral for D. Composition
then follows by adding path integrals, with no Fourier-space approximation.
Changing a Q-dependent profile during this differentiation would add chain-rule
terms and is excluded from the canonical claim.

## 2. Canonical boundary and ultraviolet endpoint

Set `s=b^2`, `s0=b0^2/mu0^2`. The smooth scale equals `b0/sqrt(s)` for
`s<=s0/2` and `mu0` for `s>=s0`, independently of final endpoints. Define

    D(s;mu)=d(a_c)+2 int_mu_c^mu Gamma(a) dln(nu).

Evolution from `(mu_c,mu_c^2)` to `(mu,zeta)` gives exactly the implemented
`log_v_tc`. A fresh endpoint-anchored RK45 reference integrates
`gamma-2 Gamma ln(Q/nu)` directly; it does not reuse G,N,J.

Asymptotic freedom gives `a~1/(2 beta0 ln(nu/Lambda))`. Consequently,
`ln v = -(Gamma0/beta0) ln(mu_c/Lambda) + O((ln ln mu_c)^2)`.
The cusp sign is negative for evolution down from arbitrarily large mu_c.
For any epsilon>0 below p=Gamma0/(2 beta0), `v=O(s^(p-epsilon))`.
Thus `v(0)=0`, not a finite scale assignment. For s>=s0, v is constant.

The implementation assigns the exact infrared constant at `tc(infinity)` before
taking a logarithm and rejects NaN inputs. This avoids an indeterminate endpoint
expression while preserving the finite-s prescription.

## 3. Endpoint domain and contacts

For `w=1+ln(1+tau/s)`, the dense domain `intersection_m Dom(w(T)^m)` contains
Schwartz inputs. Finite endpoint ratios grow at most as a finite power of w;
multiplication preserves this common domain. The bounded v maps L2 into it.
These statements concern T=-Delta_k, whose spectral coordinate has units GeV^-2;
s is not measured k_T^2.

If `M-M_infinity` is bounded and compactly supported in s, its inverse Fourier
image is L2. Thus `M(T) delta = M_infinity delta + g` with g in L2. This is a
distributional extension: delta is not incorrectly placed in L2. With the
declared Fourier pair, `||g_error||_2^2=(1/(4 pi)) int |E(s)|^2 ds`.
A weight `(1+s/tau)` gives the stated `tau/(4 pi)` bound. Finite-bin bounds
follow from Cauchy--Schwarz, including an explicit origin convention.

The Bernstein construction proves existence and weighted contact convergence.
The new shared Chebyshev recurrence is a practical polynomial evaluator checked
against eigendecomposition and the separately organized Clenshaw recurrence.
Its observed continuum refinement must not be promoted to a rigorous global
supremum certificate or a Bernstein rate.

## 4. Two ordinary legs

For L2 legs, convolution is bounded: `||f*g||_infinity <= ||f||_2 ||g||_2`.
If e_A and e_B are their errors, the pointwise convolution error is at most
`||e_A||_2 ||f_B||_2 + ||f_A||_2 ||e_B||_2 + ||e_A||_2 ||e_B||_2`.
Multiplying by the finite bin area gives its bin bound. This does not assert
that arbitrary L2 convolutions are integrable over the whole plane.

For the new benchmark the two normalized transverse boundary legs are Student
densities (p=3, q=4) or Gaussians. The Student convolution is evaluated natively
by a Feynman-parameter integral over the two-dimensional momentum convolution.
Its reference is the product of two analytic Bessel-K transforms. This tests
an actual two-leg, algebraic-tail input, not a relabeled Gaussian closure.

## 5. Retained logarithmic towers and remainder

At `(mu,zeta)=(Q,Q^2)`, in the canonical perturbative region with mu_c<=Q,
`L=ln(Q/mu_c)` and `a_-(t)=a(Q exp(-t))` give

    E = int_0^L [gamma(a_-)-2t Gamma(a_-)]dt - L d(a_-(L)).

Induction on the running-coupling ODE gives degree(A_n)<=n-1 in
`a_-=sum a_Q^n A_n(t)`. First omitted cusp a^5 contributes degree<=n-3
after multiplication by t and integration. First omitted gamma or d at a^4
also contributes degree<=n-3. First omitted beta*a^6 changes a by a_Q^6*t,
again reaching degree<=n-3 in the cusp integral. Products preserve these bounds.
The four retained towers are therefore n+1,n,n-1,n-2: N3LL **evolution**.

For uniform remainder estimates one additionally needs local derivative and
remainder bounds, an interval below the first running singularity, and
`a_-(t)<=C a_Q`. Gronwall gives Delta a=O(a_Q^5), and the exponent error is
O(a_Q^3) at fixed a_Q L. These hypotheses do not prove convergence of QCD's
perturbation series. Nor is this a uniform statement about arbitrary infrared
completions, zeros of physical channels, or high-qT inclusive observables.

## 6. Matching and hard-factor normalization

Collins--Rogers (74) gives canonical cqq=CF[-pi^2 delta/6+2(1-z)] and
cqg=2z(1-z), with a=alpha_s/(4pi). In production, 4TF=2 in cqg.
The whole-function plus convention for Pqq is consistent with
`df/dln(mu)=4a P*f`. Differentiating the one-loop general-scale insertion
gives `-4a P*f + a CF(6-4 ell)f`; adding PDF running gives the required
TMD anomalous dimension. The rapidity derivative gives `-4a CF Lb f=-D1 f`.
The derivative of the hard insertion cancels the two TMD virtuality anomalous
dimensions at zeta_A*zeta_B=Q^4 to this order.

The physical benchmark evaluates canonical matching at
`b_eval=b0/mu_c(s)`, so Lb=ell=0 at the boundary. In the small-s canonical
region this is the actual b. In the infrared it is the declared completion;
nonperturbative physics is not inferred from this choice.
Six Born/hard/qq/qg terms are kept. Hard-times-OPE and OPE-times-OPE are
excluded by strict order counting. Antiparticle reversal, both q/antiquark beam
orientations and electric charge squares are independently assembled.

## 7. Finite PDF table versus an uncapped scale

The pinned NNPDF set has five flavors and Qmax=1e5 GeV. The massless nf5
benchmark uses a separate four-loop coupling; it is not a massive treatment
of crossing the physical top mass. Production raises if a spectral node would
request an untabulated PDF; it does not cap mu_c or extrapolate secretly.

The continuum reference integrates down to bcut=b0/Qmax and bounds the missing
full-symbol contribution. For the explicit sufficient class
`|f(x,mu)| <= |f(x,Qmax)| (1+ln(mu/Qmax))^4`, including convolution arguments,
the two-leg luminosity grows by at most the eighth power. Signed matching is
bounded using absolute integrals, not a canceled value of cqq.

For mu_c>=Qmax>Q, monotone decay of v follows by proving on the whole interval
`0<a<=a(Qmax)` that `(d-gamma)/a<=0` and `(Gamma+B*d')/a>0`. These are
finite polynomials; extrema are evaluated at endpoints and all stationary
roots. With `|q J1(qb)|<=q^2 b/2` and normalized positive input transforms
bounded by one, the omitted bin is bounded component by component by

    Delta(q^2)/2 * bcut^2 * v(Qmax)^2 * A_component
      * integral_0^infinity exp(-2u) (1+u)^8 du.

This is a conditional bound on a declared continuation class, not evidence that
the published finite PDF table specifies the unique infinite-scale continuation.
Its measured size is reported against the bin peak and the numerical budget.

The actual bound sharpens this constant-v expression by dividing u=ln(mu/Qmax)
into intervals on [0,24], bounding v^2 by its left endpoint on each interval,
and bounding the final infinite remainder by the last value. Polynomial moments
are integrated with incomplete gamma functions. A numerical ODE safety margin
is included; this is not interval arithmetic. The resulting high-scale class
bound is 1.6683e-6 of the peak; a universal 1e-6 bound is not asserted.

Two explicitly completed input fixtures are then integrated over the missing
interval with the uncapped running/evolution: endpoint PDFs constant above Qmax,
or multiplied by (1+u)^4. Remaining quadrature/tail estimates are below 1e-17
of the peak and the two inputs differ by at most 1.532e-6. They are sensitivity
fixtures, not DGLAP solutions. The broader class bound stays in the original
primary bin, component, derivative, contact and total-budget criteria, all of
which pass. This distinguishes the achieved class bound from the much smaller
integration error for a specified completed input.

The derivative bounds have an independent input-level justification. For the
p=3 Student leg, the absolute log-scale derivative is at most four times the
positive density, so its L1 norm is at most four. The Gaussian derivative also
has L1 norm at most four. Multiplying by the other normalized leg's transform
therefore gives the conservative factor four used in the class-wide derivative
bound. For the much smaller explicit-fixture residual, both transforms admit
a positive Gaussian scale-mixture representation: their log-scale derivative
is an average of -2w exp(-w), whose magnitude is at most 2/e < 1. Thus the
ordinary-input UV tail bound also bounds that derivative remainder.

## 8. Signed outputs are required, not a failed positivity test

The full admitted canonical symbol also tends to zero at s=0 when the
longitudinal boundary has at most polynomial logarithmic growth. If its output
is integrable, its total transverse integral is therefore zero. A nonzero
output cannot then be nonnegative everywhere. If integrability is absent, the
finite-positive-measure version of the obstruction still applies. This is not
the inclusive Drell--Yan cross section: the W term needs separate high-qT/Y
matching. The audit retains signed bins out to 80 GeV as an operator stress
test and makes physical low-qT statements only in the stated region.

## Relation to numerical evidence

The reconstruction supports the stated canonical evolution theorem under its
declared hypotheses. Full numerical acceptance is determined by the new audit
[published reports](VALIDATION_RESULTS.md), not by this prose. External scientific
review remains distinct from this technical audit.
