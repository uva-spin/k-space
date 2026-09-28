# Formalism and computational construction

This guide explains the implemented distribution calculus and unexpanded operator,
including their hypotheses, conventions, defining equations, and executable checks.
The source lineage is recorded in `manifests/source_provenance.json`. Numerical
checks supplement the analytic construction under the assumptions stated below.

## 1. Parent scheme and conventions

The parent operator definition is Collins/CSS2. The transverse-momentum
representation changes the calculus used to evaluate that definition.

```math
a=\frac{\alpha_s}{4\pi},\quad
\frac{da}{d\ln\mu}=-2\sum_{n\ge0}\beta_n a^{n+2},\quad
D_i=-\widetilde K_i^{\rm Collins},
```

```math
\frac{\partial F_i}{\partial\ln\sqrt\zeta}=-\mathcal D_i\otimes_T F_i,
\qquad
\frac{d\mathcal D_i}{d\ln\mu}=2\Gamma_i\,\delta^{(2)},
```

```math
\gamma_F^i(\mu,\zeta)
=\gamma_i(a(\mu))-\Gamma_i(a(\mu))\ln\frac{\zeta}{\mu^2}.
```

The forward Fourier convention used in **proofs and coefficient dictionaries** is
`F_tilde(b) = integral d²k exp(+i b·k) F(k)`. The inverse carries `(2π)⁻²`.
The evaluator does not transform fitted arrays. The process hard function multiplies
the two-leg observable; it is outside each universal TMD definition.

## 2. Distribution calculus and origin contacts

Define conventional radial plus distributions with the logarithm
`ln(mu²/k_T²)`. For a Schwartz test function,

```math
\langle\mathcal L_n(\mu),\varphi\rangle
=\int\frac{d^2k}{\pi k_T^2}\ln^n\frac{\mu^2}{k_T^2}
\left[\varphi(k)-\varphi(0)\theta(\mu-k_T)\right].
```

Their finite-bin actions are cumulants. For a bin beginning at the origin,

```math
\mathcal B_n(0,K;\mu)
=-\frac{1}{n+1}\ln^{n+1}\frac{\mu^2}{K^2},\qquad
\mathcal B_\delta(0,K)=1.
```

For a positive-lower-radius annulus the action is the difference of disk cumulants.

The generator and derivative basis satisfy

```math
\mathcal F[\mathcal G_\eta]=e^{-\eta L_b},\quad
L_b=\ln\frac{\mu^2 b_T^2}{b_0^2},\quad b_0=2e^{-\gamma_E},
\qquad
\mathcal F[\mathfrak L_n]=L_b^{n+1}.
```

The direct generator convolution identity holds in the convergent strip
`Re η>0`, `Re ξ>0`, `Re(η+ξ)<1`; its Taylor coefficients give the fixed-order
polynomial algebra. Larger multiplier products are admitted individually.

```math
\mathfrak L_m\otimes_T\mathfrak L_n=\mathfrak L_{m+n+1},\qquad
\mathcal L_0\otimes_T\mathcal L_0=-2\mathcal L_1,
```

```math
\mathcal L_0\otimes_T\mathcal L_1
=-\frac32\mathcal L_2-2\zeta_3\delta^{(2)}.
```

### Weak scale setting

Set `t=k_T²`, let `ψφ(t)` be the angular average of `φ(k)`, and let `C_p(t)` be
the cumulative after inserting the declared scale profile. The contact-preserving
weak definition is

```math
\langle\mathcal S_p A,\varphi\rangle
=-\int_0^\infty C_p(t)\psi_\varphi'(t)\,dt.
```

It represents the radial lift of `d[θ(t) C_p(t)]/dt`. With the stated local
integrability/growth conditions, its disk action is `C_p(K²)`. If the origin
limit is finite and the cumulative has the required absolute continuity,

```math
\mathcal S_p A
=C_p(0^+)\delta^{(2)}(k)+\frac1\pi C_p'(k_T^2).
```

A constant cumulative is therefore a delta distribution, not zero. The structural
validator includes this exact negative control, contact-plus-density first-bin
quadratures, logarithmically singular origin cases, and bin additivity.

## 3. Fixed-order boundaries and flavor ordering

The canonical disk moments of the derivative basis through degree six are

```math
(J_0,\ldots,J_6)
=(1,0,0,-4\zeta_3,0,-48\zeta_5,160\zeta_3^2).
```

They convert fixed-order parent polynomials into momentum-disk boundary constants.
These constants have a specific role: they are required by the momentum-disk
boundary convention. The **spectral-canonical construction below uses the parent
boundary at L_b=0 instead**. Its later projection generates the disk contacts;
adding them to that parent boundary a second time would double count them.

Matching logarithms satisfy ordered flavor-matrix equations. Anomalous dimensions
act on the outgoing TMD row, while collinear splitting kernels act on the right.
In the stated DGLAP convention,

```math
P_\mu=4aP^{(0)}+8a^2P^{(1)}+16a^3P^{(2)}+\cdots.
```

The structural suite checks the complete ordered recurrence through third order
with noncommuting rational matrix fixtures. It also checks every partition
`r+s+t=n` of the hard/two-leg expansion through `n=3`. Physical finite matching
coefficients remain independent inputs, separate from this logarithmic recursion.

## 4. Explicit spectral-canonical completion

Let `T=-Δ_k` be the nonnegative self-adjoint momentum-space Laplacian. Its scalar
spectral variable `s` has units of inverse momentum squared. It represents the
same transverse resolution as `b_T²` in the analytic proof.

For fixed `mu0>0`, set `s0=b0²/mu0²` and `y=s/s0`. Define the smooth switch

```math
\chi(u)=\begin{cases}
0,&u\le0,\\
\dfrac{e^{-1/u}}{e^{-1/u}+e^{-1/(1-u)}},&0<u<1,\\
1,&u\ge1.
\end{cases}
```

The implementation uses the numerically stable logistic equivalent. The scale is

```math
t_c(s)=-\frac12[1-\chi(2y-1)]\ln y,
\qquad \mu_c(s)=\mu_0 e^{t_c(s)}.
```

It is canonical for `s≤s0/2`, constant for `s≥s0`, and independent of hadron,
longitudinal momentum, hard process, and evolution endpoint. There is no ultraviolet
cap. The chosen beta polynomial is required to remain on its asymptotically free
branch throughout the integrated interval.

## 5. Unexpanded CS and endpoint evolution

The CS symbol is first reconstructed by an RG integral:

```math
D_i^c(s;\mu)
=d_{i,[3]}\!\left(a(\mu_c(s))\right)
+2\int_{\mu_c(s)}^\mu\frac{d\nu}{\nu}\Gamma_{i,[4]}(a(\nu)).
```

For endpoints `e_a=(mu_a,zeta_a)` and `e_b=(mu_b,zeta_b)`,

```math
R_{i;ba}=\int_{\mu_a}^{\mu_b}\frac{d\nu}{\nu}
\left[\gamma_{i,[3]}(a(\nu))-
\Gamma_{i,[4]}(a(\nu))\ln\frac{\zeta_a}{\nu^2}\right],
\qquad \rho_{ba}=\ln\sqrt{\frac{\zeta_b}{\zeta_a}},
```

```math
u_{i;ba}(s)=\exp\left[R_{i;ba}-\rho_{ba}D_i^c(s;\mu_b)\right],
\qquad \mathbb U_i(e_b;e_a)=u_{i;ba}(T).
```

This solves the retained-system integrability and endpoint composition equations.
The endpoint symbol may be unbounded at spectral zero, with at most logarithmic
power growth. A common domain is

```math
\mathscr D_{\log}=\bigcap_{m\ge0}\mathrm{Dom}\,w(T)^m,
\qquad w(s)=1+\ln(1+\tau/s).
```

The analytic argument in the source construction proves density and invariance
on this domain. The finite endpoint tests exercise identities on specified samples;
they do not determine the functional-analytic domain by testing it numerically.

### Complete canonical leg

The canonical leg is

```math
v_i(s;\mu,\zeta)=\exp\left[
-d_i(a_c)\ln\frac{\sqrt\zeta}{\mu_c}
+\int_{\mu_c}^{\mu}\frac{d\nu}{\nu}
\left(\gamma_i(a(\nu))-
\Gamma_i(a(\nu))\ln\frac{\zeta}{\nu^2}\right)\right].
```

It has a bounded completion, tends to zero at spectral zero, and is constant
beyond the infrared completion. It satisfies `u_ba=v_b/v_a`. With admitted
canonical matching/PDF multipliers attached, this gives the bounded operator used
for extraction. The longitudinal assumptions include continuity and controlled
logarithmic growth at fixed x or admitted Mellin moments; endpoint uniformity at
x=0 or x=1 is a separate problem.

## 6. Momentum-resolvent evaluation

For `tau>0` with units of inverse momentum squared,

```math
R_\tau=\tau(\tau+T)^{-1},\quad
r_\tau(k)=\frac{\tau}{2\pi}K_0(\sqrt\tau\,k_T),\quad
\int d^2k\,r_\tau(k)=1.
```

For a completed canonical symbol `M`, write

```math
r=\frac{\tau}{\tau+s},\quad
m(r)=M\!\left(\tau\frac{1-r}{r}\right),\quad
h(r)=\frac{m(r)-M_\infty}{r}.
```

The continuous endpoint values are supplied analytically. Bernstein polynomials
of `h` define

```math
\mathbb M_N=M_\infty I+R_\tau B_N[h](R_\tau),\qquad
\epsilon_N=\|B_N[h]-h\|_{\infty}\to0.
```

The momentum-space implementation uses operator-valued de Casteljau recurrence and
sparse solves with `tau*I+T`. It preserves the constant contact coefficient exactly.
The source theorem gives

```math
\|(\mathbb M_N-M(T))f\|_2\le\epsilon_N\|f\|_2,
\qquad
\|(\mathbb M_N-M(T))\delta^{(2)}\|_2
\le\sqrt{\frac{\tau}{4\pi}}\epsilon_N,
```

where the second norm is the regular error after the identical contact terms
cancel. A finite-bin integral is bounded by its area's square root times the
regular L² error. The continuum bound uses the true uniform symbol error.
The shipped grid experiment computes a bound for its finite matrix, not an
estimate of the continuum supremum or grid truncation error.

## 7. Logarithmic counting

At final scales `(Q,Q²)`, set `L=ln(Q/mu_c)` and `a_-(t)=a(Q exp(-t))`. The
retained exponent is

```math
E_i(L,a_Q)=\int_0^Ldt\,[\gamma_{i,[3]}(a_-(t))
-2t\Gamma_{i,[4]}(a_-(t))]
-Ld_{i,[3]}(a_-(L)).
```

The source argument uses the running-coupling recurrence at arbitrary order. First
omitted cusp, noncusp, rapidity-boundary, and beta coefficients change only

```math
\Delta E_i\in\sum_{n\ge4}a_Q^n\,\mathbb C[L]_{\le n-3}.
```

The retained exponent towers have degrees `n+1`, `n`, `n−1`, and `n−2`. This is
the declared N³LL evolution convention. Physical third-order hard and canonical
matching boundaries supply N³LL′. The tests examine the recurrence through order
nine with explicit rational fixtures and reject an `a⁴L⁵` counterexample; the
all-order induction is the analytic support, rather than the finite test count.

Uniform resummed-remainder statements additionally require the weak-coupling
remainder/derivative assumptions stated in the source argument and a compact
resummation window away from the RG singularity.

## 8. Nonperturbative modules and observables

A sufficient module class is a localized radial L¹∩L² smearing with controlled
moments. The operator acts on it without selecting a neural architecture.
Admissibility of an extra nonperturbative CS contribution requires the declared
product domain and the same representation/universality conditions.

The process hard factor acts on the two-leg observable. The universal TMD tail
comes from matching coefficients; the observable-level remainder is
`Y = fixed order − expanded perturbative singular W`. The statistical fit,
full-spectrum matching, heavy-flavor thresholds, and tensor harmonics each carry
their own inputs and tests.

## Code entry points

- `src/resummed_operator.py`: `Coefficients`, `CanonicalRG`,
  `bernstein_scalar`, `apply_resolvent_polynomial`.
- `validators/validate_unexpanded.py`: symbolic recurrence fixtures, endpoint
  tests, and the independent finite-matrix oracle.
- `validators/validate_v1p4.py`: contact-complete distribution calculus,
  finite boundary conversions, matrix ordering, and negative controls.
- `validators/resolvent_operator_demo.py`: analytic bounded test-multiplier demo.

See [VERIFICATION.md](VERIFICATION.md) for the complete test map and
[REFERENCES.md](REFERENCES.md) for the imported theoretical sources.
