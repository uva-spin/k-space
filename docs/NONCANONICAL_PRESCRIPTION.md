# Noncanonical prescription and measurement readiness

This extension implements the power-16 matching coordinate with the corrected
Collins/CSS2 matching logarithms and a three-loop general-scale CS boundary. It
also pins the 54 selected E605 continuum points to their primary HEPData tables.
It does not yet define a full threshold-crossing experimental prediction.

## General-scale CS boundary

For a=alpha_s/(4*pi), L=ln(mu^2*b_eval^2/b0^2), D=-K, d1=0, the boundary is

```math
D(a,L)=a\Gamma_0L
+a^2\left(\frac{\beta_0\Gamma_0}{2}L^2+\Gamma_1L+d_2\right)
+a^3\left[\frac{\beta_0^2\Gamma_0}{3}L^3+
\left(\frac{\beta_1\Gamma_0}{2}+\beta_0\Gamma_1\right)L^2+
(\Gamma_2+2\beta_0d_2)L+d_3\right].
```

The coefficients follow by solving dD/dln(mu)=2*Gamma with the declared beta
function and canonical constants. A symbolic recurrence checks the result
independently. At a fixed b_eval and anchor scale, `BoundaryEvolution` transports
this boundary using the unexpanded RG integrals. Both scale derivatives and
arbitrary two-scale composition hold. Re-expanding at a different anchor is a
finite-order change; it must not be mistaken for exact transport of the same
boundary. Matching and evolution use the same b_eval.

The coordinate is `(bstar^16+(b0/(C5*Q))^16)^(1/16)`, with C5=1. Its Q dependence
remains an explicit process prescription. A moving-Q derivative of this map is
not asserted to be the universal fixed-coordinate CS equation. The old canonical
remapping remains available as a separate regression.

## Low-scale and threshold scope

The pinned NNPDF4.0 member has Qmin=1.65 GeV and a bottom threshold of 4.92 GeV.
The old floor=1.3 GeV falls outside its tabulated domain. No extrapolation or
replacement values have been introduced. The current implementation permits
only a declared fixed-flavor interval; it rejects a path across that threshold.
The complete soft floor/cap is evaluated before choosing the RG lower limit:
a soft cap can lie slightly below the nominal floor. The interval is extended to
the actual endpoint, preserving the profile, rather than clipping its values.

Two massless fixed-flavor diagnostics are independently checked:

| nf | Q (GeV) | Nominal floor (GeV) | Coupling anchor |
|---|---:|---:|---|
| 4 | 4.5 | 1.65 | Official PDF alpha_s(4.5), then four-loop nf=4 evolution |
| 5 | 7.5 | 5 | alpha_s(91.1876)=0.1180024, four-loop nf=5 evolution |

These use the real PDFs, power-16 coordinate, unexpanded evolution and strict
NLO insertions. Neither is a complete massive-quark calculation. The four-flavor
case explores the available low-scale PDF interval; the five-flavor case shares
the candidate E605 mass center but retains the higher floor. NP CS is explicitly
zero and the Gaussian-mixture input is a test input, not fitted physics.

Heavy-quark factorization depends on qT relative to the quark mass as well as Q.
Changing nf in beta and cusp coefficients is insufficient for a claim of complete
massive-TMD accuracy. The primary source 1703.09702v2, Sections 2.2–2.4 and 4,
contains collinear and soft mass-mode factors, their rapidity evolution, and
mass-dependent beam/soft functions. These contributions have not been imported
or converted to this CSS2 implementation. The threshold-crossing gate therefore
stays closed. No accuracy claim is inferred from the massless 4/4/3/3 vectors.

## E605 measurement contract

[The immutable point manifest](../manifests/e605_primary_points.json) maps all
54 selected rows uniquely to HEPData record 22831 version 1, Tables 17–21
(original publication Table 10). Values, mass windows, source hashes and error
labels are retained. The observable is E*d3sigma/dp^3 in pb/GeV^2/nucleon at
xF=-0.1 and the published qT points, for proton-copper collisions at sqrt(s)=38.8 GeV.
The primary data do not provide qT or xF bin edges. The 15% normalization and 5%
point-to-point systematics remain distinct from the quoted absolute errors.

The CSV's CS column is a derived A/PreFactor quantity. The adapter uses the
primary invariant A directly. It does not need that analysis-specific rescaling.
For every selected row, the stored y_Low/y_High equal the leading-power
fixed-xF rapidities at the two mass endpoints. They are not an additional
independent rapidity acceptance interval. A double integration over that range
at each Q would describe a different observable.

`InvariantPoint` provides a mass-window integral of the invariant prediction at
fixed xF and qT. Collins--Rogers eq:fact.CSS2 supplies the prefactor; converting
dqT^2 to d^2qT divides it by pi. With a charge-summed W in GeV^-2, the integrand
is `(4*pi*alpha_em^2/(9*s*Q^2))*W*dQ^2`, converted by hbar*c squared to pb/GeV^2.
It is not divided by the mass-window width. The exact pair geometry uses
`y=asinh(xF*sqrt(s)/(2*sqrt(Q^2+qT^2)))`; a separately named leading-power option
uses Q instead of the transverse mass. This geometrical option does not supply
missing power corrections to the TMD factorization theorem.

The adapter requires an explicit W callback, electromagnetic coupling and
kinematics choice. The callback must contain the selected target mixture and NP
input. No copper proton fraction, nuclear correction or covariance has been
invented. Full publication-level interpretation of the mass-window correction
and acceptance still needs confirmation against the original paper. The manifest
therefore rejects production prediction. A conservative candidate W-only core
qT/Q_low<=0.2 retains 52 of the 54 points; this is a declared new diagnostic cut,
not an assertion about the reference likelihood.

## Validation and matching scope

The noncanonical benchmark solves a separate coupled ODE for alpha_s and the
TMD logarithm and generates D by an independent coefficient recurrence. It
uses Bessel quadrature only for the validation path. Native evaluations use
radial elliptic solves; the point's longitudinal kernels are cached separately
and checked against direct official-PDF convolutions. Pole count, spacing, box,
input gradients and noiseless closure are refined. Tolerances are the preceding
resolution run's thresholds: bin/gradient errors below 3e-4 of the peak bin,
closure bias below 5e-4, box change below 1e-7 and degree change below 1e-5.

The nf=4 W-only diagnostic has a negative [3,4] GeV bin, which reaches qT/Q=0.89.
Numerical agreement does not validate that region as a physical cross section.
No positivity clipping or Y-term replacement is used.

The next matched extension must use Y=FO minus the expansion of the same defined
W prescription. Existing external reference grids cannot be added to the new
strict-NLO CSS2 result without order, scheme, target, observable and subtraction
checks. The reference's interpolation `(1-p)*W+p*FO` is separately labeled and
is not the conventional additive remainder. The execution plan places this full
spectrum extension after a stable measured W-only baseline.

Reproduction uses `scripts/benchmark_noncanonical.py`, the previous run's pinned
LHAPDF data/library and the new local bridge with alpha_s anchor access. Official
primary files and failed attempts remain in the external run directory.


## Pinned source and run records

[The source manifest](../manifests/noncanonical_sources.json) records exact paths,
versions and SHA-256 values. The inherited `css2_sources.json` describes the
previous resolution run; its files remain under the sibling
`m1_resolution_2026-09-30/sources`, not this run's new sources directory.
[The historical E605 readiness manifest](../manifests/e605_diagnostic_readiness.json)
preserves the unresolved copper-target diagnostic contract. The later E866 pp
demonstrator and canonical methods-paper validation are summarized in the
[current readiness manifest](../manifests/production_readiness.json).

Primary references: [Collins and Rogers](https://arxiv.org/abs/1705.07167v3),
[Pietrulewicz et al.](https://arxiv.org/abs/1703.09702v2), and
[E605 data version 1](https://www.hepdata.net/record/ins302822?version=1).
The external run's `reproduce.sh` supplies all runtime paths and commands.
