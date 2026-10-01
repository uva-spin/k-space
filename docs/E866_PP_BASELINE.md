# E866 proton-only forward baseline

This is a new, explicitly specified massless proton-proton W-only calculation.
It removes the nuclear target from the first measured-bin study. It does not
reproduce the earlier proton-copper fit or claim complete massive-quark accuracy.

## Data and observable

The primary input is Jason Webb's thesis, arXiv:hep-ex/0301031v1, Section 5.2,
Eq. (5.1), Tables 5.3-5.6. The source supplies 279 hydrogen-target transverse
bins, not deuterium or nuclear-target measurements. The extraction retains the
original source line, all bin bounds, mean columns, values, separate statistical
and point-to-point systematic errors, and common 6.5% normalization uncertainty.
HEPData record 41970 version 1 contains only Tables 1-32 for mass/xF spectra;
it does not contain the transverse Table 40 referenced by a Rivet implementation.
The absent-table response is preserved as an error, not accepted as data.

The candidate core has 41 bins with Q_low>=6.2 GeV, xF in [-0.05,0.15] or
[0.15,0.35], and qT_high/Q_low<=0.2. The mass windows avoid the upsilon region.
This is a new declared selection. No claim is made that it reproduces a prior fit.

The published invariant is I=E d3sigma/d3q=dsigma/(pi dy dqT^2). For a
charge-summed transverse structure function W in GeV^-2,

```math
I(Q,x_F,q_T)=\frac{4\pi\alpha_{\rm em}^2}{9sQ^2}
             W(Q,x_F,q_T)\,(\hbar c)^2.
```

The forward baseline integrates dQ^2 over the mass window, averages uniformly
over xF and qT, and returns pb/GeV^2/nucleon. It does not divide by the mass
width, add a target factor, or multiply an invariant by an extra rapidity
Jacobian. The radial operator uses its proper 2*pi*k*dk measure; the final
measurement functional is a qT average of the invariant, not an annular yield.
An explicitly separate annular average is evaluated as a bin-convention
sensitivity, not silently interchanged with the default.

At leading power, xA*xB=Q^2/s and xA-xB=xF. The exact relation
xF=2*sqrt(Q^2+qT^2)*sinh(y)/sqrt(s) differs by power corrections. Those corrections
are not included in this baseline. Alpha_em=1/137.035999084 is fixed; only
photon exchange is included. sqrt(s)=38.8 GeV.

## A source-metadata limitation

253 of the 279 tabulated mean energies exceed the maximum possible c.m. pair
energy within their stated mass, xF and qT bounds, using the thesis definition
xF=2*pL/sqrt(s). For example, the first bin lists <E>=6.277 GeV whereas the bin
allows at most sqrt(5.2^2+0.25^2+(38.8*0.15/2)^2)=5.964109 GeV. Rounding cannot
explain this difference. The source values are retained without correction.
This check concerns the interpretation of the energy column; it does not prove
that the published invariant cross sections themselves are wrong.

The forward calculation uses the published observable equation and bin edges,
not these inconsistent mean energies. The source does not fully specify the
finite-bin conversion from its underlying cross-section estimator to the
reported invariant. A uniformly averaged invariant, an annular average and an
estimator with a mean-energy conversion need not coincide. This is a remaining
measurement-systematics question for a precision absolute extraction. It does
not prevent a declared forward-operator benchmark or conditional pseudodata
closure. No undocumented rescaling of the measured values is applied.

## Internal normalization audit

An independent transcription of the 184 mass-spectrum rows in Table 5.1
permits 18 completely covered coarse mass/xF comparisons. Undoing the scaling
uses 1000*S*dM*dxF/M_reference^3; undoing the transverse invariant uses
pi*sqrt(s)*I*q_reference*dq*dxF/E_reference. Both mean and center choices are
retained. This is an approximate internal comparison, not an independent
chi-square: the tables share events, their cross-covariance is absent, and
rebinnings, mean conversions and unmeasured transverse tails matter.

For xF=0.15-0.35, the 6.2-7.2 and 7.2-8.7 GeV windows give transverse-to-mass
integral ratios 0.9829 and 1.0189 when the listed energies and mean momenta are
used. Replacing only the energies by their kinematic values at the listed means
gives 1.4417 and 1.3914. Using centers for mass and qT in the first comparison
gives 0.9764 and 1.0125. Several other windows show the same pattern, although
low-xF and sparse high-mass windows do not close and remain inconclusive.

This is evidence that the energy inconsistency may enter the reported
normalization, beyond a mislabeled column. It does not establish an exact
correction. In particular, the energy of mean kinematics need not equal a mean
energy, and averaging and rebinning need not commute. No cross section, bin or
uncertainty is corrected on this basis. The declared forward calculation is
validated independently of this unresolved experimental interpretation.

## Perturbative and nonperturbative prescription

Use the pinned NNPDF40_nnlo_as_01180 member 0 through official LHAPDF, with f
rather than x*f. The massless evolution has nf=5, four-loop beta/cusp and
three-loop noncusp/CS boundary. Alpha_s(91.1876)=0.1180024 is evolved with that
four-loop beta function; the PDF set itself remains its published NNLO set.
Matching and hard coefficients use strict NLO insertions, with no hard*OPE or
OPE*OPE products and with unexpanded evolution.

This baseline chooses a small perturbative b_max=b0/(5 GeV), distinct from the
previous noncanonical diagnostic. With s_b=b^2,

```math
b_*=(1/s_b+1/b_{\max}^2)^{-1/2},\qquad
b_{\rm eval}=[b_*^{16}+(b_0/Q)^{16}]^{1/16},\qquad
\mu_i=b_0/b_{\rm eval}.
```

The matching logarithm Lb is exactly zero. Every matching/evolution scale in the
selected mass windows remains above the pinned bottom threshold of 4.92 GeV.
The implementation verifies the actual endpoint and rejects crossings. No
below-domain PDF extrapolation, nf switch, or massive-TMD matching is invoked.
This resolves the running-domain problem for this declared massless baseline;
it does not calculate omitted heavy-mass corrections. A conservative small
b_max shifts more physics into the boundary and CS model. It cannot be combined
with a nonperturbative factor taken unaltered from a different b_max prescription.
Collins--Rogers explicitly discuss this dependence of the g-functions.

Each proton leg shares a normalized Gaussian boundary with variance lambda0 at
sqrt(zeta)=Q0=5 GeV. The nonperturbative single-leg CS contribution is
D_NP(b)=g*b^2, with D=-K. The two-leg result is consequently

```math
\exp[-\lambda_0 b^2/2-2g b^2\ln(Q/Q_0)],\qquad
v_{\rm pair}(Q)=2\lambda_0+8g\ln(Q/Q_0).
```

lambda0 and g have units GeV^2. Native evaluation takes a Gaussian directly in
kT with this pair variance and propagates its analytic parameter derivatives
through the same shifted elliptic solves. No fitted transverse input is Fourier
transformed. The independent validation path uses the analytic Gaussian
transform and Bessel quadrature. The nominal validation input is lambda0=0.7,
g=0.05; these are registered test parameters, not values inferred from E866.

## Validation and identifiability

The benchmark keeps all six Born/hard/qq/qg components separate under mass and
xF integration. It checks independently integrated evolution, adaptive PDF
convolutions, interpolation, transverse integration, mass/xF quadrature, radius,
pole count, gradients and noiseless two-parameter closure. The other xF interval
is held out from closure fitting. The xF integral can be taken before the radial
operator only because this particular input is x/flavor independent.

A real identifiability restriction remains. If one additionally releases
lambda(x)=a+b*ln(1/x), the pair width at fixed beam energy is

```math
v=2a+b\ln(s/Q^2)+8g\ln(Q/Q_0).
```

Its Jacobian has rank two, not three. A change in b can be compensated exactly
by a change in a and g. A single-energy E866 study cannot separately extract
all three, even with perfect data or a more powerful optimizer. The conditional
closure fixes b=0. A broader extraction needs additional information, such as
a second collision energy, an independently constrained CS kernel, or a stated
restriction on the input's x dependence.

No Y term or positivity clipping is introduced. A matched full spectrum still
requires FO minus the expansion of this same W and separate validation.

## Actual data comparison

After numerical validation, a restricted Gaussian diagnostic fits lambda0, g
and one common normalization N. Statistical and point-to-point systematic errors
are added in quadrature and treated as diagonal; N has a 6.5% Gaussian prior.
This is an explicit likelihood approximation. It includes no PDF or scale
uncertainty and does not establish confidence-interval coverage.

The uniform-qT fit returns lambda0=0.302206 GeV^2, g=0.252158 GeV^2 and
N=1.216835, with chi2=119.045 for 41 points plus one prior and three parameters
(39 nominal residual degrees of freedom). The annular fit gives chi2=120.091
and N=1.222725. No bounds are active. Independent best-point forward checks
pass in both cases. Optimizer convergence is not physical adequacy: chi2/dof
is about 3.05, and the normalization shift is 3.34 times its prior width.
The restricted model and measurement assumptions do not adequately describe
these data; the fitted widths must not be presented as a reliable extraction.

Changing the qT averaging convention changes nominal predictions by at most
2.664% and does not repair the poor fit. The existing results do not isolate
whether a more flexible transverse input or an experimental conversion issue
dominates. Resolve or control the estimator before widening the model for a
precision absolute fit. A different-energy, non-nuclear measurement is also
needed to break the exact x/CS degeneracy described above if both are released.

Primary sources: [E866 thesis](https://arxiv.org/abs/hep-ex/0301031v1),
[Collins--Rogers](https://arxiv.org/abs/1705.07167v3).
The finite-bin interpretation is registered as a choice with sensitivity tests,
not attributed to a nonexistent primary prescription.
