# Corrected CSS2 diagnostic prescription

This branch resolves the scientific discrepancies identified in the September 30
reference audit. It does not reproduce a historical FK2608 cache or claim a fit to
experimental data. The original repository and source packet are unchanged.

## Sourced ingredients and corrections

`src/css2_coefficients.py` supplies massless SU(3) quark beta/cusp/noncusp/CS-boundary
vectors at 4/4/3/3 loops. The convention is a=alpha_s/(4*pi), Gamma=gamma_K/2,
gamma=gamma_j, D=-K; da/dln(mu)=-2 sum beta_n a^(n+2). The complete evolution
content is N3LL. Boundary assembly retains strict NLO hard and OPE insertions;
this is not an N3LL-prime calculation. The original demo payload stays available.
Sources, exact versions and external archive hashes are in
[the source manifest](../manifests/css2_sources.json). The four-loop cusp is
implemented from von Manteuffel--Panzer--Schabinger and checked independently
against the Henn--Korchemsky--Mistlberger ancillary expression for nf=3,4,5,6.
The cross-check caught and corrected a transcription error in a quartic color
invariant before the physical benchmark was accepted.

For L=ln(Q/mu_i), the fixed-coupling one-loop single-leg logarithm is
`a*(6*CF*L-4*CF*L^2)`. Therefore the pair's Sudakov exponent is
`a*(8*CF*L^2-12*CF*L)`. The audited source integrates half of this exponent.
The new path uses the CSS2 anomalous dimensions directly and includes the
canonical CS boundary d(a_i) through three loops.

With Lb=ln(mu^2*b^2/b0^2), ell=ln(zeta/mu^2), and the whole-function plus
Pqq=CF*[(1+z^2)/(1-z)]_+, the quark matching insertion is
`a*(Cqq_finite - 2*Lb*Pqq + CF*((3-2*ell)*Lb-Lb^2)*delta(1-z))`.
The audited general-scale implementation omitted `3*CF*Lb*delta(1-z)`.
The finite canonical constants and hard coefficient follow Collins--Rogers;
the general-scale logarithms are derived from their matching evolution equations.
Both virtuality and rapidity equations are tested, including the LO PDF
splitting contribution. PDF functions return f, after explicit conversion from xf.

Hard factors remain outside universal TMDs. `strict_nlo` returns separate Born,
hard, qq and qg terms for both legs. No hard-times-OPE or OPE-times-OPE insertion
is included. Evolution remains unexpanded.

## Defined physical test point

This first bounded diagnostic uses proton--proton xA=0.1, xB=0.2, Q=20 GeV,
photon charge weights, both orientations of all five quark flavors, NNPDF4.0
NNLO member 0, and no nuclear/experimental prefactor or Y term. Official
LHAPDF 6.5.6 performs interpolation. It is built outside the repository. PDF
failures and extrapolation requests raise errors. The files are versioned and
hashed in [the PDF manifest](../manifests/css2_pdf.json).

The evolution coupling runs at four loops with fixed nf=5, anchored at
alpha_s(91.1876 GeV)=0.1180024. It is distinguished from the PDF set's published
NNLO evolution and tabulated three-loop alpha_s. All matching/evolution scales
for this diagnostic lie between 5 and 20 GeV, above the PDF's bottom threshold.
No continuity through a heavy-flavor threshold is being tested here.

For s=b^2, the profile uses raw=b0*sqrt(1/s+1/bmax^2), bmax=1.5 GeV^-1,
floor=5 GeV, and softplus floor/cap widths 0.12 GeV, with cap Q. Its matching
coordinate is explicitly b_eval=b0/mu_profile(s), and zeta_i=mu_profile(s)^2.
Thus its initial matching is canonical. This differs from the audited reference's
1.3 GeV floor and power-16 matching coordinate. The corrected general-scale
matching API is tested but is not required by this canonical diagnostic.

This Q-dependent remapping is a process prescription composed with universal
CSS2 evolution. The RG/CS tests hold b_eval and its boundary fixed. Differentiating
Q while also moving the profile is not asserted to be universal TMD evolution.
NP CS is zero for this specified diagnostic; it is not inferred from a missing
historical cache setting.

The capped symbol has nonzero M(0) and finite M(infinity). The latter acts as an
exact identity/contact term. The rational fit enforces sum(c)=M(0)-M(infinity).
It does not inherit the uncapped M(0)=0 or a zero 1/s-tail coefficient. The existing
uncapped implementation remains separate; its previously failing all-zero input
case is repaired and tested.

## Independent validation and limits

The production native-kT path uses the radial finite-volume Laplacian and signed
shared-shift elliptic solves. It never transforms the fitted transverse input.
Scalar spectral samples build the fixed operator. The longitudinal PDF matching
cache is checked against direct values and doubled quadrature order.

The benchmark uses an independently integrated four-loop coupling and Sudakov
integral, separately assembled strict luminosity, and an analytic Gaussian input
transform followed by adaptive Bessel quadrature. Its PDF matching is independently
checked with adaptive integration split at the official PDF x-grid knots. The
b-space benchmark shares physical coefficients and the checked PDF cache with the
native path; it does not share the radial discretization or rational approximation.

Each input leg is theta*G(0.5)+(1-theta)*G(1.4), where
G(v)=exp(-k^2/v)/(pi*v), and theta=0.37. It is a smooth, rapidly decaying
controlled input, not a compact-support realization or a DNN parameterization.
Bins use 2*pi*q*dq and edges [0,0.5,1,2,3,4,6,8] GeV. The noiseless closure
checks numerical bias; it is not a data fit or evidence of identifiability.

Rational degree, spacing and box size are refined separately. Sampled symbol
errors are not certified continuum suprema. The benchmark requires bin and
gradient errors below 3e-4 relative to the peak reference bin, closure bias below
5e-4, box change below 1e-7 and degree change below 1e-5. These are diagnostic
acceptance thresholds, not experimental uncertainty estimates.

## Reproduction

Run all repository tests and the scientific suite normally. For this new path,
use the pinned compatible Python environment and the external source/runtime
locations created for the resolution run:

```bash
python -m pytest tests/test_css2_resolution.py -q
python scripts/check_css2_sources.py --sources ../sources --out ../evidence/cusp_crosscheck.json
python scripts/benchmark_css2_point.py --runtime ../runtime --out ../evidence/physical_point
```

`setup_lhapdf_local.py` reproduces the isolated macOS LHAPDF build from the pinned
manifest. It neither modifies system Python nor downloads anything during tests.
The run directory retains the failed initial box study and the final refinements,
source checks, logs, runtime hashes and a handoff. The original prepared-source
ledger documents the base import; its optional exact-import check is not an
identity claim for these deliberately changed files.

The next physical extension requires a declared low-scale/threshold prescription,
noncanonical boundary transport if the old bOPE profile is adopted, explicit NP-CS
input, and resolved experimental row/target/unit metadata. Historical-paper
reproduction remains a separate unresolved task. No performance or training claim
follows from this fixed-point validation.

## Observed physical limitation

For the declared test input the independent W-only bins [4,6] and [6,8] GeV
are negative (approximately -0.01798 and -0.02956 in the diagnostic luminosity
normalization). Both evaluators reproduce this. These bins reach qT/Q=0.4;
no Y term is included. Agreement therefore validates the numerical map for
this signed target, not positivity or the adequacy of a full measured spectrum.
A physical extension must examine the low-scale/profile choice, matching to the
nonsingular contribution, and selected experimental range. No assumption is made
that a missing contribution will automatically cure this behavior.
