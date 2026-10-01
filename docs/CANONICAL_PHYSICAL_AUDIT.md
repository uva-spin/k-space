# Canonical physical assembly and independent validation

The canonical physical benchmark connects the uncapped spectral-canonical
evolution prescription directly to actual PDFs and all six strict-NLO
Born/hard/qq/qg components. It does not use the capped E866 process profile.
Its final endpoints are Q=8, 20, and 91.1876 GeV, with distinct xA/xB values and
both Gaussian and algebraic-tail Student inputs on two legs.

## Independent calculations

Production uses the existing sourced four-loop/three-loop CSS2 RG system,
canonical matching from `css2_matching`, and a cache of longitudinal blocks.
Native two-leg Student convolution is a Feynman-parameter integral in momentum
space. A shared Chebyshev recurrence applies sparse elliptic resolvents and
stores only bin moments. Component predictions are contractions of those moments
with the separately calculated scalar symbol coefficients.

The reference imports none of the production RG, matching, physical-symbol or
flavor-assembly routines. It shares the numerical physical coefficients and the
official PDF provider. It uses endpoint-anchored RK45 integration of the full
anomalous dimension, canonical matching in log(xi) rather than z, an independent
charge/antiparticle loop, and adaptive Bessel integration of a product of analytic
input transforms. Direct adaptive longitudinal calculations check its cache.
The logarithmic input-scale derivative and a separate delta/contact input are
tested along with ordinary bins.

The derivation validator uses exact arithmetic through order eight and deliberate
missing-loop and incorrect-product controls. Finite witnesses support the
all-order degree argument; they do not replace it.

## Scope and conditions

- Accuracy means N3LL evolution with strict-NLO hard and OPE boundaries.
- The benchmark is a massless fixed-nf5 photon W-term calculation, including a
  formal five-flavor ultraviolet continuation. It does not implement physical
  top-threshold matching, the Z contribution, or a Y term.
- The canonical matching coordinate is b_eval=b0/mu_c. It equals b in the
  ultraviolet canonical region and follows the declared completion elsewhere.
- The PDF table ends at 1e5 GeV; production never extrapolates it or caps mu_c.
  The first reference omits b<b0/1e5 and bounds it for a specified longitudinal
  continuation class, |f(x,mu)| <= |f(x,1e5)|[1+ln(mu/1e5)]^4, including the
  convolution arguments. The bound uses monotone RG decay on successive
  intervals. Its numerical evaluation is not an interval-arithmetic certificate.
  A separate completion auditor then integrates this interval for two explicitly
  defined benchmark inputs: endpoint PDFs constant above the table, or multiplied
  by [1+ln(mu/1e5)]^4. Neither is asserted to solve DGLAP above the table. Coupling
  and RG evolution remain uncapped. The broader class bound is retained in the
  bin, component, gradient, contact and observed-budget tests.
- Signed bins out to 80 GeV stress the numerical operator. The high-qT W term
  is not interpreted as the complete physical spectrum. Positivity is not imposed.
- Reported errors are normalized to each case's peak reference total. Per-bin
  relative errors are retained separately. An observed refinement budget is not
  a rigorous global operator-norm certificate.

## Reproduction

```sh
python -m pytest -q tests/test_canonical_physical.py
python validators/audit_canonical_derivation.py --out outputs/derivation
python scripts/benchmark_canonical_physical.py \
  --library /path/to/libkt_lhapdf.dylib \
  --pdf-data /path/to/LHAPDF \
  --extended --out outputs/canonical_physical
python scripts/audit_canonical_uv_completion.py \
  --input outputs/canonical_physical/canonical_physical_audit.json \
  --library /path/to/libkt_lhapdf.dylib \
  --pdf-data /path/to/LHAPDF --out outputs/completed_physical
python scripts/audit_canonical_support.py \
  --physical-report outputs/completed_physical/completed_physical_audit.json \
  --out outputs/canonical_support
```

The extended finite-table command returns exit code 1 because its auxiliary
1e-6 target for the entire unspecified ultraviolet class fails at the high
scale (1.67e-6). This failed report must be preserved. The completion auditor
rejects any other baseline failure and tests both specified inputs with the same
1e-6 target for unresolved integration, now below 1e-17. Their maximum difference,
1.54e-6 of the peak, is input dependence. It also requires the broader class-wide
errors to meet the original bin/component/contact limits of 5e-4, gradient limit
of 1e-3, and observed budget limit of 1e-3. No primary tolerance was relaxed.
Run these commands sequentially; a wrapper may handle the documented baseline
exit 1 only if it then requires a successful completion and support audit.

The first high-scale run at radius 1280 GeV failed its bin, component and
refinement-budget tests. The extended configuration enlarges the high-scale
domain to 5120 GeV and separately refines degree, spacing and radius. The
acceptance thresholds stay fixed. These large radii are numerical domains,
not physical kinematic claims or recommended future fit settings.

The benchmark and the separately coded derivation audit were performed within
one development effort. They must not be described as independent human peer
review or as validation of an unrestricted phenomenological fit.
