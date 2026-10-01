# Verification map

Run `python run_all_checks.py --list` to obtain the authoritative suite IDs.
The JSON registry is `manifests/verification_suite.json`; changing a filename
requires updating that registry and rerunning the suite.

## Scientific and bookkeeping layers

| ID | Script | Tested content | Interpretation |
|---|---|---|---|
| `one_loop` | [kt_kcss_one_loop_validator.py](../validators/kt_kcss_one_loop_validator.py) | Plus cumulants, additivity, diagonal normalization | Fixed-order bin regression |
| `toy_x` | [kt_kcss_toypdf_xconv_validator.py](../validators/kt_kcss_toypdf_xconv_validator.py) | Longitudinal plus prescription, leg interchange, flavor sums | Analytic toy PDFs |
| `generator` | [kt_kcss_generator_convolution_validator.py](../validators/kt_kcss_generator_convolution_validator.py) | Fourier polynomial dictionary and convolution table | Symbolic distribution algebra |
| `fixed_evolution` | [kt_kcss_one_loop_evolution_validator.py](../validators/kt_kcss_one_loop_evolution_validator.py) | Endpoint, path, and semigroup identities | Fixed-coupling one-loop model |
| `two_loop` | [kt_kcss_two_loop_ingredients_validator.py](../validators/kt_kcss_two_loop_ingredients_validator.py) | CS coefficients, projection, scalar matching logs | Two-loop structural identities |
| `running_expansion` | [kt_kcss_running_coupling_evolution_validator.py](../validators/kt_kcss_running_coupling_evolution_validator.py) | Integrated CS identity and endpoint expansion | Running-coupling expansion through second order |
| `w_expansion` | [kt_kcss_binned_wterm_prototype_validator.py](../validators/kt_kcss_binned_wterm_prototype_validator.py) | Constant auxiliary-scale cancellation and bin assembly | Fixed-order W prototype |
| `counting` | [kt_kcss_accuracy_ingredient_validator.py](../validators/kt_kcss_accuracy_ingredient_validator.py) | Loop-order manifest | Bookkeeping |
| `mock_import` | [kt_kcss_n3lo_matching_import_validator.py](../validators/kt_kcss_n3lo_matching_import_validator.py) | Flavor coverage, degree bounds, import contract | Synthetic coefficient records |
| `structural` | [validate_v1p4.py](../validators/validate_v1p4.py) | Contacts, canonical disk shifts, ordered matrix recurrence, counterexamples | Audited structural calculations |
| `resolvent_demo` | [resolvent_operator_demo.py](../validators/resolvent_operator_demo.py) | Bernstein operator bound against matrix exponential | Analytic non-QCD test multiplier, same grid |
| `unexpanded` | [validate_unexpanded.py](../validators/validate_unexpanded.py) | RG-integral identities, log-degree fixtures, contact-aware resolvent | Lower-order QCD and explicitly scoped finite-grid tests |
| `physical_canonical` | [validate_physical_canonical.py](../validators/validate_physical_canonical.py) | Sourced nf=5 4/4/3/3-loop coefficients, RG equations, independent quadrature, non-Gaussian action | N³LL evolution with strict NLO boundaries; no N³LL′ spectrum |
| `canonical_derivation` | [audit_canonical_derivation.py](../validators/audit_canonical_derivation.py) | Exact logarithmic degrees through order eight, omitted-loop controls, matching/hard identities | Finite witnesses supporting the analytic argument; not an all-order proof by testing |

All scripts are in `validators/`. Their outputs are in separate subdirectories of
`outputs/latest/`, rather than one shared directory with ambiguous historical files.

The additional [canonical physical audit](CANONICAL_PHYSICAL_AUDIT.md) requires
external, hash-verified LHAPDF data. It tests all strict-NLO components, contacts
and input derivatives against a separately coded continuum reference at three
scales. It is a separate executable workflow, not one of the offline suites.

## Theory-to-test connections

| Construction | Main executable evidence | Essential scope |
|---|---|---|
| Plus/generator algebra | `generator`, `structural` | Formal polynomial coefficients and admitted generator strip |
| Contact-preserving weak profiles | `structural` | Pure delta, contact-plus-density first bins, weak actions, additivity |
| Canonical disk shifts through degree six | `structural` | Independent triangular projection and generator expansion |
| Ordered flavor-matrix recurrence | `structural` | Noncommuting rational matrix fixtures, not finite physical N³LO tables |
| RG-completed unexpanded CS kernel | `unexpanded`, `physical_canonical` | Retained-system identities plus the sourced physical N³LL evolution vector |
| Endpoint domain and bounded canonical leg | Analytic argument in `FORMALISM.md`; sampled diagnostics in `unexpanded` | Samples supplement the written domain argument |
| Logarithmic-degree argument | `unexpanded` formal fixtures through ninth order | Finite checks of the all-order recurrence, not induction by test count |
| Resolvent convergence | `resolvent_demo`, `unexpanded` | Operator approximation on the same specified momentum grid |
| Exact contact retention | `structural`, `unexpanded` | Contact coefficients are not inferred by integrating only a regular density |
| Physical N³LL evolution input | `physical_canonical`, `tests/test_css2_resolution.py` | 4/4/3/3-loop evolution only; boundaries are strict NLO |
| Physical N³LL′ spectrum | Separate boundary/observable program | Not claimed by the evolution-input tests |

## Run counts

Read current counts from `outputs/latest/run_all_checks_summary.json` and the
individual suite reports. Reports exposing only aggregate metrics have no invented
individual record count. The component suites overlap with parts of the structural
and unexpanded suites, so quote counts by category. The `counting` and `mock_import`
suites check bookkeeping and synthetic records; they do not validate physical
high-order coefficient values.

## Numerical criteria

Each detailed JSON/CSV stores the actual residual and its test-specific tolerance.
The structural weak quadratures use 55 decimal digits in mpmath. Double-precision
RG and same-grid tests use the tolerances declared beside their calculations.
The aggregate runner requires both a successful subprocess and a new JSON success
field, and checks individual `passed` records wherever they are provided.

An independent matrix exponential or diagonalization serves as the test oracle on
the same momentum grid. The production recurrence itself uses sparse elliptic
solves. Neither same-grid comparison estimates the finite-box or discretization
error relative to the continuum.

## Contact regression

The historical profile script lost the origin contact in its first-bin formula.
It is superseded by [validate_v1p4.py](../validators/validate_v1p4.py), which tests
`C(t)=c -> c delta²(k)` and retains the contact explicitly in numerical first bins.
The check acts directly on scientific cumulants and test functions.

## Extending verification

When adding a new comparison, identify its input coefficients, reference source,
precision, domain, and success criterion. Give it a distinct suite or test ID.
Preserve a failing negative control. Separate analytic assumptions from numerical
checks, and preserve a distinction between a numerical bound on a finite matrix
and a continuum estimate.
