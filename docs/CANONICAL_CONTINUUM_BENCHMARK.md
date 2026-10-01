# Physical canonical continuum benchmark

`scripts/benchmark_canonical_continuum.py` compares the native radial evaluator
with an independent transform-based continuum reference for the identical
physical canonical symbol. The symbol uses fixed-nf=5 N³LL evolution. No hard or
OPE boundary enters this operator-only comparison.

## Non-Gaussian input and reference

The normalized input is

```text
f(k;lambda,p) = (p-1)/(pi lambda^2) * (1+k^2/lambda^2)^(-p)
```

with `lambda=1.2 GeV` and `p=3`. It has an algebraic tail and an analytic Hankel
transform proportional to `(lambda*b)^(p-1) K_(p-1)(lambda*b)`. Its analytic
derivative with respect to `lambda` supplies a separate gradient benchmark.

The continuum reference retains the constant spectral endpoint as an exact
identity/contact action and integrates only the compact nonconstant remainder by
adaptive Bessel quadrature. The native path uses a cell-centered radial
finite-volume Laplacian and a Chebyshev-Lobatto polynomial of the resolvent
`R=tau*(tau+T)^(-1)`. Clenshaw recurrence applies the polynomial using sparse
elliptic solves. Both spectral endpoints are interpolation nodes and are retained
to floating-point accuracy.

## Accepted refinement

The accepted final calculation uses degree 1024, radial spacing 0.0125 GeV, and
radius 1280 GeV. Relative to the peak reference bin:

| Quantity | Observed value |
|---|---:|
| Maximum bin discrepancy | 4.51e-5 |
| Maximum input-gradient discrepancy | 1.12e-4 |
| Degree 512 to 1024 change | 2.90e-5 |
| Spacing 0.025 to 0.0125 GeV change | 1.32e-4 |
| Radius 640 to 1280 GeV change | 2.17e-5 |
| Conservative additive bin-error budget | 2.28e-4 |

The bin and conservative-budget acceptance threshold is 4e-4; the gradient
threshold is 7e-4. The adaptive reference reports quadrature error below 3e-13.

The error budget is an observed bin-level refinement budget. It is not a rigorous
global operator-norm certificate. The older shifted-pole least-squares prototype
stalled near a 1.3% sampled weighted symbol error for this uncapped target; it is
not used for the accepted benchmark. The failed development path is recorded in
the method-paper run rather than hidden by the accepted result.

Run:

```bash
python scripts/benchmark_canonical_continuum.py \
  --out outputs/canonical_continuum
```
