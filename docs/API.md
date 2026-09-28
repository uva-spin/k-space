# Operator API

The public computational module is [`src/resummed_operator.py`](../src/resummed_operator.py).
It is usable directly from a checkout or through an editable install:

```bash
python -m pip install --no-deps -e .
```

The installation exposes the module `resummed_operator`. The examples modify
their local import path explicitly, so an editable install is optional for
`python examples/...` commands.

## Coefficients

```python
from resummed_operator import Coefficients
coefficients = Coefficients.load("data/coefficients_demo.json")
```

Fields `beta`, `cusp`, `noncusp`, and `d` are tuples of numerical coefficients.
The beta polynomial begins at a²; the other series begin at a. The CS boundary
array stores d₁, d₂, d₃ in order, with d₁=0 required by this implementation.
The supported maximum lengths are 4, 4, 3, and 3, respectively. The loader checks
finite values, positive leading beta/cusp coefficients, and d₁=0.

The demonstration provides two beta and cusp terms, one noncusp term, and d₁=0,d₂.
`physical_n3llprime_payload` is descriptive metadata, not independent evidence of
coefficient provenance or physical accuracy. The reference examples explicitly
stay within the demonstration scope. See [COEFFICIENTS.md](COEFFICIENTS.md).

## CanonicalRG

```python
from resummed_operator import CanonicalRG
rg = CanonicalRG(coefficients, mu0=2.0, a0=0.02)
```

The constructor solves a retained, unexpanded ODE system in
`t=ln(nu/mu0)` using DOP853 and dense output. Defaults are `tmax=1200`,
`rtol=2e-12`, and `atol=2e-14`. The beta polynomial must remain positive on the
selected asymptotically free branch up to a₀. ODE tolerances are solver settings,
not a global physical uncertainty bound. Requests outside the solved logarithmic
interval are rejected.

| Method | Meaning / return value |
|---|---|
| `values(t)` | Array of a(t), integrated cusp, integrated noncusp, and t-weighted integrated cusp |
| `tc(s)` | ln(μ_c(s)/μ₀); s must be positive |
| `D(s, mu)` | Completed CS symbol at positive spectral s |
| `v(s, mu, zeta)` | Complete canonical leg; supports the analytic value v(0)=0 |
| `log_u(s, target, initial)` | Log endpoint evolution; each endpoint is `(mu, zeta)`; use positive s |
| `h_nodes(degree, tau, mu, zeta)` | `(M_infinity, h_nodes)` for the contact-preserving resolvent polynomial |

Use finite real inputs and admitted scale/spectral domains. s and tau have units
GeV⁻², μ has units GeV, and ζ has units GeV². At s=0 the endpoint operator's domain
is subtler than the bounded canonical leg's endpoint: use `v(0)`, not an arbitrary
finite substitute for `D(0)` or `log_u(0)`.

The formulas, infrared completion, and common operator domain are documented in
[FORMALISM.md](FORMALISM.md). Functions returning scalar symbols do not transform
a fitted momentum array.

## Apply the resolvent polynomial

```python
from resummed_operator import apply_resolvent_polynomial
M_infinity, nodes = rg.h_nodes(128, tau, mu, zeta)
evolved = apply_resolvent_polynomial(T, input_array, tau, M_infinity, nodes)
```

`T` is a square sparse discretization of the nonnegative momentum-space operator
−Δ_k. `input_array` is a real vector of matching length; `tau` is positive.
The caller supplies a discretization satisfying the intended self-adjointness and
boundary assumptions. The routine checks basic shape and tau consistency; it does
not prove those operator assumptions from an arbitrary matrix.

The recurrence evaluates `M_infinity*I + R*B_N[h](R)` with
`R=tau*(tau*I+T)^(-1)`. It reuses a sparse factorization and uses no FFT or
eigensolver. It does not silently estimate continuum error from a finite grid.
The grid example shows a valid small Dirichlet setup and an optional independent
same-grid comparison.

`bernstein_scalar(nodes, r)` evaluates the same scalar polynomial by de Casteljau
recursion. It is used in numerical diagnostics; evaluation of the fitted array
remains through sparse momentum-space solves.

## Stable scientific source

The core module, 12 scientific validators, and numerical coefficient fixture in
this package preserve the supplied source bytes. New runnable examples and
software tests call that module. Check the imported hashes with:

```bash
python tools/check_repository_policy.py --check-prepared-hashes
```

A future scientific edit should update tests, the source ledger, and interpretation
of results rather than retaining a byte-identity claim after a formula changes.
