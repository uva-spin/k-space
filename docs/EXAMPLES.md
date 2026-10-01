# Runnable examples

The two examples are thin, executable clients of `src/resummed_operator.py`.
Their parameters and output are described here so the README commands can be
used without reading a validator first. Run them from the repository root with
the environment in [REPRODUCIBILITY.md](REPRODUCIBILITY.md) activated.

## Common inputs and units

| Option | Default | Meaning |
|---|---:|---|
| `--coefficients` | `data/coefficients_demo.json` | Explicit lower-order demonstration input |
| `--mu0` | 2 | Infrared completion scale, GeV |
| `--a0` | 0.02 | α_s(μ₀)/(4π), dimensionless |
| `--mu` | 20 | Final virtuality scale, GeV; μ ≥ μ₀ |
| `--zeta` | 400 | Final rapidity scale, GeV² |
| `--out` | Example-specific child of `outputs/` | Output directory; files from that example may be overwritten on a repeat run |

The examples require finite parameters in the implementation's domain and use
its declared fixed-flavor completion. They do not activate a physical N³LL′
claim if a JSON flag is changed. See [COEFFICIENTS.md](COEFFICIENTS.md).

To exercise the sourced physical evolution coefficients, pass
`--coefficients data/coefficients_css2_nf5_n3ll.json`. The generated JSON then
records `evolution_accuracy: N3LL`, `hard_and_ope_order: strict NLO`, the source
manifest hash, and `physical_N3LLprime_prediction_validated: false`.

### Canonical-symbol table

```bash
python examples/canonical_evolution.py --mu 20 --zeta 400 --s 0 0.01 0.1 1
python examples/canonical_evolution.py --coefficients data/coefficients_css2_nf5_n3ll.json
python examples/canonical_evolution.py --help
```

`--s` takes a list of nonnegative spectral resolutions in GeV⁻². It is not a
list of transverse momenta. The example calls `CanonicalRG.v` and, for positive
s, `CanonicalRG.tc` and `CanonicalRG.D`.

It writes `canonical_evolution.csv` and `canonical_evolution.json` below
`outputs/examples/canonical/` by default. CSV columns are:

| Column | Interpretation |
|---|---|
| `s_GeV_minus2` | Spectral resolution |
| `log_mu_c_over_mu0` | Canonical-scale logarithm; absent at s=0 |
| `mu_c_GeV` | Canonical/completed scale; absent at s=0 |
| `D_c` | Resummed CS symbol; absent at s=0 |
| `canonical_leg` | Dimensionless canonical evolution symbol; analytic value 0 at s=0 |

The unassigned fields at s=0 are empty in CSV and null in JSON. They are not
replaced by fictitious finite endpoint values. The JSON records coefficient and
operator hashes, dependency versions, and the input scales.

### Apply evolution to a finite momentum grid

```bash
# Momentum-space recurrence only.
python examples/momentum_grid.py

# Add an independent same-grid validation reference.
python examples/momentum_grid.py --reference

# An alternative numerical setup, still a demonstration.
python examples/momentum_grid.py --grid 13 --halfwidth 8 --degrees 16 64 \
  --reference --out outputs/examples/grid13
```

The default 17×17 grid contains interior nodes of a square with half-width 8 GeV.
The spacing is `2*halfwidth/(grid+1)` and the finite-difference Laplacian has
Dirichlet boundaries. The Gaussian `exp(−k²/2)` is normalized in the weighted
L² norm, `norm(f)*spacing=1`. This is a normalization for the numerical error
comparison, not a TMD-to-collinear matching condition.

For each degree N, `CanonicalRG.h_nodes` constructs the scalar polynomial data
and `apply_resolvent_polynomial` applies it using sparse solves. The default
resolvent parameter is `tau=b0²/mu0²`, with units GeV⁻². Without `--reference`,
there is no diagonalization and no numerical error-bound claim. With that flag,
an independent dense diagonalization of the same finite matrix is used to check
the result and its matrix-polynomial error bound. The dense reference is limited
to 25 nodes per axis in this example to keep it a small diagnostic.

Outputs below `outputs/examples/grid/`:

| File | Contents |
|---|---|
| `degree_convergence.csv` | Degree, input/output L² norms, and optional reference errors/bounds |
| `momentum_grid.csv` | k_x, k_y, original input, evolved array at the largest requested degree |
| `momentum_grid.json` | Setup, hashes, scope, and all degree records |

The error `L2_same_grid_error` compares against the unapproximated function of
the **same** finite Laplacian. `recurrence_oracle_residual` compares the sparse
recurrence against the polynomial evaluated by the same-grid diagonalization.
`finite_matrix_bound` bounds the former polynomial error for this matrix and
normalized input. No quantity in this table measures finite-box or continuum
spacing error. The full registered `unexpanded` suite provides more degree and
endpoint checks than the small example.

## Error handling and reruns

Invalid scales, negative/nonfinite s, unsupported grid/degree sizes, nonfinite
outputs, or a failed reference comparison produce a nonzero exit. Outputs must
be in a child of this checkout's `outputs/` directory and must not use symlinks.
Examples never recursively delete directories. A repeat run replaces their named
CSV/JSON outputs; use a different `--out` to keep a second parameter choice.

This example layer supplies convenient execution and interpretation. It does not
turn the N³LL evolution input into an N³LL′ boundary or process-spectrum claim.
