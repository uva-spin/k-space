# Frozen numerical evidence

This directory publishes the computational evidence for the canonical CSS2
methods paper. The algebra/software evidence comes from audited source commit
`57870ad031ea413576a80bdd5b9ef51f5244c67f`. The full physical calculation was
rerun for publication with the declared 2e-6 conditional UV criterion; exact
executable and report hashes distinguish the refresh. Fresh runs go to
the ignored `outputs/` directory; these records are a versioned reference.

Start with [the results and interpretations](../docs/VALIDATION_RESULTS.md) and
[the reproduction recipe](../docs/REPRODUCE_METHOD_PAPER.md).

| Record | Contents |
|---|---|
| [Scientific suites](reference/2026-10-01/scientific_suites.json) | All 14 required suites, categories, source hashes and report links |
| [Per-suite records](reference/2026-10-01/suites) | Individual algebra, RG, contact, order-counting and numerical results |
| [Ordinary-input continuum](reference/2026-10-01/canonical_continuum.json) | Algebraic-tail input, bins, derivatives and refinement budget |
| [Finite-table physical validation](reference/2026-10-01/canonical_physical_finite_table.json) | Three scales, six components, contacts, gradients and passing 2e-6 conditional UV criterion |
| [Completed physical inputs](reference/2026-10-01/canonical_physical_completed.json) | Two explicit UV inputs, class-wide bounds and original unchanged primary limits |
| [Support audit](reference/2026-10-01/support_audit.json) | Independent RG tolerance, omitted-term and normalization negative controls |
| [Coefficient source comparison](reference/2026-10-01/primary_cusp_audit.json) | Independent four-loop cusp evaluation for nf=3,4,5,6 |
| [Demonstrators](reference/2026-10-01/demonstrators) | Conditional E866 forward/fit and restricted neural closure records |
| [Tables](reference/2026-10-01/tables) | Numerical values underlying the paper figures and summary |

`SNAPSHOT_MANIFEST.json` records published and original report hashes, numerical
source identities, and the publication transformation. Private host paths were
replaced by symbolic locations; computed values and all test outcomes are
unchanged. Every published physical stage passes. The manifest lists the exact
transformation and full original report hash. Development history,
runtime binaries and execution logs are not part of this snapshot. It does not
contain author correspondence, agent instructions or operational handoffs.

The neural demonstration includes both the historical coarse kernel report and
`refined_kernel_validation.json`. The closure/data records reference the latter
by its original hash. Adding this missing refined record repairs the earlier
bundle's provenance gap; no fit prediction or accepted numerical result changed.

```sh
python tools/check_validation_snapshot.py
```

This checks integrity, source identity and the intended scope of the stored
results. It does not rerun the calculations. The completed-input report passes
after resolving specified UV inputs and checking all original primary criteria.
The broader conditional UV bound is 1.6683e-6 of the peak in the most demanding
case, passes the 2e-6 criterion, and stays in the budget.
There is no invented aggregate count of independent physics results and no
assertion of independent human review.
