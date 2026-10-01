# Reproducibility

## Installation and execution

Use Python 3.13 with NumPy 2.3.5, SciPy 1.17.0, SymPy 1.14.0, and mpmath 1.3.0.
The scientific dependencies are pinned in
`requirements.txt`; software-contract tests additionally use pytest 9.0.2.

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements-dev.txt
python tools/check_repository_policy.py
python tools/check_documentation.py
python -m pytest -q
python run_all_checks.py --clean
```

After dependency installation, the checks run offline. No TeX executable,
manuscript source, external TMD table, API key, or PDF is needed.

## Runnable example clients

See [EXAMPLES.md](EXAMPLES.md) for two executable clients and [API.md](API.md)
for argument domains and units. These clients call the shared operator module.
Their outputs are under `outputs/examples/` and do not change the 13-suite registry.
Software tests of example arguments, output shapes, and source consistency are
reported by pytest, separately from the scientific suites.

## Output ownership

The runner creates a fresh subdirectory for every executed suite, even without
`--clean`. It captures stdout/stderr into `console.log` and requires a newly
written structured report. A nonzero process exit, malformed report, missing
success field, failed individual test, skipped required work, or unsupported
accuracy claim causes a nonzero aggregate exit. A subprocess adapter binds each report digest to its suite
and fresh run UUID; the runner checks that envelope before accepting the report.
The manifest lists required IDs and must agree with discovery of validator scripts.
Selected runs explicitly record unselected suites and `full_required_run: false`;
all registered suites are required in the advertised full run.

```text
outputs/latest/
  run_all_checks_report.txt
  run_all_checks_summary.json
  structural/
    audit_results.json
    audit_tests.csv
    contact_complete_bins.csv
    console.log
  unexpanded/
    unexpanded_validation.json
    unexpanded_tests.csv
    log_tower_tests.csv
    numeric_scope.json
    resolvent_grid_results.csv
    console.log
  physical_canonical/
    physical_canonical_validation.json
    console.log
  ...
```

Use a distinct, repository-local output directory to preserve another run:

```bash
python run_all_checks.py --output-dir outputs/repeat --clean
```

The runner restricts cleanup to a child of `outputs/` and rejects symlinked cleanup
paths. It records the Python executable environment, package versions, thread
settings, execution time, input hashes, and Git commit/dirty state when available.

## Source provenance

`manifests/source_provenance.json` identifies the recovered source archives and
maps each imported scientific file to its original archive member and SHA256.
The earlier component validators remain byte-identical imports. The operator
retains its imported canonical construction and adds accuracy metadata plus the
Chebyshev-Lobatto/Clenshaw resolvent evaluator used by the continuum benchmark.
The structural and unexpanded validator ports replace manuscript hash references
with the source ledger and remove manuscript output; their scientific test bodies
are retained. These adaptations are explicitly recorded.

```bash
python tools/check_repository_policy.py --check-prepared-hashes
```

This last command verifies the recorded imported scientific files. After a deliberate
scientific edit, retain the source lineage and record the new hash/adaptation rather than
claiming byte identity. The ordinary policy/CI check does not require every future
revision to preserve the original imported hashes.

## Generated reports

Run reports are generated locally under `outputs/` and are excluded from Git.
CI uploads reports and logs as workflow artifacts even when a check fails. Read
the structured residuals and test-specific tolerances rather than expecting
bit-for-bit equality of platform-dependent floating-point results.

A changed coefficient value or tolerance requires a scientific justification.
A changed test total should be explained by the added or removed tests.

Repository policy excludes LaTeX/BibTeX sources, draft PDFs, build products, nested
archives, agent instructions, automation state, and generated report files from
the tracked tree. Scientific text data/reports and PDF plots under `figures/` are
allowed. Review plot contents before publication. The required dependency lists
remain `requirements.txt` and `requirements-dev.txt`; explanatory documentation
uses Markdown.

## Same-grid and continuum comparisons

The lower-order unexpanded demo uses a 17×17 interior grid, Dirichlet boundary
conditions, a half-width of 8 GeV, and polynomial degrees 8, 16, 32, 64, 128.
Its independent diagonalization and recurrence operate on the **same matrix**.
A separate program must vary box size and spacing and derive/measure a continuum
error before publishing a continuum numerical result.

## CI and fresh-clone acceptance

CI installs dependencies, checks repository policy, runs software-contract tests,
executes all scientific/component suites, and retains reports as artifacts.
A fresh-clone acceptance run must succeed with no file copied from the writing
workspace. Inspect the tracked inventory and run the repository policy checker
before publishing changes.
