# Reproduce the methods-paper checks

Use Python 3.13 and the pinned dependency files. No author workstation files are
needed for the source tree or offline checks.

```sh
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements-dev.txt
python scripts/run_method_validation.py
```

This runs repository/document checks, snapshot integrity, software tests, all 14
registered suites and the non-Gaussian ordinary-input continuum benchmark.
It explicitly reports that the actual-PDF audit was not requested. Reports and
stage logs go to `outputs/method_validation/` and `outputs/method_suites/`.

## Full canonical physical benchmark

The external requirements are a C++ compiler, make, curl and official LHAPDF
6.5.6 with the pinned NNPDF4.0 member. The explicit setup command downloads and
hash-checks the sources and data, builds outside the repository, and writes
`resources.json`. It does not change a system installation.

```sh
python scripts/setup_lhapdf_local.py --work ../kt-physical-runtime
```

Use the resulting paths; macOS uses `.dylib`, Linux `.so`:

```sh
python scripts/run_method_validation.py --physical \
  --library ../kt-physical-runtime/runtime/libkt_lhapdf.dylib \
  --pdf-data ../kt-physical-runtime/runtime/lhapdf/share/LHAPDF
```

The published reference ran on macOS/Apple silicon. The build recipe supports
Linux, but the archived numerical results are not a claim of a completed Linux
acceptance run. For an existing installation, compile `src/lhapdf_bridge.cc`
against that installation and supply the library/data paths explicitly. The
physical runner checks the PDF bytes against `manifests/css2_pdf.json`.

The full command performs all offline steps plus the degree/spacing/radius study
at three scales, explicit UV input completions and support checks. Expect several
minutes and substantial memory for the 409600-cell high-scale grid. The workflow
is a validation calculation, not an optimized fitting configuration.

The runner first constructs a finite-table intermediate, then integrates the
explicit UV continuations and checks the completed results. The intermediate
retains a historical auxiliary 1e-6 all-class diagnostic, which yields exit 1
at the high scale because its achieved bound is 1.6683e-6. This is separate from
the final acceptance criteria. The next auditor rejects every other failure,
checks source hashes, and requires the unchanged primary criteria and both
completed inputs to pass. Only then can `reproduction.json` report full success.
Intermediate records stay in ignored `outputs/`; accepted results are in
`completed/`. The published snapshot presents the accepted result and achieved
bounds, with its projection recorded in the snapshot manifest.

## Primary-source cusp cross-check

```sh
python scripts/fetch_css2_sources.py --out ../kt-coefficient-sources
python scripts/check_css2_sources.py --sources ../kt-coefficient-sources \
  --out outputs/primary_cusp_audit.json
```

Source archives are identified by versioned URLs and hashes. A changed archive
is rejected, not silently substituted. No downloaded TeX or executable code is
run; the independent ancillary expression is parsed with a restricted grammar.

## Conditional E866 and neural demonstrations

These are separate from the canonical theorem. They require the same explicit
PDF paths and retain the measurement, profile and statistical qualifications in
the [E866](E866_PP_BASELINE.md) and [neural](NATIVE_DNN_PILOT.md) documents.
Use a dedicated output directory and substitute your runtime library path:

```sh
python scripts/benchmark_e866.py --out outputs/e866/forward \
  --library ../kt-physical-runtime/runtime/libkt_lhapdf.dylib \
  --pdf-data ../kt-physical-runtime/runtime/lhapdf/share/LHAPDF
python scripts/e866_data_diagnostic.py --out outputs/e866/fit \
  --validation outputs/e866/forward/e866_forward.json \
  --library ../kt-physical-runtime/runtime/libkt_lhapdf.dylib \
  --pdf-data ../kt-physical-runtime/runtime/lhapdf/share/LHAPDF
python scripts/build_native_dnn_kernel.py --out outputs/pilot/kernel_dense \
  --nQ 8 --nx 16 \
  --library ../kt-physical-runtime/runtime/libkt_lhapdf.dylib \
  --pdf-data ../kt-physical-runtime/runtime/lhapdf/share/LHAPDF
python scripts/run_native_dnn.py --stage closure --seeds 11 --maxiter 60000 \
  --evidence-root outputs/pilot --out outputs/pilot/closure
python scripts/run_native_dnn.py --stage data --seeds 7,19,31,43,59 --maxiter 20000 \
  --evidence-root outputs/pilot --out outputs/pilot/data_initial
python scripts/polish_dnn.py --input outputs/pilot/data_initial \
  --evidence-root outputs/pilot --out outputs/pilot/data
python scripts/audit_dnn_pilot.py --evidence-root outputs/pilot
```

The data stage requires closure to pass first. The continuation step preserves
the objective and thresholds when the first fit is iteration-limited. Seeds,
optimizer status, source hashes and numerical gradient comparisons are recorded.
If continuation is still iteration-limited, repeat `polish_dnn.py` from its last
saved output into a new directory. Accept it as `outputs/pilot/data` only when
every optimizer and numerical gate passes; do not overwrite the failed stages.
These stages can be much more expensive than the canonical benchmark. Their
frozen reference outputs are already included in `validation/`.

## Figures, records and CI

Scientific plot PDFs/PNGs are under `figures/`; CSV values and all signed
components are under `validation/reference/2026-10-01/`. The plotting command is:

```sh
python -m pip install matplotlib
python scripts/plot_method_validation.py
```

Default CI runs software tests, the offline suites and snapshot checks on each
push/PR. The manual **full canonical physical validation** workflow builds the
external PDF runtime and runs the full command. Its availability does not imply
that it has run or passed; inspect the specific Actions result. Both workflows
retain generated results as downloadable artifacts.

Numerical tolerances, not byte-identical floating-point results, define acceptance
on another platform. The snapshot manifest preserves exact historical identities.
