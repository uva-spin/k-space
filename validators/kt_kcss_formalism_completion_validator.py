#!/usr/bin/env python3
"""Completion checks for the paper-ready pure-kT formalism note."""
from __future__ import annotations

# Repository-local defaults. Override generated-output location with KCSS_OUTPUT_DIR.
import os as _os_for_repo
from pathlib import Path as _PathForRepo
_REPO_ROOT = _PathForRepo(__file__).resolve().parents[1]
_OUTPUT_DIR = _PathForRepo(_os_for_repo.environ.get('KCSS_OUTPUT_DIR', str(_REPO_ROOT / 'outputs'))).resolve()
_OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


import csv
import hashlib
import json
from pathlib import Path

ROOT = _REPO_ROOT
TEX = ROOT / 'docs' / 'formalism_note_v1p0.tex'
PDF = ROOT / 'docs' / 'formalism_note_v1p0.pdf'

text = TEX.read_text()
checks = []

def add(name: str, passed: bool, detail: str = '') -> None:
    checks.append({'check': name, 'passed': bool(passed), 'detail': detail})

required_strings = {
    'version_1p0': r'\date{Version 1.0 -- June 14, 2026}',
    'generic_process_section': r'\section{Generic process-level $W$ operator}',
    'model_agnostic_np_section': r'\section{Model-agnostic nonperturbative module}',
    'regulated_norm_section': r'\section{Regulated normalization and transverse moments}',
    'full_spectrum_matching_section': r'\section{Full-spectrum matching and perturbative uncertainties}',
    'claim_table_section': r'\section{Paper-ready claim table and extraction workflow}',
    'validation_section': r'\section{Validation sequence}',
    'closure_section': r'\section{Closure of the formalism note}',
    'n3ll_certificate': r"\mathrm{ClaimAllowed}_{N^3\mathrm{LL}'}",
    'formal_vs_implementation': 'Formal prescription label',
    'np_not_accuracy_source': r'DNNs are a flexible implementation of $\mathcal M_{\NP}$',
    'tmd_tail_vs_y_tail': r'\boxed{\text{TMD large-}k_T\text{ tail}}',
    'wy_formula': r'\dd\sigma_Y',
    'regulated_cumulants': r'\mathcal I_w[F_i]',
    'hybrid_numerical_decision': 'The paper-level recommendation is the hybrid approach',
}
for name, needle in required_strings.items():
    add(name, needle in text, f'found={needle in text}')

add('no_active_todo_calls', r'\todo{' not in text, 'no active \\todo{...} boxes')
add('no_immediate_next_tasks_section', r'\section{Immediate next tasks}' not in text, 'formal note closed')
add('physical_payload_not_overclaimed', r"Implementation-level N$^3$LL$'$ claim: \textbf{blocked" in text, 'implementation gate remains blocked')
add('pdf_exists', PDF.exists(), str(PDF))
add('tex_exists', TEX.exists(), str(TEX))

sha_tex = hashlib.sha256(TEX.read_bytes()).hexdigest() if TEX.exists() else None
sha_pdf = hashlib.sha256(PDF.read_bytes()).hexdigest() if PDF.exists() else None
passed = all(row['passed'] for row in checks)
summary = {
    'passed': passed,
    'num_checks': len(checks),
    'num_passed': sum(row['passed'] for row in checks),
    'tex_sha256': sha_tex,
    'pdf_sha256': sha_pdf,
    'formalism_note_status': 'paper-ready formal prescription' if passed else 'incomplete',
    'implementation_N3LLprime_claim_allowed': False,
    'implementation_gate_reason': 'physical NNLO/N3LO coefficient payloads and a_s^3 expansion gate are intentionally not completed in this formalism note',
}

(_OUTPUT_DIR / 'kt_kcss_formalism_completion_validator_output.json').write_text(json.dumps({'summary': summary, 'checks': checks}, indent=2))
with (_OUTPUT_DIR / 'kt_kcss_formalism_completion_validator_tests.csv').open('w', newline='') as f:
    writer = csv.DictWriter(f, fieldnames=['check', 'passed', 'detail'])
    writer.writeheader()
    writer.writerows(checks)

lines = [
    'KCSS formalism completion validator',
    '====================================',
    f"passed: {passed}",
    f"checks: {summary['num_passed']}/{summary['num_checks']}",
    f"tex_sha256: {sha_tex}",
    f"pdf_sha256: {sha_pdf}",
    '',
    'Interpretation:',
    '- The formalism note is ready to serve as the basis for a first paper if all checks pass.',
    '- The implementation-level N3LL-prime extraction claim remains blocked until the physical coefficient payloads and a_s^3 expansion gate are completed.',
    '',
    'Failed checks:'
]
failed = [row for row in checks if not row['passed']]
if failed:
    for row in failed:
        lines.append(f"- {row['check']}: {row['detail']}")
else:
    lines.append('- none')
(_OUTPUT_DIR / 'kt_kcss_formalism_completion_validator_report.txt').write_text('\n'.join(lines) + '\n')
print(json.dumps(summary, indent=2))
