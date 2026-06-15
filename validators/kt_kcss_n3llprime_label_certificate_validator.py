#!/usr/bin/env python3
"""KCSS N3LL' label certificate validator.

This script distinguishes a formal prescription certificate from an implementation-level
physics claim. The formal prescription can be checked from the note/manifest. The
implementation-level N3LL' label is allowed only when source-stamped physical coefficient
files are present and mock payloads are not used.
"""
from __future__ import annotations

# Repository-local defaults. Override generated-output location with KCSS_OUTPUT_DIR.
import os as _os_for_repo
from pathlib import Path as _PathForRepo
_REPO_ROOT = _PathForRepo(__file__).resolve().parents[1]
_OUTPUT_DIR = _PathForRepo(_os_for_repo.environ.get('KCSS_OUTPUT_DIR', str(_REPO_ROOT / 'outputs'))).resolve()
_OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

import csv, hashlib, json, os
from pathlib import Path
from datetime import datetime, timezone

ROOT = _REPO_ROOT
NOTE = ROOT / 'docs' / 'formalism_note_v1p0.tex'
SOURCE_DIR = ROOT / 'coefficient_sources'

EXPECTED = [
    # arxiv, version, relative_dir, filename, mandatory_for_dy_pdf
    ('2006.05329', 'v2', '2006.05329v2', 'TMDPDF.m', True),
    ('2006.05329', 'v2', '2006.05329v2', 'PTBF.m', False),
    ('2006.05329', 'v2', '2006.05329v2', 'PTBF_ZExpansion.m', False),
    ('2006.05329', 'v2', '2006.05329v2', 'PTBF_ZbExpansion.m', False),
    ('2006.05329', 'v2', '2006.05329v2', 'PT_SoftFunction.m', False),
    ('2006.05329', 'v2', '2006.05329v2', 'PT_SoftFunction_Renormalized.m', False),
    ('2006.05329', 'v2', '2006.05329v2', 'info.txt', False),
    ('2012.03256', 'v2', '2012.03256v2', 'BeamFunction.m', True),
    ('2012.03256', 'v2', '2012.03256v2', 'BeamfunctionN.m', False),
    ('2012.03256', 'v2', '2012.03256v2', 'TMDPDF.m', True),
    ('2012.03256', 'v2', '2012.03256v2', 'TMDPDFN.m', False),
    ('2012.03256', 'v2', '2012.03256v2', 'TMDFF.m', False),
    ('2012.03256', 'v2', '2012.03256v2', 'TMDFFN.m', False),
    ('2012.03256', 'v2', '2012.03256v2', 'FFMatchingKernels.m', False),
    ('2012.03256', 'v2', '2012.03256v2', 'FFMatchingKernelsN.m', False),
    ('2012.03256', 'v2', '2012.03256v2', 'ancillary_readme.pdf', False),
    ('2012.03256', 'v2', '2012.03256v2', 'softfunction.m', False),
]

FORMAL_TOKENS = [
    "\\section{N3LL-prime label certificate}",
    "\\mathrm{ClaimAllowed}_{N^3\\mathrm{LL}'}",
    "\\mathrm{ParentScheme}_{\\KCSS}",
    "\\mathrm{ExpandedTo}\\{a_s^3\\}",
    "N$^3$LL$'_{\\KCSS}",
]


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open('rb') as f:
        for chunk in iter(lambda: f.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest()


def main() -> None:
    text = NOTE.read_text(errors='replace') if NOTE.exists() else ''
    formal_checks = []
    for tok in FORMAL_TOKENS:
        formal_checks.append({'check': f'note_contains:{tok[:40]}', 'passed': tok in text})
    formal_ready = all(c['passed'] for c in formal_checks)

    rows = []
    present = 0
    mandatory_present = 0
    mandatory_total = 0
    hashes = []
    for arxiv, version, d, fname, mandatory in EXPECTED:
        p = SOURCE_DIR / d / fname
        exists = p.exists() and p.stat().st_size > 0
        if exists:
            present += 1
            digest = sha256(p)
            hashes.append((str(p), digest))
        else:
            digest = ''
        if mandatory:
            mandatory_total += 1
            if exists:
                mandatory_present += 1
        rows.append({
            'arxiv': arxiv,
            'version': version,
            'relative_dir': d,
            'filename': fname,
            'mandatory_for_dy_tmdpdf': mandatory,
            'present': exists,
            'sha256': digest,
        })

    physical_payload_ready = mandatory_present == mandatory_total and mandatory_total > 0
    mock_payload_used = False
    expansion_gate_a3_passed = False  # cannot be true until physical parser/expansion is run
    no_numerical_b_transform = True
    np_separated = 'NPSeparated' in text or 'NPSeparated' in text.replace(' ', '')
    implementation_claim_allowed = all([
        formal_ready,
        physical_payload_ready,
        not mock_payload_used,
        expansion_gate_a3_passed,
        no_numerical_b_transform,
        np_separated,
    ])

    tests = list(formal_checks) + [
        {'check': 'formal_prescription_ready', 'passed': formal_ready},
        {'check': 'physical_mandatory_payload_present', 'passed': physical_payload_ready,
         'details': f'{mandatory_present}/{mandatory_total} mandatory files present'},
        {'check': 'mock_payload_not_used_for_physics', 'passed': not mock_payload_used},
        {'check': 'a_s3_expansion_gate_passed', 'passed': expansion_gate_a3_passed,
         'details': 'requires parsed physical coefficients and reference expansion comparison'},
        {'check': 'no_numerical_b_transform_required', 'passed': no_numerical_b_transform},
        {'check': 'np_module_separated_from_accuracy_label', 'passed': np_separated},
        {'check': 'implementation_claim_allowed', 'passed': implementation_claim_allowed},
    ]

    out = {
        'validator': 'kt_kcss_n3llprime_label_certificate_validator',
        'timestamp_utc': datetime.now(timezone.utc).isoformat(),
        'note': str(NOTE),
        'source_dir': str(SOURCE_DIR),
        'formal_prescription_ready': formal_ready,
        'physical_payload_ready': physical_payload_ready,
        'a_s3_expansion_gate_passed': expansion_gate_a3_passed,
        'implementation_claim_allowed': implementation_claim_allowed,
        'present_file_count': present,
        'expected_file_count': len(EXPECTED),
        'mandatory_present_count': mandatory_present,
        'mandatory_file_count': mandatory_total,
        'blocked_reason': None if implementation_claim_allowed else 'Implementation-level N3LL-prime claim remains blocked until physical coefficient files are staged and the O(a_s^3) expansion gate passes.',
        'tests': tests,
        'files': rows,
    }

    (_OUTPUT_DIR / 'kt_kcss_n3llprime_label_certificate_output.json').write_text(json.dumps(out, indent=2))
    with (_OUTPUT_DIR / 'kt_kcss_n3llprime_label_certificate_tests.csv').open('w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=['check','passed','details'])
        w.writeheader()
        for t in tests:
            w.writerow({'check': t.get('check'), 'passed': t.get('passed'), 'details': t.get('details','')})
    with (_OUTPUT_DIR / 'kt_kcss_n3llprime_label_certificate_files.csv').open('w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=['arxiv','version','relative_dir','filename','mandatory_for_dy_tmdpdf','present','sha256'])
        w.writeheader(); w.writerows(rows)
    with (_OUTPUT_DIR / 'kt_kcss_n3llprime_label_certificate_report.txt').open('w') as f:
        f.write("KCSS N3LL' label certificate report\n")
        f.write('====================================\n')
        f.write(f"formal_prescription_ready: {formal_ready}\n")
        f.write(f"physical_payload_ready: {physical_payload_ready}\n")
        f.write(f"a_s3_expansion_gate_passed: {expansion_gate_a3_passed}\n")
        f.write(f"implementation_claim_allowed: {implementation_claim_allowed}\n")
        f.write(f"mandatory files present: {mandatory_present}/{mandatory_total}\n")
        f.write(f"all expected files present: {present}/{len(EXPECTED)}\n")
        if out['blocked_reason']:
            f.write(f"blocked_reason: {out['blocked_reason']}\n")
        f.write('\nTest details:\n')
        for t in tests:
            f.write(f"- {t.get('check')}: {t.get('passed')} {t.get('details','')}\n")

if __name__ == '__main__':
    main()
