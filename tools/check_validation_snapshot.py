#!/usr/bin/env python3
"""Check accepted public evidence, numerical identity and declared accuracy scope."""
from __future__ import annotations
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
manifest = json.loads((ROOT/'validation/SNAPSHOT_MANIFEST.json').read_text())
for name, entry in manifest['files'].items():
    path = (ROOT/name).resolve()
    if not path.is_relative_to(ROOT) or not path.is_file():
        raise SystemExit('Missing or unsafe evidence path: '+name)
    if hashlib.sha256(path.read_bytes()).hexdigest() != entry['sha256']:
        raise SystemExit('Published evidence hash mismatch: '+name)
for name, expected in manifest['numerical_source_sha256'].items():
    if hashlib.sha256((ROOT/name).read_bytes()).hexdigest() != expected:
        raise SystemExit('Numerical source differs from published audit; refresh evidence: '+name)
base = ROOT/'validation/reference/2026-10-01'
load = lambda name: json.loads((base/name).read_text())
suite = load('scientific_suites.json')
assert suite['full_required_run'] and suite['all_suites_passed']
assert suite['suites_passed'] == suite['suites_total'] == 14
completed = load('canonical_physical_completed.json')
assert completed['passed'] and len(completed['cases']) == 3
for case in completed['cases']:
    assert all(case['tests'].values())
    limits = {'bin':5e-4, 'component':5e-4, 'gradient':1e-3, 'contact':5e-4}
    assert all(case['class_wide_error_bounds_over_peak'][key] < limit for key, limit in limits.items())
    assert case['observed_budget'] < 1e-3
    assert len(case['UV_completion_models']) == 2
    assert all(m['passed'] and all(m['tests'].values()) for m in case['UV_completion_models'])
for name in ['canonical_continuum.json', 'support_audit.json', 'primary_cusp_audit.json']:
    assert load(name)['passed'], name
kernel_name='validation/reference/2026-10-01/demonstrators/native_dnn/refined_kernel_validation.json'
kernel_origin=manifest['files'][kernel_name]['origin_sha256']
assert load('demonstrators/native_dnn/refined_kernel_validation.json')['passed']
for name in ['closure.json','data.json']:
    fit=load('demonstrators/native_dnn/'+name)
    assert fit['passed'] and fit['kernel_validation_sha256']==kernel_origin, name
print('PASS: published snapshot hashes, numerical source identity and accepted scoped results')
print('Historical software tests: 121; suites: 14. Fresh execution is a separate check.')
