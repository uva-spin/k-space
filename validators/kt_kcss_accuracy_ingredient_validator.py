#!/usr/bin/env python3
"""
kt_kcss_accuracy_ingredient_validator.py

Manifest validator for the accuracy convention used in
kt_perturbative_formalism_note_v0p13.tex.

It does not validate external perturbative coefficients.  It validates that the internal
bookkeeping convention maps N^kLL and N^kLL' labels to the loop orders stated in the note.
"""
from __future__ import annotations

# Repository-local defaults. Override generated-output location with KCSS_OUTPUT_DIR.
import os as _os_for_repo
from pathlib import Path as _PathForRepo
_REPO_ROOT = _PathForRepo(__file__).resolve().parents[1]
_OUTPUT_DIR = _PathForRepo(_os_for_repo.environ.get('KCSS_OUTPUT_DIR', str(_REPO_ROOT / 'outputs'))).resolve()
_OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


import csv
import json
from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Dict, List


@dataclass(frozen=True)
class AccuracyOrder:
    label: str
    k: int
    primed: bool
    cusp_loops: int
    beta_loops: int
    noncusp_loops: int
    cs_kernel_loops: int
    hard_boundary_loops: int
    tmd_matching_loops: int
    dglap_max_P_index: int | None
    fixed_order_Y_power: int | None


def boundary_order(k: int, primed: bool) -> int:
    """Loop order of hard/TMD boundary functions for the convention in v0.13."""
    if primed:
        return k
    return max(k - 1, 0)


def manifest(label: str, k: int, primed: bool) -> AccuracyOrder:
    bnd = boundary_order(k, primed)
    # To reconstruct an m-loop matching coefficient's mu dependence one needs
    # P^(0) ... P^(m-1).  If m=0, no DGLAP kernel is part of the boundary logs.
    dglap_max = None if bnd == 0 else bnd - 1
    return AccuracyOrder(
        label=label,
        k=k,
        primed=primed,
        cusp_loops=k + 1,
        beta_loops=k + 1,
        noncusp_loops=k,
        cs_kernel_loops=k,
        hard_boundary_loops=bnd,
        tmd_matching_loops=bnd,
        dglap_max_P_index=dglap_max,
        fixed_order_Y_power=bnd if bnd > 0 else None,
    )


def main() -> None:
    outdir = _OUTPUT_DIR
    rows: List[AccuracyOrder] = [
        manifest('LL', 0, False),
        manifest('NLL', 1, False),
        manifest("NLL'", 1, True),
        manifest('NNLL', 2, False),
        manifest("NNLL'", 2, True),
        manifest('N3LL', 3, False),
        manifest("N3LL'", 3, True),
    ]

    target = manifest("N3LL'", 3, True)
    expected: Dict[str, object] = {
        'cusp_loops': 4,
        'beta_loops': 4,
        'noncusp_loops': 3,
        'cs_kernel_loops': 3,
        'hard_boundary_loops': 3,
        'tmd_matching_loops': 3,
        'dglap_max_P_index': 2,
        'fixed_order_Y_power': 3,
    }

    tests = []
    for key, val in expected.items():
        got = getattr(target, key)
        tests.append({'test': f'N3LLprime_{key}', 'expected': val, 'got': got, 'passed': got == val})

    # Generic consistency checks for the full table.
    for row in rows:
        tests.append({
            'test': f'{row.label}_cusp_one_above_noncusp',
            'expected': 1,
            'got': row.cusp_loops - row.noncusp_loops,
            'passed': (row.cusp_loops - row.noncusp_loops) == 1,
        })
        tests.append({
            'test': f'{row.label}_beta_equals_cusp_loop_count',
            'expected': row.cusp_loops,
            'got': row.beta_loops,
            'passed': row.beta_loops == row.cusp_loops,
        })
        tests.append({
            'test': f'{row.label}_hard_equals_matching_boundary',
            'expected': row.hard_boundary_loops,
            'got': row.tmd_matching_loops,
            'passed': row.hard_boundary_loops == row.tmd_matching_loops,
        })
        if row.hard_boundary_loops == 0:
            expected_p = None
        else:
            expected_p = row.hard_boundary_loops - 1
        tests.append({
            'test': f'{row.label}_dglap_index_from_boundary',
            'expected': expected_p,
            'got': row.dglap_max_P_index,
            'passed': row.dglap_max_P_index == expected_p,
        })

    passed = all(t['passed'] for t in tests)

    # CSV manifest.
    manifest_path = outdir / 'kt_kcss_accuracy_ingredient_manifest.csv'
    with manifest_path.open('w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=list(asdict(rows[0]).keys()))
        writer.writeheader()
        for row in rows:
            writer.writerow(asdict(row))

    tests_path = outdir / 'kt_kcss_accuracy_ingredient_validator_tests.csv'
    with tests_path.open('w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=['test', 'expected', 'got', 'passed'])
        writer.writeheader()
        writer.writerows(tests)

    output = {
        'passed': passed,
        'target': asdict(target),
        'expected_N3LLprime': expected,
        'num_tests': len(tests),
        'num_failed': sum(1 for t in tests if not t['passed']),
        'tests': tests,
    }
    json_path = outdir / 'kt_kcss_accuracy_ingredient_validator_output.json'
    json_path.write_text(json.dumps(output, indent=2))

    report_path = outdir / 'kt_kcss_accuracy_ingredient_validator_report.txt'
    lines = [
        'KCSS accuracy ingredient manifest validator',
        '================================================',
        f'passed: {passed}',
        f'num_tests: {len(tests)}',
        f'num_failed: {output["num_failed"]}',
        '',
        "N3LL' target manifest:",
    ]
    for key, val in asdict(target).items():
        lines.append(f'  {key}: {val}')
    lines.extend(['', 'Failed tests:'])
    failed = [t for t in tests if not t['passed']]
    if failed:
        for t in failed:
            lines.append(f'  {t["test"]}: expected {t["expected"]}, got {t["got"]}')
    else:
        lines.append('  none')
    report_path.write_text('\n'.join(lines) + '\n')

    print(report_path.read_text())


if __name__ == '__main__':
    main()
