#!/usr/bin/env python3
"""
Symbolic validator for the v0.14 N^3LO matching-import interface and manifest-driven
N^3LL' K_CSS2 kernel object.

The mock matching coefficients are placeholders.  This script validates the interface,
projection map, manifest completeness gates, and the three-loop Collins-Soper-kernel RG
identity used by the kernel object.
"""
from __future__ import annotations

# Repository-local defaults. Override generated-output location with KCSS_OUTPUT_DIR.
import os as _os_for_repo
from pathlib import Path as _PathForRepo
_REPO_ROOT = _PathForRepo(__file__).resolve().parents[1]
_OUTPUT_DIR = _PathForRepo(_os_for_repo.environ.get('KCSS_OUTPUT_DIR', str(_REPO_ROOT / 'outputs'))).resolve()
_OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


from dataclasses import dataclass
from pathlib import Path
import csv
import json
from typing import Dict, Iterable, List, Tuple

import sympy as sp

OUT_DIR = _OUTPUT_DIR
PARTON_BASIS = ('g', 'u', 'ubar', 'd', 'dbar', 's', 'sbar')
LEG_TYPES = ('PDF', 'FF')


@dataclass(frozen=True)
class AccuracyManifest:
    label: str
    log_order: int
    prime: bool
    scheme: str = 'CSS2'

    @property
    def required_orders(self) -> Dict[str, int]:
        boundary = self.log_order if self.prime else max(self.log_order - 1, 0)
        return {
            'cusp_loops': self.log_order + 1,
            'beta_loops': self.log_order + 1,
            'gammaF_loops': self.log_order,
            'CS_kernel_loops': self.log_order,
            'hard_boundary_loops': boundary,
            'matching_boundary_loops': boundary,
            'DGLAP_kernel_count': boundary,  # P^(0) ... P^(boundary-1)
        }


@dataclass(frozen=True)
class CoefficientRecord:
    scheme: str
    operator_type: str
    loop_order: int
    to_parton: str
    from_parton: str
    Lb_power: int
    ellzeta_power: int
    z_coefficient: str
    normalization: str = 'a_s=alpha_s/(4pi); Lb=ln(mu^2 b_T^2/b0^2)'
    source: str = 'mock-css2-compatible-schema'

    def project_basis(self) -> str:
        if self.Lb_power == 0:
            return 'deltaT'
        return f'frakL_{self.Lb_power - 1}'


@dataclass
class KernelObject:
    manifest: AccuracyManifest
    coefficient_records: List[CoefficientRecord]
    hard_boundary_loops_available: int
    cusp_loops_available: int
    beta_loops_available: int
    gammaF_loops_available: int
    CS_kernel_loops_available: int
    dglap_kernel_count_available: int
    used_numerical_b_transform: bool = False

    def validate(self) -> Tuple[bool, List[Tuple[str, bool, str]]]:
        req = self.manifest.required_orders
        tests: List[Tuple[str, bool, str]] = []
        max_matching = max((r.loop_order for r in self.coefficient_records), default=-1)
        schemes = {r.scheme for r in self.coefficient_records}
        tests.extend([
            ('scheme_css2', schemes == {'CSS2'}, f'schemes={sorted(schemes)}'),
            ('no_numerical_b_transform', not self.used_numerical_b_transform,
             f'used_numerical_b_transform={self.used_numerical_b_transform}'),
            ('cusp_loop_gate', self.cusp_loops_available >= req['cusp_loops'],
             f"have={self.cusp_loops_available}, need={req['cusp_loops']}"),
            ('beta_loop_gate', self.beta_loops_available >= req['beta_loops'],
             f"have={self.beta_loops_available}, need={req['beta_loops']}"),
            ('gammaF_loop_gate', self.gammaF_loops_available >= req['gammaF_loops'],
             f"have={self.gammaF_loops_available}, need={req['gammaF_loops']}"),
            ('CS_kernel_loop_gate', self.CS_kernel_loops_available >= req['CS_kernel_loops'],
             f"have={self.CS_kernel_loops_available}, need={req['CS_kernel_loops']}"),
            ('hard_boundary_gate', self.hard_boundary_loops_available >= req['hard_boundary_loops'],
             f"have={self.hard_boundary_loops_available}, need={req['hard_boundary_loops']}"),
            ('matching_boundary_gate', max_matching >= req['matching_boundary_loops'],
             f"have={max_matching}, need={req['matching_boundary_loops']}"),
            ('dglap_gate', self.dglap_kernel_count_available >= req['DGLAP_kernel_count'],
             f"have={self.dglap_kernel_count_available}, need={req['DGLAP_kernel_count']}"),
        ])
        return all(ok for _, ok, _ in tests), tests


def make_mock_package(max_loop: int = 3) -> List[CoefficientRecord]:
    records: List[CoefficientRecord] = []
    for operator_type in LEG_TYPES:
        for n in range(max_loop + 1):
            for to_parton in PARTON_BASIS:
                for from_parton in PARTON_BASIS:
                    for m in range(2 * n + 1):
                        records.append(CoefficientRecord(
                            scheme='CSS2',
                            operator_type=operator_type,
                            loop_order=n,
                            to_parton=to_parton,
                            from_parton=from_parton,
                            Lb_power=m,
                            ellzeta_power=0,
                            z_coefficient=f'c_{operator_type}_{n}_{to_parton}_from_{from_parton}_Lb{m}(z)',
                        ))
    return records


def check_power_bounds(records: Iterable[CoefficientRecord]) -> Tuple[bool, str]:
    bad = [r for r in records if not (0 <= r.Lb_power <= 2 * r.loop_order)]
    return len(bad) == 0, f'bad_records={len(bad)}'


def check_coverage(records: List[CoefficientRecord], max_loop: int) -> Tuple[bool, str]:
    have = {(r.operator_type, r.loop_order, r.to_parton, r.from_parton, r.Lb_power) for r in records}
    missing = []
    for operator_type in LEG_TYPES:
        for n in range(max_loop + 1):
            for to_parton in PARTON_BASIS:
                for from_parton in PARTON_BASIS:
                    for m in range(2 * n + 1):
                        key = (operator_type, n, to_parton, from_parton, m)
                        if key not in have:
                            missing.append(key)
    return len(missing) == 0, f'missing={len(missing)}'


def check_projection(records: Iterable[CoefficientRecord]) -> Tuple[bool, str]:
    bad = []
    for r in records:
        expected = 'deltaT' if r.Lb_power == 0 else f'frakL_{r.Lb_power - 1}'
        if r.project_basis() != expected:
            bad.append((r, expected))
    return len(bad) == 0, f'bad_projection={len(bad)}'


def check_flavor_matrix(records: Iterable[CoefficientRecord]) -> Tuple[bool, str]:
    channels = {(r.to_parton, r.from_parton) for r in records if r.operator_type == 'PDF' and r.loop_order == 3}
    needed = {(i, j) for i in PARTON_BASIS for j in PARTON_BASIS}
    missing = needed - channels
    return len(missing) == 0, f'missing_channels={len(missing)}'


def max_direct_frak_index(records: Iterable[CoefficientRecord], loop_order: int) -> int:
    values = [r.Lb_power - 1 for r in records if r.loop_order == loop_order and r.Lb_power > 0]
    return max(values) if values else -1


def check_three_loop_cs_kernel_rg() -> Tuple[bool, str]:
    a, L = sp.symbols('a L')
    G0, G1, G2 = sp.symbols('Gamma0 Gamma1 Gamma2')
    b0, b1 = sp.symbols('beta0 beta1')
    d2, d3 = sp.symbols('d2 d3')
    Ddelta = a**2 * d2 + a**3 * d3
    D0 = a * G0 + a**2 * G1 + a**3 * (G2 + 2 * b0 * d2)
    D1 = a**2 * b0 * G0 / 2 + a**3 * (b0 * G1 + b1 * G0 / 2)
    D2 = a**3 * b0**2 * G0 / 3
    D_b = Ddelta + D0 * L + D1 * L**2 + D2 * L**3
    beta = -2 * b0 * a**2 - 2 * b1 * a**3
    derivative = sp.diff(D_b, L) * 2 + sp.diff(D_b, a) * beta
    target = 2 * (a * G0 + a**2 * G1 + a**3 * G2)
    residual = sp.expand(derivative - target)
    # Truncate at O(a^4) and higher.
    residual_truncated = sp.Poly(residual, a).terms()
    low_terms = []
    for (power_tuple, coeff) in residual_truncated:
        power = power_tuple[0]
        if power <= 3 and sp.simplify(coeff) != 0:
            low_terms.append((power, sp.simplify(coeff)))
    return len(low_terms) == 0, f'low_order_residual_terms={low_terms}'


def write_manifest_csv(path: Path, manifest: AccuracyManifest) -> None:
    with path.open('w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=['label', 'ingredient', 'required'])
        writer.writeheader()
        for k, v in manifest.required_orders.items():
            writer.writerow({'label': manifest.label, 'ingredient': k, 'required': v})
        writer.writerow({'label': manifest.label, 'ingredient': 'DGLAP_indices', 'required': 'P^0,P^1,P^2'})


def write_projection_csv(path: Path, records: List[CoefficientRecord]) -> None:
    chosen = [r for r in records if r.operator_type == 'PDF' and r.to_parton == 'u' and r.from_parton in ('u', 'g')]
    with path.open('w', newline='') as f:
        fields = ['scheme', 'operator_type', 'loop_order', 'to_parton', 'from_parton',
                  'Lb_power', 'ellzeta_power', 'projected_basis', 'z_coefficient']
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        for r in chosen:
            writer.writerow({
                'scheme': r.scheme,
                'operator_type': r.operator_type,
                'loop_order': r.loop_order,
                'to_parton': r.to_parton,
                'from_parton': r.from_parton,
                'Lb_power': r.Lb_power,
                'ellzeta_power': r.ellzeta_power,
                'projected_basis': r.project_basis(),
                'z_coefficient': r.z_coefficient,
            })


def main() -> None:
    manifest_prime = AccuracyManifest(label='N3LLprime_KCSS', log_order=3, prime=True)
    manifest_unprimed = AccuracyManifest(label='N3LL_KCSS', log_order=3, prime=False)
    records_full = make_mock_package(max_loop=3)
    records_nnlo = make_mock_package(max_loop=2)

    tests: List[Tuple[str, bool, str]] = []
    expected_req = {
        'cusp_loops': 4,
        'beta_loops': 4,
        'gammaF_loops': 3,
        'CS_kernel_loops': 3,
        'hard_boundary_loops': 3,
        'matching_boundary_loops': 3,
        'DGLAP_kernel_count': 3,
    }
    tests.append(('manifest_N3LLprime_orders', manifest_prime.required_orders == expected_req,
                  f'req={manifest_prime.required_orders}'))
    tests.append(('manifest_unprimed_boundary_order',
                  manifest_unprimed.required_orders['matching_boundary_loops'] == 2,
                  f"unprimed_matching={manifest_unprimed.required_orders['matching_boundary_loops']}"))
    tests.append(('power_bounds_full_package', *check_power_bounds(records_full)))
    tests.append(('coverage_full_package', *check_coverage(records_full, max_loop=3)))
    tests.append(('projection_rule', *check_projection(records_full)))
    max_frak3 = max_direct_frak_index(records_full, loop_order=3)
    tests.append(('N3LO_Lb6_maps_to_frakL5', max_frak3 == 5, f'max_frak_index={max_frak3}'))
    tests.append(('flavor_matrix_coverage', *check_flavor_matrix(records_full)))
    tests.append(('three_loop_CS_kernel_RG_identity', *check_three_loop_cs_kernel_rg()))

    kernel_prime = KernelObject(
        manifest=manifest_prime,
        coefficient_records=records_full,
        hard_boundary_loops_available=3,
        cusp_loops_available=4,
        beta_loops_available=4,
        gammaF_loops_available=3,
        CS_kernel_loops_available=3,
        dglap_kernel_count_available=3,
    )
    ok_prime, prime_tests = kernel_prime.validate()
    tests.append(('kernel_object_N3LLprime_accepts_complete_package', ok_prime,
                  'all ingredient gates complete'))
    tests.extend((f'kernel_{name}', ok, msg) for name, ok, msg in prime_tests)

    kernel_missing_prime = KernelObject(
        manifest=manifest_prime,
        coefficient_records=records_nnlo,
        hard_boundary_loops_available=2,
        cusp_loops_available=4,
        beta_loops_available=4,
        gammaF_loops_available=3,
        CS_kernel_loops_available=3,
        dglap_kernel_count_available=3,
    )
    ok_missing_prime, _ = kernel_missing_prime.validate()
    tests.append(('missing_three_loop_boundary_rejected_for_prime', not ok_missing_prime,
                  'N3LLprime must reject NNLO boundary package'))

    kernel_unprimed = KernelObject(
        manifest=manifest_unprimed,
        coefficient_records=records_nnlo,
        hard_boundary_loops_available=2,
        cusp_loops_available=4,
        beta_loops_available=4,
        gammaF_loops_available=3,
        CS_kernel_loops_available=3,
        dglap_kernel_count_available=2,
    )
    ok_unprimed, _ = kernel_unprimed.validate()
    tests.append(('NNLO_boundary_accepted_for_unprimed_N3LL', ok_unprimed,
                  'unprimed N3LL boundary gate is loop 2'))

    kernel_bad_transform = KernelObject(
        manifest=manifest_prime,
        coefficient_records=records_full,
        hard_boundary_loops_available=3,
        cusp_loops_available=4,
        beta_loops_available=4,
        gammaF_loops_available=3,
        CS_kernel_loops_available=3,
        dglap_kernel_count_available=3,
        used_numerical_b_transform=True,
    )
    ok_bad_transform, _ = kernel_bad_transform.validate()
    tests.append(('numerical_b_transform_rejected', not ok_bad_transform,
                  'production KCSS object must be coefficient-level'))

    passed = all(ok for _, ok, _ in tests)
    paths = {
        'manifest_csv': OUT_DIR / 'kt_kcss_n3lo_matching_import_manifest.csv',
        'projection_csv': OUT_DIR / 'kt_kcss_n3lo_matching_import_projection.csv',
        'json': OUT_DIR / 'kt_kcss_n3lo_matching_import_output.json',
        'report': OUT_DIR / 'kt_kcss_n3lo_matching_import_validator_report.txt',
        'tests_csv': OUT_DIR / 'kt_kcss_n3lo_matching_import_tests.csv',
    }
    write_manifest_csv(paths['manifest_csv'], manifest_prime)
    write_projection_csv(paths['projection_csv'], records_full)

    with paths['tests_csv'].open('w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=['name', 'passed', 'message'])
        writer.writeheader()
        for name, ok, msg in tests:
            writer.writerow({'name': name, 'passed': ok, 'message': msg})

    result = {
        'passed': passed,
        'num_tests': len(tests),
        'full_mock_records': len(records_full),
        'nnlo_mock_records': len(records_nnlo),
        'parton_basis': PARTON_BASIS,
        'leg_types': LEG_TYPES,
        'N3LLprime_manifest': manifest_prime.required_orders,
        'N3LL_manifest_unprimed': manifest_unprimed.required_orders,
        'max_direct_frak_index_loop3': max_frak3,
        'tests': [{'name': name, 'passed': ok, 'message': msg} for name, ok, msg in tests],
        'output_files': {k: str(v) for k, v in paths.items()},
    }
    paths['json'].write_text(json.dumps(result, indent=2))

    lines = [
        'KCSS N3LO matching-import validator',
        '=' * 54,
        f'passed: {passed}',
        f'num_tests: {len(tests)}',
        f'full_mock_records: {len(records_full)}',
        f'nnlo_mock_records: {len(records_nnlo)}',
        f'parton_basis: {", ".join(PARTON_BASIS)}',
        f'N3LLprime required orders: {manifest_prime.required_orders}',
        f'max direct frakL index at loop 3: {max_frak3}',
        '',
        'Tests:',
    ]
    for name, ok, msg in tests:
        lines.append(f'  [{"PASS" if ok else "FAIL"}] {name}: {msg}')
    paths['report'].write_text('\n'.join(lines) + '\n')
    print('\n'.join(lines))


if __name__ == '__main__':
    main()
