#!/usr/bin/env python3
"""Continuum Hankel benchmark of physical canonical N3LL evolution.

The production path applies the rational resolvent to a radial momentum grid.
The independent reference uses the analytic Hankel transform of a normalized
Student-type radial input and adaptive Bessel quadrature for identical bins.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from math import gamma
from pathlib import Path
import sys

import numpy as np
from scipy.integrate import quad_vec
from scipy.special import j1, kv

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'src'))
from css2_coefficients import qcd_coefficients  # noqa: E402
from rational_radial import RadialGrid  # noqa: E402
from resummed_operator import CanonicalRG, apply_chebyshev_resolvent  # noqa: E402


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def student_input(k, scale, power):
    u = (np.asarray(k)/scale)**2
    return (power-1)/(np.pi*scale**2)*(1+u)**(-power)


def student_input_derivative(k, scale, power):
    k = np.asarray(k)
    u = (k/scale)**2
    return student_input(k, scale, power)*(-2+2*power*u/(1+u))/scale


def hankel_input(b, scale, power):
    b = np.asarray(b, dtype=float)
    z = scale*b
    value = np.ones_like(z)
    positive = z > 0
    coefficient = 2**(2-power)/gamma(power-1)
    value[positive] = coefficient*z[positive]**(power-1)*kv(power-1, z[positive])
    return value


def hankel_input_derivative(b, scale, power):
    b = np.asarray(b, dtype=float)
    z = scale*b
    value = np.zeros_like(z)
    positive = z > 0
    coefficient = 2**(2-power)/gamma(power-1)
    value[positive] = (-coefficient*b[positive]*z[positive]**(power-1)
                       * kv(power-2, z[positive]))
    return value


def input_bins(edges, scale, power):
    a, b = edges[:-1], edges[1:]
    return (1+(a/scale)**2)**(1-power) - (1+(b/scale)**2)**(1-power)


def input_bin_derivative(edges, scale, power):
    def primitive_derivative(radius):
        u = (radius/scale)**2
        return 2*(power-1)*radius**2/scale**3*(1+u)**(-power)
    return primitive_derivative(edges[:-1])-primitive_derivative(edges[1:])


def continuum_bins(symbol, s0, edges, scale, power, derivative=False):
    infinity = float(symbol(np.array([2*s0]))[0])
    contact = input_bin_derivative(edges, scale, power) if derivative else input_bins(edges, scale, power)
    transform = hankel_input_derivative if derivative else hankel_input
    a, b = edges[:-1], edges[1:]
    def integrand(radial_b):
        kernel = b*j1(b*radial_b)-a*j1(a*radial_b)
        return ((float(symbol(np.array([radial_b**2]))[0])-infinity)
                * float(transform(np.array([radial_b]), scale, power)[0])*kernel)
    transition = np.sqrt(s0/2)
    value, error = quad_vec(integrand, 0., np.sqrt(s0), points=[transition],
                            epsabs=2e-12, epsrel=2e-12, limit=300)
    return infinity*contact+value, float(error)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--scale', type=float, default=1.2, help='Student input scale in GeV')
    parser.add_argument('--power', type=int, default=3, help='Student tail power, at least 3')
    args = parser.parse_args(argv)
    if args.scale <= 0 or args.power < 3:
        parser.error('Require scale>0 and power>=3')
    args.out.mkdir(parents=True, exist_ok=True)

    coefficients = qcd_coefficients(5)
    rg = CanonicalRG(coefficients)
    symbol = lambda s: rg.v(s, 20., 400.)
    edges = np.array([0., .5, 1., 2., 3., 4., 6., 8.])
    reference, reference_error = continuum_bins(
        symbol, rg.s0, edges, args.scale, args.power)
    gradient_reference, gradient_reference_error = continuum_bins(
        symbol, rg.s0, edges, args.scale, args.power, derivative=True)
    normalization = max(float(np.max(np.abs(reference))), 1e-14)

    configurations = [
        (256, 1280., .0125), (512, 1280., .0125), (1024, 1280., .0125),
        (1024, 1280., .05), (1024, 1280., .025),
        (1024, 320., .0125), (1024, 640., .0125),
    ]
    polynomials = {}
    rows = []
    for degree, radius, spacing in configurations:
        coefficients = polynomials.setdefault(
            degree, rg.chebyshev_coefficients(degree, rg.s0, 20., 400.))
        grid = RadialGrid(radius, round(radius/spacing))
        values = np.stack([
            student_input(grid.k, args.scale, args.power),
            student_input_derivative(grid.k, args.scale, args.power),
        ], axis=1)
        weighted = np.sqrt(grid.weights)[:, None]*values
        evolved = apply_chebyshev_resolvent(
            grid.matrix(), weighted, rg.s0, coefficients)/np.sqrt(grid.weights)[:, None]
        bins = grid.bins(evolved, edges)
        prediction, gradient = bins[:, 0], bins[:, 1]
        row = {
            'degree': degree, 'radius_GeV': radius, 'spacing_GeV': spacing,
            'bins': prediction.tolist(), 'gradient_bins': gradient.tolist(),
            'max_bin_error_over_peak': float(np.max(np.abs(prediction-reference))/normalization),
            'max_gradient_error_over_peak': float(np.max(np.abs(gradient-gradient_reference))/normalization),
            'polynomial_basis': 'Chebyshev-Lobatto in r=tau/(tau+s)',
        }
        rows.append(row)
        print('REFINE', degree, radius, spacing,
              row['max_bin_error_over_peak'], row['max_gradient_error_over_peak'], flush=True)

    by_key = {(r['degree'], r['radius_GeV'], r['spacing_GeV']): r for r in rows}
    final = by_key[(1024, 1280., .0125)]
    degree_change = float(np.max(np.abs(
        np.asarray(final['bins'])-np.asarray(by_key[(512, 1280., .0125)]['bins'])))/normalization)
    spacing_change = float(np.max(np.abs(
        np.asarray(final['bins'])-np.asarray(by_key[(1024, 1280., .025)]['bins'])))/normalization)
    box_change = float(np.max(np.abs(
        np.asarray(final['bins'])-np.asarray(by_key[(1024, 640., .0125)]['bins'])))/normalization)

    # Verify the analytic input derivative before interpreting evolved gradients.
    step = 2e-5
    grid = RadialGrid(20., 800)
    finite_difference = (student_input(grid.k, args.scale+step, args.power)
                         - student_input(grid.k, args.scale-step, args.power))/(2*step)
    derivative_residual = float(np.max(np.abs(
        finite_difference-student_input_derivative(grid.k, args.scale, args.power))))
    # A conservative observed error budget adds the final reference discrepancy
    # and the last independent degree, spacing, and box changes. This is an
    # empirical bin-level budget, not a rigorous continuum operator norm.
    observed_error_budget = (final['max_bin_error_over_peak'] + degree_change
                             + spacing_change + box_change)

    tests = [
        {'name': 'reference quadrature', 'passed': max(reference_error, gradient_reference_error) < 2e-10,
         'value': max(reference_error, gradient_reference_error), 'threshold': 2e-10},
        {'name': 'native bins versus continuum Hankel reference', 'passed': final['max_bin_error_over_peak'] < 4e-4,
         'value': final['max_bin_error_over_peak'], 'threshold': 4e-4},
        {'name': 'native analytic input gradient versus continuum reference', 'passed': final['max_gradient_error_over_peak'] < 7e-4,
         'value': final['max_gradient_error_over_peak'], 'threshold': 7e-4},
        {'name': 'degree refinement', 'passed': degree_change < 5e-5,
         'value': degree_change, 'threshold': 5e-5},
        {'name': 'spacing refinement', 'passed': spacing_change < 2e-4,
         'value': spacing_change, 'threshold': 2e-4},
        {'name': 'box refinement', 'passed': box_change < 5e-5,
         'value': box_change, 'threshold': 5e-5},
        {'name': 'conservative observed bin error budget', 'passed': observed_error_budget < 4e-4,
         'value': observed_error_budget, 'threshold': 4e-4},
        {'name': 'input derivative finite difference', 'passed': derivative_residual < 2e-9,
         'value': derivative_residual, 'threshold': 2e-9},
    ]
    report = {
        'schema': 'k-space-canonical-continuum-benchmark/1',
        'passed': all(row['passed'] for row in tests),
        'scope': 'Physical nf=5 N3LL canonical evolution acting on a non-Gaussian input; strict-NLO boundaries are not used in this operator-only benchmark',
        'input': {'family': 'normalized Student radial density', 'scale_GeV': args.scale,
                  'power': args.power, 'analytic_hankel_transform': True,
                  'analytic_scale_derivative': True},
        'endpoint': {'symbol_zero': float(symbol(np.array([0.]))[0]),
                     'symbol_infinity': float(symbol(np.array([2*rg.s0]))[0]),
                     'contact_retained_analytically_in_reference': True},
        'bin_edges_GeV': edges.tolist(),
        'reference_bins': reference.tolist(),
        'reference_gradient_bins': gradient_reference.tolist(),
        'reference_quadrature_errors': [reference_error, gradient_reference_error],
        'refinements': rows,
        'changes': {'degree_512_to_1024': degree_change,
                    'spacing_0p025_to_0p0125': spacing_change,
                    'radius_640_to_1280': box_change},
        'observed_bin_error_budget_over_peak': observed_error_budget,
        'tests': tests,
        'sha256': {
            'script': digest(Path(__file__)),
            'operator': digest(ROOT/'src/resummed_operator.py'),
            'radial_evaluator': digest(ROOT/'src/rational_radial.py'),
            'coefficients': digest(ROOT/'data/coefficients_css2_nf5_n3ll.json'),
        },
        'limitations': {
            'physical_N3LLprime_prediction_validated': False,
            'full_process_cross_section': False,
            'rigorous_global_rational_supremum_certified': False,
        },
    }
    (args.out/'canonical_continuum_benchmark.json').write_text(
        json.dumps(report, indent=2, allow_nan=False)+'\n', encoding='utf-8')
    print('PASSED', report['passed'], flush=True)
    return 0 if report['passed'] else 1


if __name__ == '__main__':
    raise SystemExit(main())
