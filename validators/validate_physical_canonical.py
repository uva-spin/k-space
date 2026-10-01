#!/usr/bin/env python3
"""Exercise the sourced physical N3LL coefficients in the canonical operator.

This validates the evolution input and its implementation. The hard and OPE
boundaries remain strict NLO, so the report never claims an N3LL-prime spectrum.
The finite-grid checks use a non-Gaussian input and do not certify continuum
box or spacing errors.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import sys

import numpy as np
from scipy.integrate import quad
from scipy.linalg import eigh
from scipy.sparse import diags, eye, kron

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'src'))
from css2_coefficients import qcd_coefficients  # noqa: E402
from resummed_operator import (  # noqa: E402
    CanonicalRG, Coefficients, apply_chebyshev_resolvent,
    apply_resolvent_polynomial, bernstein_scalar,
)


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, default=ROOT/'outputs/physical_canonical')
    args = parser.parse_args(argv)
    args.out.mkdir(parents=True, exist_ok=True)

    tests = []
    def check(name, passed, residual=None, tolerance=None, detail=''):
        tests.append({'name': name, 'passed': bool(passed), 'residual': residual,
                      'tolerance': tolerance, 'detail': detail})

    coefficient_path = ROOT/'data/coefficients_css2_nf5_n3ll.json'
    coefficient = Coefficients.load(coefficient_path)
    independent = qcd_coefficients(5)
    lengths = tuple(map(len, (coefficient.beta, coefficient.cusp,
                              coefficient.noncusp, coefficient.d)))
    check('N3LL coefficient orders are 4/4/3/3', lengths == (4, 4, 3, 3),
          detail=str(lengths))
    check('accuracy metadata is N3LL evolution', coefficient.evolution_accuracy == 'N3LL')
    check('boundary metadata is strict NLO', coefficient.hard_and_ope_order == 'strict NLO')
    check('payload does not claim N3LL-prime', not coefficient.physical_n3llprime_payload)
    check('fixed-flavor metadata is nf=5', coefficient.nf == 5)
    for name in ('beta', 'cusp', 'noncusp', 'd'):
        residual = float(np.max(np.abs(np.asarray(getattr(coefficient, name))
                                       - np.asarray(getattr(independent, name)))))
        check(f'JSON and formula implementation agree: {name}', residual < 2e-12,
              residual=residual, tolerance=2e-12)

    source_manifest = ROOT/(coefficient.source_manifest or '')
    source_data = json.loads(source_manifest.read_text()) if source_manifest.is_file() else {}
    archive_records = [row for row in source_data.get('records', []) if 'archive' in row]
    check('source manifest is present and populated', len(archive_records) >= 5,
          detail=f'{len(archive_records)} archive records')
    check('source records carry URL, digest, and locator', all(
        row.get('url', '').startswith('https://arxiv.org/') and
        len(row.get('sha256', '')) == 64 and row.get('locator_and_conversion')
        for row in archive_records))

    rg = CanonicalRG(coefficient)
    s = np.geomspace(2e-4, 3., 31)
    initial, middle, final = (3., 16.), (9., 100.), (31., 1600.)
    composition = rg.log_u(s, final, initial) - (
        rg.log_u(s, final, middle) + rg.log_u(s, middle, initial))
    residual = float(np.max(np.abs(composition)))
    check('two-scale evolution composes', residual < 3e-11,
          residual=residual, tolerance=3e-11)

    ti, wi = rg.endpoint(*initial)
    tf, wf = rg.endpoint(*final)
    def virtuality_integrand(t):
        a = rg.values(t)[0]
        return coefficient.series(a, coefficient.noncusp) - 2*(wi-t)*coefficient.series(a, coefficient.cusp)
    direct = quad(virtuality_integrand, ti, tf, epsabs=2e-12, epsrel=2e-12)[0]
    direct = direct - (wf-wi)*rg.D(s, final[0])
    residual = float(np.max(np.abs(direct-rg.log_u(s, final, initial))))
    check('independent quadrature agrees with closed RG integrals', residual < 2e-10,
          residual=residual, tolerance=2e-10)

    h = 1e-5
    a_final = rg.values(tf)[0]
    expected_mu = (coefficient.series(a_final, coefficient.noncusp)
                   - coefficient.series(a_final, coefficient.cusp)
                   * np.log(final[1]/final[0]**2))
    finite_mu = (rg.log_u(s, (final[0]*np.exp(h), final[1]), initial)
                 - rg.log_u(s, (final[0]*np.exp(-h), final[1]), initial))/(2*h)
    residual = float(np.max(np.abs(finite_mu-expected_mu)))
    check('virtuality endpoint equation', residual < 3e-8,
          residual=residual, tolerance=3e-8)
    finite_zeta = (rg.log_u(s, (final[0], final[1]*np.exp(2*h)), initial)
                   - rg.log_u(s, (final[0], final[1]*np.exp(-2*h)), initial))/(2*h)
    residual = float(np.max(np.abs(finite_zeta+rg.D(s, final[0]))))
    check('rapidity endpoint equation', residual < 3e-8,
          residual=residual, tolerance=3e-8)

    # A positive, radial, algebraically tailed input avoids a Gaussian-only test.
    n, halfwidth = 19, 10.
    spacing = 2*halfwidth/(n+1)
    one_d = diags([-np.ones(n-1), 2*np.ones(n), -np.ones(n-1)], [-1, 0, 1],
                  format='csc')/spacing**2
    laplacian = kron(one_d, eye(n, format='csc')) + kron(eye(n, format='csc'), one_d)
    axis = np.linspace(-halfwidth+spacing, halfwidth-spacing, n)
    x, y = np.meshgrid(axis, axis, indexing='ij')
    input_array = (1+x*x+y*y)**-3
    input_array = input_array.ravel()
    input_array /= np.linalg.norm(input_array)*spacing
    eigenvalues, eigenvectors = eigh(laplacian.toarray())
    exact_symbol = rg.v(eigenvalues, 20., 400.)
    reference = eigenvectors@(exact_symbol*(eigenvectors.T@input_array))
    degree_rows = []
    for degree in (16, 64, 256):
        minimum, nodes = rg.h_nodes(degree, rg.s0, 20., 400.)
        approximation = apply_resolvent_polynomial(
            laplacian, input_array, rg.s0, minimum, nodes)
        r = rg.s0/(rg.s0+eigenvalues)
        polynomial = minimum+r*bernstein_scalar(nodes, r)
        oracle = eigenvectors@(polynomial*(eigenvectors.T@input_array))
        recurrence_residual = float(np.linalg.norm(approximation-oracle)*spacing)
        same_grid_error = float(np.linalg.norm(approximation-reference)*spacing)
        bound = float(np.max(np.abs(polynomial-exact_symbol)))
        degree_rows.append({'degree': degree, 'same_grid_L2_error': same_grid_error,
                            'finite_matrix_bound': bound,
                            'recurrence_oracle_residual': recurrence_residual})
        check(f'non-Gaussian recurrence oracle N={degree}', recurrence_residual < 3e-10,
              residual=recurrence_residual, tolerance=3e-10)
        check(f'non-Gaussian finite-matrix bound N={degree}', same_grid_error <= bound+3e-10,
              residual=same_grid_error, tolerance=bound+3e-10)
    check('non-Gaussian degree refinement improves',
          degree_rows[-1]['same_grid_L2_error'] < degree_rows[0]['same_grid_L2_error'],
          detail=f"{degree_rows[0]['same_grid_L2_error']} -> {degree_rows[-1]['same_grid_L2_error']}")

    chebyshev = rg.chebyshev_coefficients(128, rg.s0, 20., 400.)
    chebyshev_result = apply_chebyshev_resolvent(
        laplacian, input_array, rg.s0, chebyshev)
    chebyshev_symbol = np.polynomial.chebyshev.chebval(
        2*rg.s0/(rg.s0+eigenvalues)-1, chebyshev)
    chebyshev_oracle = eigenvectors@(chebyshev_symbol*(eigenvectors.T@input_array))
    residual = float(np.linalg.norm(chebyshev_result-chebyshev_oracle)*spacing)
    check('Chebyshev resolvent agrees with matrix-polynomial oracle', residual < 3e-10,
          residual=residual, tolerance=3e-10)
    endpoint_residual = float(max(
        abs(np.polynomial.chebyshev.chebval(1., chebyshev)),
        abs(np.polynomial.chebyshev.chebval(-1., chebyshev)
            - rg.v(np.array([2*rg.s0]), 20., 400.)[0])))
    check('Chebyshev interpolant preserves both spectral endpoints', endpoint_residual < 2e-13,
          residual=endpoint_residual, tolerance=2e-13)

    demonstration = CanonicalRG(Coefficients.load(ROOT/'data/coefficients_demo.json'))
    physical_values = rg.v(s, 20., 400.)
    demo_values = demonstration.v(s, 20., 400.)
    residual = float(np.max(np.abs(physical_values-demo_values)))
    check('physical coefficient input changes the numerical symbol', residual > 1e-6,
          residual=residual, tolerance='>1e-6')

    report = {
        'schema': 'k-space-physical-canonical-validation/1',
        'passed': all(item['passed'] for item in tests),
        'tests': tests,
        'coefficient_orders': {'beta': 4, 'cusp': 4, 'noncusp': 3, 'd': 3},
        'coefficient_metadata': {
            'description': coefficient.description,
            'scheme': coefficient.scheme,
            'evolution_accuracy': coefficient.evolution_accuracy,
            'hard_and_ope_order': coefficient.hard_and_ope_order,
            'nf': coefficient.nf,
        },
        'degree_refinement': degree_rows,
        'chebyshev_degree': 128,
        'input': 'positive radial algebraic tail (1+k_T^2)^-3 on a finite Dirichlet grid',
        'sha256': {
            'coefficient_input': digest(coefficient_path),
            'source_manifest': digest(source_manifest),
            'operator': digest(ROOT/'src/resummed_operator.py'),
            'coefficient_formula': digest(ROOT/'src/css2_coefficients.py'),
        },
        'scope': {
            'physical_N3LL_evolution_input_exercised': True,
            'hard_and_ope_boundaries': 'strict NLO',
            'physical_N3LLprime_prediction_validated': False,
            'continuum_grid_error_certified': False,
            'data_fit_performed': False,
        },
    }
    (args.out/'physical_canonical_validation.json').write_text(
        json.dumps(report, indent=2, allow_nan=False)+'\n', encoding='utf-8')
    print(f"physical canonical validation: {sum(t['passed'] for t in tests)}/{len(tests)}")
    return 0 if report['passed'] else 1


if __name__ == '__main__':
    raise SystemExit(main())
