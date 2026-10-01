#!/usr/bin/env python3
"""Direct RG tolerance refinement, finite-bin conventions and failure controls."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

import numpy as np
from scipy.integrate import solve_ivp

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'src'))
from canonical_physical import physical_rg


def direct_log_leg(c, Q, mu, tolerance):
    series = lambda a, key, start=1: sum(v*a**(i+start) for i, v in enumerate(c[key]))
    beta = lambda a: -2*series(a, 'beta', 2)
    anchor = solve_ivp(lambda t, y: [beta(y[0])], (np.log(91.1876), np.log(Q)),
                       [.1180024/(4*np.pi)], method='RK45', rtol=tolerance,
                       atol=tolerance*.01)
    def rhs(t, y):
        a = y[0]
        return [beta(a), series(a, 'noncusp')-2*(np.log(Q)-t)*series(a, 'cusp')]
    answer = solve_ivp(rhs, (np.log(Q), np.log(mu)), [anchor.y[0, -1], 0.],
                       method='RK45', rtol=tolerance, atol=tolerance*.01)
    if not (anchor.success and answer.success):
        raise RuntimeError('Independent RG solve failed')
    a, integral = answer.y[:, -1]
    return float(-integral-series(a, 'd')*np.log(Q/mu))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--physical-report', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    c = json.loads((ROOT/'data/coefficients_css2_nf5_n3ll.json').read_text())
    rg = physical_rg()
    rows = []
    for Q in (8., 20., 91.1876):
        for mu in (5., 6.4, 13., 35., 150., 1900., 28000., 1e5):
            loose = direct_log_leg(c, Q, mu, 2e-9)
            tight = direct_log_leg(c, Q, mu, 2e-12)
            prod = float(rg.log_v_tc(np.log(mu/rg.mu0), Q, Q*Q))
            rows.append({'Q': Q, 'mu_c': mu, 'refinement': abs(tight-loose),
                         'production_tight_reference': abs(prod-tight)})
    report = json.loads(args.physical_report.read_text())
    negative = []
    units = []
    for case in report['cases']:
        comp = np.asarray(case['reference_components'])
        peak = case['normalization_peak']
        # Omitting one matching leg or inserting another factor of pi must be
        # visible relative to the declared numerical acceptance threshold.
        omit_matching = float(np.max(abs(comp[:, 2]))/peak)
        wrong_measure = float(np.max(abs((np.pi-1)*comp.sum(axis=1)))/peak)
        negative.append({'case': case['case'], 'omitted_qq_leg_A_error_over_peak': omit_matching,
                         'extra_pi_error_over_peak': wrong_measure,
                         'both_detected': omit_matching > 5e-4 and wrong_measure > 5e-4})
        # Source d/dqT^2 versus d/d^2qT differ by pi. Check the annulus factor
        # against the scalar source prefactor, then convert GeV^-2 -> pb.
        Q, xA, xB = case['Q'], case['xA'], case['xB']
        energy_squared = Q*Q/(xA*xB)
        alpha = 1/137.035999084
        prefactor_dqt2 = 4*np.pi**2*alpha**2/(9*Q**2*energy_squared)
        prefactor_d2qt = case['photon_prefactor_GeVminus4']
        units.append({'case': case['case'], 'prefactor_identity_residual': abs(prefactor_dqt2/np.pi-prefactor_d2qt),
                      'binned_dsigma_dQ2_dy_pb_per_GeV2': (prefactor_d2qt*3.893793721e8*comp.sum(axis=1)).tolist()})
    tests = {
        'RG_tolerance_refinement': max(x['refinement'] for x in rows) < 2e-8,
        'production_vs_tight_independent_RG': max(x['production_tight_reference'] for x in rows) < 2e-8,
        'physical_negative_controls': all(x['both_detected'] for x in negative),
        'annulus_normalization': max(x['prefactor_identity_residual'] for x in units) < 1e-22,
        'physical_cases_all_pass': report['passed'],
    }
    result = {'schema': 'canonical-support-audit/1', 'passed': all(tests.values()),
              'tests': tests, 'RG_refinement': rows, 'negative_controls': negative,
              'physical_units': units,
              'units_scope': 'photon-only d sigma/(dQ^2 dy), integrated over full-azimuth qT annulus; source conventions and leading-power xA xB s=Q^2',
              'UV_gradient_bound': 'A conservative log-scale-derivative tail bound is four times the reported ordinary-input UV bound for these normalized Student/Gaussian fixtures'}
    (args.out/'support_audit.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({'passed': result['passed'], 'tests': tests,
                      'max_RG_refinement': max(x['refinement'] for x in rows)}, indent=2))
    return 0 if result['passed'] else 1


if __name__ == '__main__':
    raise SystemExit(main())
