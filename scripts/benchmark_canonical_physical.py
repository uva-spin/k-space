#!/usr/bin/env python3
"""Same-canonical-prescription physical strict-NLO continuum audit."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import sys
import time

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'src'))
from canonical_physical import CanonicalPhysicalSymbol, physical_rg, bin_chebyshev_moments, COMPONENTS
from canonical_independent_reference import IndependentCanonical, TransversePair, continuum_components, ultraviolet_envelope
from rational_radial import RadialGrid
from lhapdf_provider import LHAPDFProvider


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def polynomial_range(coefficients, upper):
    p = np.polynomial.Polynomial(coefficients)
    points = [0., upper]+[z.real for z in p.deriv().roots()
                           if abs(z.imag) < 1e-12 and 0 < z.real < upper]
    return float(min(p(points))), float(max(p(points)))


def ultraviolet_monotonicity(ref):
    # d(log v)/d(log mu_c) = d-gamma + 2 ln(Q/mu_c)*(Gamma+B*d'),
    # with B=sum beta_n a^(n+2). Prove signs on the whole interval [0,a_max]
    # by finding all real polynomial stationary points, not sampling a grid.
    from numpy.polynomial import Polynomial as P
    amax = float(ref.running(np.log(ref.pdf_max))[0, 0])
    d, g, cusp = [P([0.]+ref.c[k]) for k in ('d', 'noncusp', 'cusp')]
    B = P([0., 0.]+ref.c['beta'])
    dg = polynomial_range((d-g).coef[1:], amax)
    bracket = polynomial_range((cusp+B*d.deriv()).coef[1:], amax)
    return {'passed': bool(dg[1] <= 0 and bracket[0] > 0 and ref.pdf_max > ref.Q),
            'a_max': amax, '(d-gamma)/a_range': dg,
            '(Gamma+B*dprime)/a_range': bracket}


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--out', type=Path, required=True)
    p.add_argument('--library', type=Path, required=True)
    p.add_argument('--pdf-data', type=Path, required=True)
    p.add_argument('--quick', action='store_true')
    p.add_argument('--extended', action='store_true', help='Extended high-scale box/degree refinement after the preserved first-run failure')
    p.add_argument('--case', type=int, choices=[0, 1, 2])
    args = p.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    start = time.monotonic()
    pdf = LHAPDFProvider(args.library, args.pdf_data)
    pdf_manifest = json.loads((ROOT/'manifests/css2_pdf.json').read_text())
    for row in pdf_manifest['files']:
        if row['name'].endswith(('.info', '.dat')):
            if digest(args.pdf_data/pdf_manifest['set']/row['name']) != row['sha256']:
                raise RuntimeError('PDF source hash mismatch')
    rg = physical_rg()
    edges = np.array([0., .5, 1., 2., 3., 4., 6., 8., 12., 20., 40., 80.])
    settings = [(8., .08, .2, False), (20., .15, .3, True), (91.1876, .02, .07, False)]
    cases = []
    for index, (Q, xA, xB, gaussian) in enumerate(settings):
        if args.case is not None and index != args.case:
            continue
        print('CASE', index, Q, xA, xB, 'building independent longitudinal/RG paths', flush=True)
        case_start = time.monotonic()
        production = CanonicalPhysicalSymbol(pdf, rg, Q, xA, xB)
        reference = IndependentCanonical(pdf, ROOT/'data/coefficients_css2_nf5_n3ll.json', Q, xA, xB)
        pair = TransversePair(.9+.8*xA, 1.2+.6*xB, gaussian=gaussian)
        print('CASE', index, 'reference bins', flush=True)
        ref_bins, quad_error = continuum_components(reference, pair, edges)
        totals = ref_bins[:, :, 0].sum(axis=1)
        peak = float(np.max(abs(totals)))
        uv = ultraviolet_envelope(reference, edges)
        monotonicity = ultraviolet_monotonicity(reference)
        # Compare independently assembled symbols through transition, canonical
        # ultraviolet region (including mu_c>Q), and constant infrared region.
        b = np.unique(np.r_[np.geomspace(2e-5, 1., 49),
                            np.sqrt(rg.s0)*np.array([.7, .707106781, .8, .9, .9999, 1.])])
        sp, sr = production.components(b*b), reference.components(b)
        symbol_error = float(np.max(abs(sp-sr))/np.max(abs(sr.sum(axis=1))))
        adaptive_errors = []
        cache_errors = []
        rg_errors = []
        for mu in [5., 6.4, 13., 35., 150., 1900., 28000., 1e5]:
            rb = reference.direct_blocks(mu, adaptive=True)
            cached = reference.cache(np.log(mu))
            scale = max(np.max(abs(rb)), 1e-14)
            adaptive_errors.append(float(np.max(abs(cached-rb))/scale))
            direct = production.direct_blocks(mu, order=512)
            cached_p = production.cache(np.log(mu/rg.mu0))
            cache_errors.append(float(np.max(abs(cached_p-direct))/max(np.max(abs(direct)), 1e-14)))
            # The reference leg is integrated from Q; production uses G,N,J.
            a, integral = reference.running(np.log(mu))[:, 0]
            log_reference = -integral-reference.series(a, reference.c['d'])*np.log(Q/mu)
            log_production = float(rg.log_v_tc(np.log(mu/rg.mu0), Q, Q*Q))
            rg_errors.append(abs(log_reference-log_production))
        if args.quick:
            configs = [(256, 320., .025)]
        else:
            configs = [(512, 1280., .0125), (1024, 1280., .0125),
                       (2048, 1280., .0125), (2048, 1280., .025),
                       (2048, 320., .0125), (2048, 640., .0125)]
            if args.extended and index == 2:
                configs = [(1024, 5120., .0125), (2048, 5120., .0125),
                           (4096, 5120., .0125), (4096, 5120., .025),
                           (4096, 1280., .0125), (4096, 2560., .0125)]
        rows = []
        polynomial_cache = {}
        for degree, radius, spacing in configs:
            grid = RadialGrid(radius, round(radius/spacing))
            native = pair.native(grid.k)
            delta = np.zeros(grid.cells)
            delta[0] = 1/grid.weights[0]
            native = np.column_stack((native, delta))
            moments = bin_chebyshev_moments(grid, native, edges, degree, rg.s0)
            if degree not in polynomial_cache:
                polynomial_cache[degree] = production.polynomials(degree, rg.s0)
            coeff = polynomial_cache[degree]
            bins = np.einsum('nbi,nc->bci', moments, coeff)
            prediction, gradient = bins[:, :, 0].sum(axis=1), bins[:, :, 1].sum(axis=1)
            delta = prediction-totals
            row = {'degree': degree, 'radius': radius, 'spacing': spacing,
                   'components': bins[:, :, 0].tolist(), 'component_gradients': bins[:, :, 1].tolist(),
                   'total_bins': prediction.tolist(), 'gradient_bins': gradient.tolist(),
                   'bin_error_over_peak': float(np.max(abs(delta))/peak),
                   'component_error_over_peak': float(np.max(abs(bins[:, :, 0]-ref_bins[:, :, 0]))/peak),
                   'gradient_error_over_peak': float(np.max(abs(gradient-ref_bins[:, :, 1].sum(axis=1)))/peak),
                   'contact_error_over_peak': float(np.max(abs(bins[:, :, 2]-ref_bins[:, :, 2]))/peak),
                   'bin_relative_errors': [float(d/t) if abs(t) > 1e-8*peak else None for d, t in zip(delta, totals)],
                   'negative_bins_retained': int(np.sum(prediction < 0))}
            rows.append(row)
            print('REFINE', index, degree, radius, spacing, row['bin_error_over_peak'],
                  row['component_error_over_peak'], row['gradient_error_over_peak'], flush=True)
            partial = {'status': 'in_progress', 'case': index, 'refinements': rows}
            (args.out/f'case_{index}_partial.json').write_text(json.dumps(partial, indent=2)+'\n')
        final = rows[0] if args.quick else rows[2]
        changes = {}
        if not args.quick:
            for name, earlier in [('degree', rows[1]), ('spacing', rows[3]), ('box', rows[5])]:
                changes[name] = float(np.max(abs(np.array(final['total_bins'])-earlier['total_bins']))/peak)
        # Direct-momentum convolution quadrature and its analytic derivative.
        ks = np.linspace(0., 90., 721)
        coarse_input, fine_input = pair.native(ks), pair.native(ks, 320)
        convolution_error = float(np.max(abs(coarse_input-fine_input))/np.max(abs(fine_input[:, 0])))
        step = 1e-5
        plus = TransversePair(pair.a*np.exp(step), pair.b, gaussian=gaussian).native(ks)[:, 0]
        minus = TransversePair(pair.a*np.exp(-step), pair.b, gaussian=gaussian).native(ks)[:, 0]
        derivative_error = float(np.max(abs((plus-minus)/(2*step)-coarse_input[:, 1]))/np.max(abs(fine_input[:, 0])))
        uv_bin_error = float(np.max(uv.sum(axis=1))/peak)
        budget = final['bin_error_over_peak']+sum(changes.values())+uv_bin_error+quad_error/peak
        tests = {
            'independent_complete_symbol': symbol_error < 2e-5,
            'independent_adaptive_longitudinal_reference': max(adaptive_errors) < 3e-6,
            'production_longitudinal_cache_refinement': max(cache_errors) < 3e-6,
            'independent_direct_RG': max(rg_errors) < 2e-8,
            'native_bins': final['bin_error_over_peak'] < 5e-4,
            'each_strict_NLO_component': final['component_error_over_peak'] < 5e-4,
            'native_input_gradient': final['gradient_error_over_peak'] < 1e-3,
            'explicit_contact_action': final['contact_error_over_peak'] < 5e-4,
            'continuum_quadrature': quad_error/peak < 1e-7,
            'PDF_table_UV_envelope': uv_bin_error < 1e-6,
            'proved_UV_evolution_monotonicity': monotonicity['passed'],
            'momentum_convolution_refinement': convolution_error < 2e-7,
            'analytic_input_derivative': derivative_error < 2e-8,
            'observed_error_budget': budget < 1e-3,
            'full_refinement_study': not args.quick,
        }
        report = {'case': index, 'Q': Q, 'xA': xA, 'xB': xB,
                  'sqrt_s_GeV': float(Q/np.sqrt(xA*xB)),
                  'photon_prefactor_GeVminus4': float(4*np.pi*(1/137.035999084)**2/(9*Q**2*(Q**2/(xA*xB)))),
                  'input': {'family': 'Gaussian pair' if gaussian else 'Student pair',
                            'scale_A': pair.a, 'scale_B': pair.b, 'p': pair.p, 'q': pair.q},
                  'reference_components': ref_bins[:, :, 0].tolist(),
                  'reference_component_gradients': ref_bins[:, :, 1].tolist(),
                  'reference_contact_components': ref_bins[:, :, 2].tolist(),
                  'reference_totals': totals.tolist(), 'normalization_peak': peak,
                  'reference_quadrature_error': quad_error,
                  'symbol_error': symbol_error, 'adaptive_reference_errors': adaptive_errors,
                  'production_cache_errors': cache_errors, 'RG_errors': rg_errors,
                  'UV_bound_components': uv.tolist(), 'UV_bound_over_peak': uv_bin_error,
                  'UV_monotonicity': monotonicity, 'momentum_convolution_error': convolution_error,
                  'input_derivative_error': derivative_error,
                  'refinements': rows, 'changes': changes, 'observed_budget': budget,
                  'tests': tests, 'passed': all(tests.values()),
                  'elapsed_seconds': time.monotonic()-case_start}
        cases.append(report)
        (args.out/f'case_{index}.json').write_text(json.dumps(report, indent=2, allow_nan=False, default=lambda v: v.item())+'\n')
        print('CASE_DONE', index, report['passed'], {k:v for k,v in tests.items() if not v}, flush=True)
    result = {'schema': 'canonical-physical-audit/1', 'cases': cases,
              'passed': all(c['passed'] for c in cases), 'quick': args.quick, 'extended': args.extended,
              'components': COMPONENTS, 'bin_edges_GeV': edges.tolist(),
              'scope': 'uncapped canonical fixed-nf5 photon W structure; strict-NLO boundaries; no massive thresholds or Y term',
              'UV_assumption': 'beyond PDF Qmax each PDF is bounded by its absolute endpoint value times (1+ln(mu/Qmax))^4; explicit conditional bin bound, no production extrapolation',
              'reference_independence': 'separate endpoint RK45 RG, log-xi matching integration and flavor assembly, Bessel-K product input, adaptive continuum quadrature; only source coefficients and PDF provider shared',
              'elapsed_seconds': time.monotonic()-start,
              'source_sha256': {str(f.relative_to(ROOT)): digest(f) for f in [
                  ROOT/'src/canonical_physical.py', ROOT/'src/resummed_operator.py',
                  ROOT/'src/css2_matching.py', ROOT/'src/css2_coefficients.py',
                  ROOT/'src/lhapdf_provider.py', ROOT/'src/rational_radial.py',
                  ROOT/'scripts/canonical_independent_reference.py', Path(__file__),
                  ROOT/'data/coefficients_css2_nf5_n3ll.json',
                  ROOT/'manifests/css2_pdf.json']}}
    (args.out/'canonical_physical_audit.json').write_text(json.dumps(result, indent=2, allow_nan=False, default=lambda v: v.item())+'\n')
    print('PASSED', result['passed'], flush=True)
    return 0 if result['passed'] else 1


if __name__ == '__main__':
    raise SystemExit(main())
