#!/usr/bin/env python3
"""Apply canonical evolution to a Gaussian on a finite transverse-momentum grid.

Run from the checkout root:
    python examples/momentum_grid.py --reference

The evaluation path uses sparse elliptic solves. Optional independent
same-grid diagonalization supplies a validation reference only. It does not
estimate box/spacing error, and the example is not a physical TMD extraction.
"""
from __future__ import annotations
import argparse
import csv
from pathlib import Path
import numpy as np
from scipy.sparse import diags, eye, kron
from _common import ROOT, add_rg_arguments, build_rg, output_directory, provenance, write_json, finite_float
from resummed_operator import apply_resolvent_polynomial, bernstein_scalar


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    add_rg_arguments(parser)
    parser.add_argument('--grid', type=int, default=17, help='Interior nodes on each axis (3 through 65)')
    parser.add_argument('--halfwidth', type=finite_float, default=8., help='Box half-width in GeV')
    parser.add_argument('--degrees', nargs='+', type=int, default=[8, 32, 128], help='Bernstein degrees (1 through 512)')
    parser.add_argument('--reference', action='store_true', help='Compute independent same-grid reference (grid<=25)')
    parser.add_argument('--out', type=Path, default=ROOT/'outputs/examples/grid')
    args = parser.parse_args(argv)
    try:
        if not 3 <= args.grid <= 65: raise ValueError('grid must be between 3 and 65')
        if args.reference and args.grid > 25: raise ValueError('Dense same-grid reference is limited to grid<=25')
        if args.halfwidth <= 0: raise ValueError('halfwidth must be positive')
        if any(n < 1 or n > 512 for n in args.degrees): raise ValueError('Each degree must be between 1 and 512')
        degrees = sorted(set(args.degrees))
        rg = build_rg(args)
        n = args.grid; dk = 2*args.halfwidth/(n+1)
        T1 = diags([-np.ones(n-1), 2*np.ones(n), -np.ones(n-1)], [-1,0,1], format='csc')/dk**2
        T = kron(T1, eye(n, format='csc')) + kron(eye(n, format='csc'), T1)
        axis = np.linspace(-args.halfwidth+dk, args.halfwidth-dk, n)
        X, Y = np.meshgrid(axis, axis, indexing='ij')
        initial = np.exp(-(X*X+Y*Y)/2).ravel()
        initial /= np.linalg.norm(initial)*dk
        eigenvalues = vectors = reference = exact_symbol = None
        if args.reference:
            from scipy.linalg import eigh
            eigenvalues, vectors = eigh(T.toarray())
            exact_symbol = rg.v(eigenvalues, args.mu, args.zeta)
            reference = vectors@(exact_symbol*(vectors.T@initial))
        rows = []
        for degree in degrees:
            minf, nodes = rg.h_nodes(degree, rg.s0, args.mu, args.zeta)
            evolved = apply_resolvent_polynomial(T, initial, rg.s0, minf, nodes)
            if not np.isfinite(evolved).all(): raise ValueError('Nonfinite evolved array')
            row = {'degree': degree, 'input_L2': float(np.linalg.norm(initial)*dk),
                   'output_L2': float(np.linalg.norm(evolved)*dk), 'L2_same_grid_error': None,
                   'finite_matrix_bound': None, 'recurrence_oracle_residual': None}
            if args.reference:
                r = rg.s0/(rg.s0+eigenvalues)
                polynomial = minf+r*bernstein_scalar(nodes, r)
                row['L2_same_grid_error'] = float(np.linalg.norm(evolved-reference)*dk)
                row['finite_matrix_bound'] = float(np.max(np.abs(polynomial-exact_symbol)))
                oracle = vectors@(polynomial*(vectors.T@initial))
                row['recurrence_oracle_residual'] = float(np.linalg.norm(evolved-oracle)*dk)
                if row['recurrence_oracle_residual'] >= 2e-10:
                    raise RuntimeError('Recurrence disagrees with its independent matrix-polynomial oracle')
                if row['L2_same_grid_error'] > row['finite_matrix_bound']+2e-10:
                    raise RuntimeError('Same-grid error exceeds its finite-matrix bound')
            rows.append(row)
        out = output_directory(args.out)
        with (out/'degree_convergence.csv').open('w', newline='', encoding='utf-8') as stream:
            writer=csv.DictWriter(stream,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
        with (out/'momentum_grid.csv').open('w', newline='', encoding='utf-8') as stream:
            writer=csv.writer(stream);writer.writerow(['kx_GeV','ky_GeV','input_L2_normalized','evolved_largest_degree'])
            writer.writerows(zip(X.ravel(),Y.ravel(),initial,evolved))
        report=provenance(args)
        report.update({'description': rg.c.description, 'grid': {'nodes_per_axis':n,'halfwidth_GeV':args.halfwidth,
            'spacing_GeV':dk,'boundary':'Dirichlet','input':'Gaussian exp(-k^2/2), normalized in discrete L2'},
            'degrees':degrees, 'tau_GeV_minus2':rg.s0,'reference_computed':args.reference,
            'reference_method':'independent diagonalization of the same finite matrix' if args.reference else None,
            'results':rows, 'array_file':'momentum_grid.csv'})
        write_json(out/'momentum_grid.json',report)
    except (ValueError, OSError, RuntimeError) as err:
        parser.error(str(err))
    print('Momentum-space resolvent example — lower-order input, finite grid')
    for row in rows:
        text=f"N={row['degree']:3d}  output L2={row['output_L2']:.8g}"
        if args.reference: text+=f"  same-grid error={row['L2_same_grid_error']:.6g}  bound={row['finite_matrix_bound']:.6g}"
        print(text)
    print(f'Wrote: {out}/degree_convergence.csv, momentum_grid.csv and momentum_grid.json')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
