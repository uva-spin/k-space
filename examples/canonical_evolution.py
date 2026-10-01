#!/usr/bin/env python3
"""Tabulate the unexpanded canonical RG symbol with the supplied coefficient input.

Run from the checkout root:
    python examples/canonical_evolution.py --mu 20 --zeta 400

The scalar spectral resolution s has units GeV^-2. This is not a fitted TMD or
cross section. At s=0 the canonical leg has its analytic zero endpoint; D and
mu_c are not evaluated there. Output files are generated under outputs/.
"""
from __future__ import annotations
import argparse
import csv
from pathlib import Path
import numpy as np
from _common import ROOT, add_rg_arguments, build_rg, finite_float, output_directory, provenance, write_json


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    add_rg_arguments(parser)
    parser.add_argument('--s', nargs='+', type=finite_float, default=[0.0, 0.01, 0.1, 1.0],
                        help='Nonnegative spectral resolutions in GeV^-2')
    parser.add_argument('--out', type=Path, default=ROOT/'outputs/examples/canonical')
    args = parser.parse_args(argv)
    try:
        if any(value < 0 for value in args.s): raise ValueError('s must be nonnegative')
        rg = build_rg(args)
        points = np.asarray(args.s)
        values = rg.v(points, args.mu, args.zeta)
        if not np.isfinite(values).all(): raise ValueError('Nonfinite symbol at the requested scales')
        rows = []
        for point, value in zip(points, values):
            if point == 0:
                tc = mu_c = kernel = None
            else:
                tc = float(rg.tc(point))
                mu_c = float(args.mu0*np.exp(tc))
                kernel = float(rg.D(point, args.mu))
                if not np.isfinite([tc, mu_c, kernel]).all():
                    raise ValueError('Nonfinite value outside the example domain')
            rows.append({'s_GeV_minus2': float(point), 'log_mu_c_over_mu0': tc,
                         'mu_c_GeV': mu_c, 'D_c': kernel, 'canonical_leg': float(value)})
        out = output_directory(args.out)
        with (out/'canonical_evolution.csv').open('w', newline='', encoding='utf-8') as stream:
            writer = csv.DictWriter(stream, fieldnames=list(rows[0])); writer.writeheader(); writer.writerows(rows)
        report = provenance(args)
        report.update({'description': rg.c.description, 'samples': rows,
                       'endpoint_convention': 'v(0)=0 analytically; no finite mu_c or D assigned at zero'})
        write_json(out/'canonical_evolution.json', report)
    except (ValueError, OSError, RuntimeError) as err:
        parser.error(str(err))
    label = rg.c.evolution_accuracy or 'lower-order demonstration'
    print(f'Canonical evolution symbol — {label}; boundaries: {rg.c.hard_and_ope_order or "demonstration only"}')
    print('s [GeV^-2]    v(s; mu, zeta)')
    for row in rows: print(f"{row['s_GeV_minus2']:11.5g}    {row['canonical_leg']:.12g}")
    print(f'Wrote: {out}/canonical_evolution.csv and .json')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
