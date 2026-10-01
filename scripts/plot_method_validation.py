#!/usr/bin/env python3
"""Regenerate scientific figures and tables from the completed-input audit."""
import argparse
import csv
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
p = argparse.ArgumentParser(description=__doc__)
p.add_argument('--reference-dir', type=Path, default=ROOT/'validation/reference/2026-10-01')
p.add_argument('--out', type=Path, default=ROOT/'outputs/method_figures')
args = p.parse_args()
report = json.loads((args.reference_dir/'canonical_physical_completed.json').read_text())
if not report['passed']:
    raise RuntimeError('Cannot publish figures for a failed audit')
figdir = args.out
tables = args.out/'tables'
figdir.mkdir(parents=True, exist_ok=True)
tables.mkdir(parents=True, exist_ok=True)
edges = np.asarray(report['bin_edges_GeV'])
centers = (edges[:-1]+edges[1:])/2
plt.rcParams.update({'font.size': 10, 'axes.spines.top': False,
                     'axes.spines.right': False, 'figure.dpi': 150})
fig, axes = plt.subplots(2, 3, figsize=(12.8, 7.), sharex='col', constrained_layout=True)
colors = ['#21618c', '#b03a2e', '#148f77', '#8e44ad', '#b9770e', '#5d6d7e']
for col, case in enumerate(report['cases']):
    peak = case['normalization_peak']
    final = case['refinements'][2]
    reference = np.asarray(case['reference_totals'])
    native = np.asarray(final['total_bins'])
    ax = axes[0, col]
    ax.axvspan(.2*case['Q'], edges[-1], color='0.95', zorder=-3)
    ax.stairs(reference/peak, edges, color='black', label='Continuum reference', linewidth=1.3)
    ax.plot(centers, native/peak, 'o', markersize=3.5, color='#148f77', label='Native resolvent')
    ax.axhline(0, color='0.65', linewidth=.7)
    ax.set_title(f"Q = {case['Q']:g} GeV\nxA = {case['xA']:g}, xB = {case['xB']:g}")
    ax.set_xscale('symlog', linthresh=1)
    ax.set_xlim(0, 80)
    ax.set_ylabel('Finite-bin W / peak reference bin' if col == 0 else '')
    if col == 0:
        ax.legend(loc='lower left', fontsize=8)
    err = (np.asarray(final['components'])-np.asarray(case['reference_components']))/peak
    ax = axes[1, col]
    for j, component in enumerate(report['components']):
        ax.plot(centers, err[:, j], '.-', color=colors[j], linewidth=.75, markersize=3,
                label=component.replace('ope_', '').replace('_', ' '))
    ax.plot(centers, err.sum(axis=1), 'k--', linewidth=1., label='total')
    ax.axhline(0, color='0.65', linewidth=.7)
    ax.set_xscale('symlog', linthresh=1)
    ax.ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
    ax.set_ylabel('(Native - reference) / peak' if col == 0 else '')
    ax.set_xlabel(r'$q_T$ [GeV]')
    ax.grid(alpha=.18)
axes[1, 1].legend(loc='upper left', fontsize=7, ncol=2)
fig.suptitle('Uncapped canonical CSS2 with physical PDFs and strict-NLO boundaries\n'
             'Shading marks qT > 0.2 Q: numerical stress region, not a full-spectrum prediction', fontsize=12)
fig.savefig(figdir/'canonical_physical_audit.png', dpi=180)
fig.savefig(figdir/'canonical_physical_audit.pdf', metadata={'CreationDate': None, 'ModDate': None})
plt.close(fig)

fig, ax = plt.subplots(figsize=(8., 4.5), constrained_layout=True)
labels = [f"Q={c['Q']:g}" for c in report['cases']]
bottom = np.zeros(3)
parts = [('bin discrepancy', lambda c: c['refinements'][2]['bin_error_over_peak']),
         ('degree change', lambda c: c['changes']['degree']),
         ('spacing change', lambda c: c['changes']['spacing']),
         ('box change', lambda c: c['changes']['box']),
         ('PDF UV class bound', lambda c: c['UV_class_envelope_over_peak'])]
for label, getter in parts:
    values = np.array([getter(c) for c in report['cases']])
    ax.bar(labels, values, bottom=bottom, label=label)
    bottom += values
ax.axhline(1e-3, color='black', linestyle='--', label='Preregistered budget limit')
ax.set_ylim(0, 1.1e-3)
ax.set_ylabel('Observed contributions / peak reference bin')
ax.set_title('Separated refinement budget, including the broader PDF continuation class')
ax.legend(fontsize=8, ncol=2)
fig.savefig(figdir/'canonical_error_budget.png', dpi=180)
fig.savefig(figdir/'canonical_error_budget.pdf', metadata={'CreationDate': None, 'ModDate': None})
plt.close(fig)

with (tables/'canonical_physical_bins.csv').open('w', newline='') as f:
    writer = csv.writer(f)
    writer.writerow(['Q_GeV', 'xA', 'xB', 'q_low', 'q_high', 'reference_total', 'native_total',
                     'reference_log_scale_gradient', 'native_log_scale_gradient'])
    for c in report['cases']:
        final = c['refinements'][2]
        grad = np.asarray(c['reference_component_gradients']).sum(axis=1)
        for i, (lo, hi) in enumerate(zip(edges[:-1], edges[1:])):
            writer.writerow([c['Q'], c['xA'], c['xB'], lo, hi, c['reference_totals'][i],
                             final['total_bins'][i], grad[i], final['gradient_bins'][i]])
with (tables/'canonical_physical_summary.csv').open('w', newline='') as f:
    writer = csv.writer(f)
    writer.writerow(['Q_GeV', 'bin_error_class_bound', 'component_error_class_bound',
                     'gradient_error_class_bound', 'contact_error_class_bound',
                     'observed_budget', 'UV_input_spread', 'UV_integration_residual'])
    for c in report['cases']:
        e = c['class_wide_error_bounds_over_peak']
        writer.writerow([c['Q'], e['bin'], e['component'], e['gradient'], e['contact'],
                         c['observed_budget'], c['UV_completion_spread_over_peak'],
                         max(m['UV_residual_bound_over_peak'] for m in c['UV_completion_models'])])
print('Wrote physical canonical figures and source tables.')


report = json.loads((args.reference_dir/'canonical_continuum.json').read_text())
edges = np.asarray(report['bin_edges_GeV'])
centers = (edges[:-1]+edges[1:])/2
reference = np.asarray(report['reference_bins'])
final = next(row for row in report['refinements']
             if row['degree'] == 1024 and row['radius_GeV'] == 1280
             and row['spacing_GeV'] == .0125)
native = np.asarray(final['bins'])
gradient_reference = np.asarray(report['reference_gradient_bins'])
gradient_native = np.asarray(final['gradient_bins'])
scale = np.max(np.abs(reference))

plt.rcParams.update({'font.size': 9, 'axes.labelsize': 10, 'legend.fontsize': 8})
fig, axes = plt.subplots(2, 2, figsize=(7.0, 5.3), constrained_layout=True)
ax = axes[0, 0]
ax.plot(centers, reference, 'o-', label='continuum Hankel reference', ms=4)
ax.plot(centers, native, 's--', label='native resolvent', ms=3)
ax.set(xlabel=r'$q_T$ bin center [GeV]', ylabel='bin action', title='Non-Gaussian canonical action')
ax.legend(frameon=False)

ax = axes[1, 0]
ax.axhline(0, color='0.7', lw=.8)
ax.plot(centers, (native-reference)/scale, 'o-', color='#b33b2e', ms=4)
ax.set(xlabel=r'$q_T$ bin center [GeV]', ylabel='difference / peak bin')
ax.ticklabel_format(axis='y', style='sci', scilimits=(-2, 2))

ax = axes[0, 1]
degree_rows = [row for row in report['refinements']
               if row['radius_GeV'] == 1280 and row['spacing_GeV'] == .0125]
ax.loglog([row['degree'] for row in degree_rows],
          [row['max_bin_error_over_peak'] for row in degree_rows], 'o-', label='bins')
ax.loglog([row['degree'] for row in degree_rows],
          [row['max_gradient_error_over_peak'] for row in degree_rows], 's--', label='gradient')
ax.axhline(4e-4, color='0.55', ls=':', label='bin threshold')
ax.set(xlabel='Chebyshev degree', ylabel='maximum error / peak bin', title='Degree refinement')
ax.legend(frameon=False)

ax = axes[1, 1]
ax.axhline(0, color='0.7', lw=.8)
ax.plot(centers, (gradient_native-gradient_reference)/scale, 'o-', color='#376aa3', ms=4)
ax.set(xlabel=r'$q_T$ bin center [GeV]', ylabel='gradient difference / peak bin',
       title=r'Analytic $\partial/\partial\lambda$ check')
ax.ticklabel_format(axis='y', style='sci', scilimits=(-2, 2))

for suffix in ('pdf', 'png'):
    fig.savefig(args.out/f'canonical_continuum.{suffix}', dpi=220)
plt.close(fig)

with (tables/'canonical_continuum_bins.csv').open('w', newline='') as stream:
    writer = csv.writer(stream)
    writer.writerow(['qT_low_GeV', 'qT_high_GeV', 'reference', 'native',
                     'difference_over_peak', 'gradient_reference', 'gradient_native'])
    for index in range(len(reference)):
        writer.writerow([edges[index], edges[index+1], reference[index], native[index],
                         (native[index]-reference[index])/scale,
                         gradient_reference[index], gradient_native[index]])

with (tables/'canonical_refinement.csv').open('w', newline='') as stream:
    fields = ['degree', 'radius_GeV', 'spacing_GeV', 'max_bin_error_over_peak',
              'max_gradient_error_over_peak']
    writer = csv.DictWriter(stream, fieldnames=fields)
    writer.writeheader()
    for row in report['refinements']:
        writer.writerow({key: row[key] for key in fields})

print('wrote canonical continuum figure and tables')
