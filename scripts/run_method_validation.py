#!/usr/bin/env python3
"""Portable reproduction of the methods-paper validation, with explicit PDF paths."""
from __future__ import annotations
import argparse
import json
import os
from pathlib import Path
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[1]

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--physical', action='store_true', help='Also run the three-scale actual-PDF continuum audit')
    parser.add_argument('--library', type=Path, help='Built libkt_lhapdf shared library')
    parser.add_argument('--pdf-data', type=Path, help='Directory containing the pinned NNPDF set')
    parser.add_argument('--out', type=Path, default=ROOT/'outputs/method_validation')
    args = parser.parse_args()
    if args.physical and (args.library is None or args.pdf_data is None):
        parser.error('--physical requires --library and --pdf-data')
    out = args.out.resolve()
    out.mkdir(parents=True, exist_ok=True)
    env = dict(os.environ, OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1', MKL_NUM_THREADS='1')
    record = {'schema': 'method-paper-portable-run/1', 'passed': False, 'physical_requested': args.physical, 'steps': []}
    def save():
        (out/'reproduction.json').write_text(json.dumps(record, indent=2)+'\n')
    def run(name, argv, allowed=(0,)):
        print('STAGE', name, flush=True)
        start = time.monotonic()
        with (out/(name+'.log')).open('w') as log:
            result = subprocess.run([sys.executable, *map(str,argv)], cwd=ROOT,
                                     env=env, stdout=log, stderr=subprocess.STDOUT)
        record['steps'].append({'name': name, 'arguments': list(map(str,argv)),
            'return_code': result.returncode, 'allowed_exit_codes': list(allowed),
            'elapsed_seconds': time.monotonic()-start})
        save()
        if result.returncode not in allowed:
            raise RuntimeError(f'{name} failed; see {out/(name+".log")}')
    run('policy', ['tools/check_repository_policy.py', '--check-prepared-hashes'])
    run('documentation', ['tools/check_documentation.py'])
    run('snapshot', ['tools/check_validation_snapshot.py'])
    run('software', ['-m', 'pytest', '-q'])
    run('suites', ['run_all_checks.py', '--output-dir', 'outputs/method_suites'])
    run('continuum', ['scripts/benchmark_canonical_continuum.py', '--out', out/'continuum'])
    if args.physical:
        pdf = ['--library', args.library.resolve(), '--pdf-data', args.pdf_data.resolve()]
        raw = out/'finite_table/canonical_physical_audit.json'
        run('physical_finite_table', ['scripts/benchmark_canonical_physical.py', *pdf,
                                      '--extended', '--out', raw.parent], allowed=(0, 1))
        # Exit 1 is not accepted on its own. The next auditor rejects all raw
        # failures except the explicitly addressed auxiliary UV class target.
        complete = out/'completed/canonical_physical_audit.json'
        run('physical_completed', ['scripts/audit_canonical_uv_completion.py', *pdf,
                                   '--input', raw, '--out', complete.parent])
        complete = complete.parent/'completed_physical_audit.json'
        run('support', ['scripts/audit_canonical_support.py', '--physical-report', complete,
                        '--out', out/'support'])
    record['passed'] = True
    record['scope'] = 'Offline suites and continuum plus full canonical physical audit' if args.physical else 'Offline suites and ordinary-input continuum only; actual-PDF physical audit was not requested'
    save()
    print('PASS:', record['scope'], flush=True)

if __name__ == '__main__':
    main()
