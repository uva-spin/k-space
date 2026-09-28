#!/usr/bin/env python3
"""Run manifest-listed, manuscript-independent verification suites.

Fresh per-suite directories, explicit JSON success fields, captured subprocess
logs, and input hashes prevent stale outputs or process exit alone from being
interpreted as successful scientific checks. No network is used by the runner.
"""
from __future__ import annotations
import argparse
import datetime as dt
import hashlib
import importlib.metadata
import json
import os
from pathlib import Path
import platform
import shutil
import subprocess
import sys
import time
import uuid

ROOT = Path(__file__).resolve().parent


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def get_key(data: dict, keys: list[str]):
    value = data
    for key in keys:
        if not isinstance(value, dict) or key not in value:
            raise ValueError(f"Missing success field: {'.'.join(keys)}")
        value = value[key]
    return value


def validate_report(data: dict, keys: list[str]) -> tuple[bool, int | None, int | None]:
    """Only JSON boolean True counts; additionally check all individual records."""
    if not isinstance(data, dict):
        raise ValueError('Report must be a JSON object')
    ok = get_key(data, keys) is True
    def inspect(value):
        if isinstance(value, dict):
            for key, item in value.items():
                if key in {'skipped', 'blocked'} and item is not False:
                    raise ValueError('Required report contains skipped/blocked work')
                if key == 'status' and item in ('skipped', 'blocked'):
                    raise ValueError('Required report contains skipped/blocked work')
                if key in {'pass', 'passed', 'passed_symbolic', 'passed_numeric'} and item is False:
                    raise ValueError('Individual assertion failed')
                if (key.endswith(('_validated', '_certified', '_certified_by_this_finite_test_run'))
                    or key == 'physical_n3llprime_payload') and item is not False:
                    raise ValueError('Unsupported accuracy claim: '+key)
                inspect(item)
        elif isinstance(value, list):
            for item in value: inspect(item)
    # Preserve the existing (False, passed, total) contract for ordinary failed
    # records while also checking nested aggregate assertions and scope metadata.
    try:
        inspect(data)
    except ValueError:
        ok = False
    tests = data.get('tests')
    total = passed = None
    if isinstance(tests, list):
        if not tests:
            return False, 0, 0
        if not all(isinstance(t, dict) and isinstance(t.get('passed'), bool) for t in tests):
            return False, None, None
        total = len(tests)
        passed = sum(t['passed'] is True for t in tests)
        ok = ok and total == passed
    # Some original component reports expose aggregate metrics rather than a
    # list of checks. Do not invent a count for those report formats.
    return ok, passed, total


def validate_execution(case: Path, report: Path, suite_id: str, run_id: str):
    if report.is_symlink() or (case/'execution.json').is_symlink():
        raise ValueError('Report/envelope must not be a symlink')
    envelope = json.loads((case/'execution.json').read_text())
    expected = dict(schema='k-space-suite-execution/1', suite_id=suite_id,
                    run_id=run_id, report_sha256=sha(report))
    if envelope != expected:
        raise ValueError('Wrong suite/run identity or report digest')
    return json.loads(report.read_text())


def validate_manifest(manifest: dict, root: Path):
    suites = manifest['suites']
    ids = [s['id'] for s in suites]
    if not ids or len(ids) != len(set(ids)):
        raise ValueError('Empty registry or duplicate suite IDs')
    required = manifest['required_suite_ids']
    if set(required) != set(ids) or len(required) != len(set(required)):
        raise ValueError('Required suites must all be registered and executed')
    scripts = [s['script'] for s in suites]
    discovered = {str(p.relative_to(root)) for p in (root/'validators').glob('*.py')}
    if len(scripts) != len(set(scripts)) or set(scripts) != discovered:
        raise ValueError('Validator discovery and suite registry disagree')
    for spec in suites:
        for name in (spec['id'], spec['output']):
            if not name or Path(name).name != name or name in {'.', '..'}:
                raise ValueError('Suite IDs/output names must be single safe path components')
    return suites


def git_state() -> dict:
    def git(*args):
        try:
            p = subprocess.run(['git', '-C', str(ROOT), *args], capture_output=True,
                               text=True, timeout=10, check=False)
            return p.stdout.strip() if p.returncode == 0 else None
        except (OSError, subprocess.TimeoutExpired):
            return None
    top = git('rev-parse', '--show-toplevel')
    if not top or Path(top).resolve() != ROOT:
        return {'commit': None, 'dirty': None, 'note': 'standalone package, not repository HEAD'}
    status = git('status', '--porcelain')
    return {'commit': git('rev-parse', 'HEAD'), 'dirty': bool(status),
            'branch': git('branch', '--show-current')}


def checked_output_dir(value: str) -> Path:
    p = Path(value).expanduser().absolute()
    allowed = (ROOT/'outputs').resolve()
    # Never follow a symlinked output path when cleaning.
    for q in [p, *p.parents]:
        if q == ROOT.parent:
            break
        if q.is_symlink():
            raise ValueError(f'Symlink in output path: {q}')
    resolved = p.resolve()
    if not resolved.is_relative_to(allowed) or resolved == allowed:
        raise ValueError('Output directory must be a child of this repository\'s outputs/ directory')
    return resolved


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--suite', action='append', help='Select a suite ID (repeatable); default: all')
    ap.add_argument('--list', action='store_true', help='List IDs without running')
    ap.add_argument('--clean', action='store_true', help='Remove only the selected output directory')
    ap.add_argument('--output-dir', default=str(ROOT/'outputs/latest'))
    ap.add_argument('--timeout', type=float, default=600, help='Per-suite timeout in seconds')
    args = ap.parse_args(argv)
    manifest = json.loads((ROOT/'manifests/verification_suite.json').read_text())
    try:
        suites = validate_manifest(manifest, ROOT)
    except (KeyError, ValueError) as e:
        ap.error(str(e))
    ids = [s['id'] for s in suites]
    if len(ids) != len(set(ids)):
        ap.error('Duplicate suite IDs')
    if args.list:
        print('\n'.join(f"{s['id']:20} {s['description']}" for s in suites))
        return 0
    if args.timeout <= 0:
        ap.error('timeout must be positive')
    if args.suite:
        missing = set(args.suite)-set(ids)
        if missing:
            ap.error(f'Unknown suite IDs: {sorted(missing)}')
        suites = [s for s in suites if s['id'] in args.suite]
    try:
        out = checked_output_dir(args.output_dir)
    except ValueError as e:
        ap.error(str(e))
    if args.clean and out.exists():
        shutil.rmtree(out)
    out.mkdir(parents=True, exist_ok=True)
    env = os.environ.copy()
    # Explicitly fix test threading; reported numerical tolerances remain the
    # criterion, not bitwise equality across BLAS implementations.
    for k in ['OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS']:
        env[k] = '1'
    env['PYTHONDONTWRITEBYTECODE'] = '1'
    started = dt.datetime.now(dt.timezone.utc).isoformat()
    run_id = uuid.uuid4().hex
    results = []
    for spec in suites:
        case = out/spec['id']
        if case.is_symlink():
            raise ValueError(f'Refusing symlink output {case}')
        if case.exists():
            shutil.rmtree(case)
        case.mkdir()
        script = (ROOT/spec['script']).resolve()
        if not script.is_relative_to((ROOT/'validators').resolve()) or not script.is_file():
            raise ValueError(f'Invalid registered script: {spec["script"]}')
        report = case/spec['output']
        command = [sys.executable, str(ROOT/'tools/execute_suite.py'),
                   spec['id'], run_id, str(report), str(script)]
        if spec['invocation'] == 'out':
            command += ['--out', str(case)]
        elif spec['invocation'] != 'environment':
            raise ValueError('Unsupported invocation')
        env['KCSS_OUTPUT_DIR'] = str(case)
        print(f"[run] {spec['id']}: {spec['description']}", flush=True)
        t0 = time.monotonic()
        rc, error, log, passed, total = -1, None, '', None, None
        ok = False
        try:
            process = subprocess.run(command, cwd=ROOT, env=env, capture_output=True,
                                     text=True, timeout=args.timeout, check=False)
            rc = process.returncode
            log = process.stdout + process.stderr
            if rc == 0:
                data = validate_execution(case, report, spec['id'], run_id)
                ok, passed, total = validate_report(data, spec['pass_key'])
                if not ok:
                    error = 'Failed assertion, skipped work, or unsupported claim metadata'
            else:
                error = f'process exit {rc}'
        except (OSError, ValueError, subprocess.TimeoutExpired) as e:
            error = f'{type(e).__name__}: {e}'
        (case/'console.log').write_text(log, encoding='utf-8')
        result = dict(id=spec['id'], category=spec['category'], passed=bool(ok),
                      check_passed=passed, check_total=total, returncode=rc,
                      error=error, elapsed_seconds=time.monotonic()-t0,
                      script=spec['script'], script_sha256=sha(script),
                      report=str((case/spec['output']).relative_to(out)))
        results.append(result)
        print(f"[{'PASS' if ok else 'FAIL'}] {spec['id']}" +
              (f' {passed}/{total}' if total is not None else '') +
              (f' — {error}' if error else ''), flush=True)
    code = [ROOT/'run_all_checks.py']
    for name in ['validators','src','tools','manifests','data']:
        code += [p for p in (ROOT/name).rglob('*') if p.is_file() and p.suffix in ['.py','.json']]
    coefficient_path = ROOT/'data/coefficients_demo.json'
    coefficient_input = json.loads(coefficient_path.read_text()) if coefficient_path.exists() else None
    if coefficient_input is not None:
        coefficient_input = dict(path='data/coefficients_demo.json',
                                 sha256=sha(coefficient_path), input=coefficient_input,
                                 coefficient_terms={k:len(coefficient_input[k])
                                     for k in ['beta','cusp','noncusp','d']})
    summary = dict(schema='k-space-verification-run/1', run_id=run_id, started_utc=started,
                   finished_utc=dt.datetime.now(dt.timezone.utc).isoformat(),
                   all_suites_passed=all(r['passed'] for r in results),
                   suites_passed=sum(r['passed'] for r in results), suites_total=len(results),
                   requested_scope='all' if not args.suite else args.suite,
                   full_required_run=not bool(args.suite),
                   required_suite_ids=manifest['required_suite_ids'],
                   unselected_suite_ids=[i for i in ids if i not in [r['id'] for r in results]],
                   optional_suites_skipped=0,
                   results=results, git=git_state(),
                   coefficient_input=coefficient_input,
                   environment={'python':sys.version,'platform':platform.platform(),
                       'packages':{n:importlib.metadata.version(n) for n in ['numpy','scipy','sympy','mpmath']},
                       'threads':{k:env[k] for k in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS']}},
                   input_sha256={str(p.relative_to(ROOT)):sha(p) for p in sorted(set(code))},
                   scientific_scope={'numerical_inputs':'declared lower-order QCD demonstration and labeled formal fixtures',
                       'physical_three_loop_payload_validated':False,
                       'physical_N3LLprime_prediction_validated':False,
                       'continuum_grid_error_certified':False,
                       'data_extraction_performed':False,
                       'interpretation':'Checks supplement the documented analytic arguments; counts are not a theorem certificate.'})
    (out/'run_all_checks_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    lines=['k-space verification run',f"all_suites_passed: {summary['all_suites_passed']}",
           f"suites: {summary['suites_passed']}/{summary['suites_total']}",'',
           'Suite counts refer to different test scopes. No grand physics pass count is inferred.']
    for r in results:
        count=f" ({r['check_passed']}/{r['check_total']} records)" if r['check_total'] is not None else ' (aggregate report)'
        lines.append(f"{'PASS' if r['passed'] else 'FAIL'} {r['id']}{count}: {r['report']}")
    lines += ['', 'physical_N3LLprime_prediction_validated: False', 'continuum_grid_error_certified: False',
              'Source hashes, environment, and categories: run_all_checks_summary.json']
    (out/'run_all_checks_report.txt').write_text('\n'.join(lines)+'\n')
    print('\n'.join(lines))
    return 0 if summary['all_suites_passed'] else 1


if __name__ == '__main__':
    raise SystemExit(main())
