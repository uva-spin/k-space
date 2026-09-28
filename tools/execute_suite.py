"""Adapt existing validator reports to a run-identified envelope, without physics edits."""
from __future__ import annotations
import hashlib
import json
from pathlib import Path
import runpy
import sys


def main():
    suite_id, run_id, report_name, script, *arguments = sys.argv[1:]
    sys.argv = [script, *arguments]
    try:
        runpy.run_path(script, run_name='__main__')
    except SystemExit as exc:
        if exc.code not in (None, 0):
            raise
    # Only a successful validator with a newly created, valid report gets an
    # envelope. The parent independently checks identity, digest and assertions.
    report = Path(report_name)
    if report.is_symlink():
        raise ValueError('Report must not be a symlink')
    raw = report.read_bytes()
    json.loads(raw)
    envelope = dict(schema='k-space-suite-execution/1', suite_id=suite_id,
                    run_id=run_id, report_sha256=hashlib.sha256(raw).hexdigest())
    (report.parent/'execution.json').write_text(json.dumps(envelope, indent=2)+'\n')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
