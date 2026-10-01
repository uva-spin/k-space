#!/usr/bin/env python3
"""Check the human-facing source tree and locally generated verification files."""
from __future__ import annotations
import argparse
import ast
import hashlib
import json
from pathlib import Path
import subprocess

ROOT=Path(__file__).resolve().parents[1]
IGNORE={'.git','.venv','venv','__pycache__','.pytest_cache','.mypy_cache','.ruff_cache'}
BANNED={'.tex','.ltx','.bib','.bbl','.aux','.toc','.fls','.fdb_latexmk','.synctex',
        '.pdf','.docx','.zip','.tar','.tgz','.gz','.7z','.rar','.dvi','.lof','.lot',
        '.nav','.snm','.vrb','.bcf','.blg','.run.xml'}
AGENT_FILES={'agents.md','claude.md','gemini.md','codex_task.md','launch_codex.sh',
             'skill.md','.cursorrules','.cursorignore','.copilot-instructions.md',
             'copilot-instructions.md','launch_codex.sh',
             'k_space_cleanup_prompt.md'}

def source_issue(rel: Path) -> str | None:
    """Rules for files published in Git, including otherwise ignored files."""
    parts=[part.lower() for part in rel.parts]
    if rel.name.lower() in AGENT_FILES or any(
        part.startswith('.codex') or part in {'.claude','.cursor','.agents','.agent','agents','skills'}
        for part in parts
    ):
        return f'agent instructions/automation artifact: {rel}'
    if parts[0] in {'outputs','reference_outputs','backups','branch_reports'}:
        return f'generated report or local work artifact must remain untracked: {rel}'
    if '__pycache__' in parts or rel.suffix.lower() in {'.pyc','.pyo','.log'}:
        return f'generated cache/log artifact: {rel}'
    return None

def prohibited_artifact(rel: Path) -> bool:
    # Scientific PDF plots may be published in figures/ or regenerated in
    # ignored outputs/. The tracked-inventory check still rejects every output
    # artifact. This does not permit tracked manuscripts or TeX build products.
    if rel.suffix.lower() == '.pdf' and rel.parts[0] in {'figures', 'outputs'}:
        return False
    return any(rel.name.lower().endswith(ext) for ext in BANNED)

def tracked_files(root: Path) -> list[Path]:
    """Do not accidentally inspect an enclosing repository's inventory."""
    try:
        top=subprocess.run(['git','-C',str(root),'rev-parse','--show-toplevel'],
                           capture_output=True,text=True,check=False,timeout=10)
        if top.returncode or Path(top.stdout.strip()).resolve()!=root.resolve():
            return []
        listing=subprocess.run(['git','-C',str(root),'ls-files','-z'],
                               capture_output=True,text=True,check=True,timeout=10)
        return [Path(p) for p in listing.stdout.split('\0') if p]
    except (OSError,subprocess.SubprocessError):
        return []

def scan(root: Path) -> list[str]:
    issues=[]
    for p in root.rglob('*'):
        rel=p.relative_to(root)
        if set(rel.parts)&IGNORE:
            continue
        if p.is_symlink():
            issues.append(f'symlink requires explicit review: {rel}')
            continue
        if not p.is_file():
            continue
        if rel.parts[0] == 'outputs':
            issue=source_issue(rel)
            if issue and 'agent instructions' in issue: issues.append(issue)
        if prohibited_artifact(rel):
            issues.append(f'manuscript/build/archive artifact: {rel}')
        # Locally generated reports are allowed under outputs/. Other source
        # files follow the same rules as the tracked inventory.
        if rel.parts[0] != 'outputs':
            issue=source_issue(rel)
            if issue: issues.append(issue)
        if p.suffix=='.py':
            try: ast.parse(p.read_text(encoding='utf-8'), filename=str(rel))
            except (SyntaxError,UnicodeError) as e: issues.append(f'invalid Python: {rel}: {e}')
        if p.suffix in {'.py','.md','.toml','.yml','.cff','.txt'}:
            raw=p.read_bytes()
            if b'\r\n' in raw: issues.append(f'CRLF rather than declared LF: {rel}')
            if raw and not raw.endswith(b'\n'): issues.append(f'missing final newline: {rel}')
    for rel in tracked_files(root):
        issue=source_issue(rel)
        if issue and issue not in issues: issues.append(issue)
        if prohibited_artifact(rel):
            issue=f'manuscript/build/archive artifact: {rel}'
            if issue not in issues: issues.append(issue)
    manifest=root/'manifests/verification_suite.json'
    if manifest.exists():
        try:
            specs=json.loads(manifest.read_text())['suites']
            ids=[s['id'] for s in specs]
            if len(ids)!=len(set(ids)): issues.append('duplicate suite IDs')
            for s in specs:
                p=(root/s['script']).resolve()
                if not p.is_relative_to(root.resolve()) or not p.is_file():issues.append('missing/unsafe registered script '+s['script'])
        except (KeyError,ValueError) as e:issues.append(f'invalid suite manifest: {e}')
    return issues

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('root',nargs='?',type=Path,default=ROOT)
    p.add_argument('--check-prepared-hashes',action='store_true',help='Verify the imported scientific files against the source ledger')
    a=p.parse_args();root=a.root.resolve();issues=scan(root)
    if a.check_prepared_hashes:
        ledger=json.loads((root/'manifests/source_provenance.json').read_text())
        for f in ledger['files']:
            path=root/f['target']
            if not path.is_file() or hashlib.sha256(path.read_bytes()).hexdigest()!=f['prepared_sha256']:
                issues.append('prepared-source hash mismatch: '+f['target'])
    if issues:
        print('\n'.join(issues));return 1
    print('PASS: human-facing source tree, untracked generated reports, valid Python, and resolvable suites')
    if a.check_prepared_hashes:print('PASS: imported scientific source hashes match the source ledger')
    return 0
if __name__=='__main__':raise SystemExit(main())
