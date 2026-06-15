#!/usr/bin/env python3
"""Sanity-check that repository text files were not uploaded with newlines stripped.

This is intentionally lightweight. It catches the failure mode where Markdown,
Python, shell, Makefile, TOML, and CFF files are collapsed into one or a few very
long lines. It does not enforce style formatting.
"""
from __future__ import annotations

import argparse
from pathlib import Path

TEXT_EXTENSIONS = {'.md', '.py', '.sh', '.toml', '.yml', '.yaml', '.cff', '.txt', '.mk'}
TEXT_NAMES = {'Makefile', 'requirements.txt'}
SKIP_PARTS = {'.git', '.venv', '__pycache__', 'outputs'}


def is_checked(path: Path) -> bool:
    if any(part in SKIP_PARTS for part in path.parts):
        return False
    return path.suffix in TEXT_EXTENSIONS or path.name in TEXT_NAMES


def main() -> int:
    parser = argparse.ArgumentParser(description='Check for newline-stripped text files.')
    parser.add_argument('root', nargs='?', default='.', help='Repository root')
    parser.add_argument('--max-line-length', type=int, default=2000)
    args = parser.parse_args()

    root = Path(args.root).resolve()
    failures: list[str] = []

    for path in sorted(p for p in root.rglob('*') if p.is_file() and is_checked(p.relative_to(root))):
        rel = path.relative_to(root)
        text = path.read_text(encoding='utf-8')
        lines = text.splitlines()
        if len(lines) <= 1 and path.name not in {'.gitkeep'}:
            failures.append(f'{rel}: only {len(lines)} line(s)')
            continue
        too_long = max((len(line) for line in lines), default=0)
        if too_long > args.max_line_length:
            failures.append(f'{rel}: max line length {too_long} > {args.max_line_length}')

    if failures:
        print('line-break check: FAIL')
        for item in failures:
            print(' -', item)
        return 1

    print('line-break check: PASS')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
