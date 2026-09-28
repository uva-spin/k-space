#!/usr/bin/env python3
"""Offline checks for local Markdown links, fences and display-math structure.

This does not emulate GitHub/MathJax rendering and is not a physics validator.
"""
from __future__ import annotations
import argparse
import re
from pathlib import Path
from urllib.parse import unquote, urlsplit

ROOT = Path(__file__).resolve().parents[1]
EXCLUDE = {'.git','.venv','venv','outputs','__pycache__','.pytest_cache','build','dist'}


def heading_id(heading: str) -> str:
    value=re.sub(r'<[^>]*>', '', heading).strip().lower()
    value=re.sub(r'[^\w\- ]', '', value)
    return value.replace(' ', '-')


def check_file(path: Path, root: Path) -> list[str]:
    text=path.read_text(encoding='utf-8'); issues=[]; plain=[]; math=[]; fence=None; math_lines=[]
    for line in text.splitlines():
        match=re.match(r'^\s*(`{3,}|~{3,})(.*)$',line)
        if match:
            marker,lang=match.groups()
            if fence is None:
                fence=marker;in_math=lang.strip()=='math';math_lines=[]
            elif marker[0]==fence[0] and len(marker)>=len(fence) and not lang.strip():
                if in_math:math.append('\n'.join(math_lines))
                fence=None
            else:
                if in_math:math_lines.append(line)
            continue
        if fence is None:plain.append(line)
        elif in_math:math_lines.append(line)
    if fence:issues.append(f'{path.relative_to(root)}: unclosed code fence')
    for block in math:
        if not block.strip():issues.append(f'{path.relative_to(root)}: empty math block')
        begins=re.findall(r'\\begin\{([^}]+)\}',block);ends=re.findall(r'\\end\{([^}]+)\}',block)
        if begins!=list(reversed(ends)) and sorted(begins)!=sorted(ends):
            issues.append(f'{path.relative_to(root)}: unmatched math environments')
    visible='\n'.join(plain)
    code_spans=[m.span() for m in re.finditer(r'(`+).*?\1',visible)]
    for match in re.finditer(r'(?<!!)\[[^\]\n]+\]\(([^\s)]+)(?:\s+"[^"]*")?\)',visible):
        if any(start <= match.start() < end for start,end in code_spans): continue
        target=match.group(1)
        parsed=urlsplit(target)
        if parsed.scheme or parsed.netloc:continue
        rel=unquote(parsed.path);resolved=(path.parent/rel).resolve() if rel else path.resolve()
        if not resolved.is_relative_to(root.resolve()) or not resolved.exists():
            issues.append(f'{path.relative_to(root)}: broken/escaping local link {target}')
        elif parsed.fragment and resolved.suffix=='.md':
            headings=[heading_id(h) for h in re.findall(r'^#{1,6}\s+(.+?)\s*#*$',resolved.read_text(encoding='utf-8'),re.M)]
            if unquote(parsed.fragment) not in headings:
                issues.append(f'{path.relative_to(root)}: unresolved local heading {target}')
    return issues


def scan(root: Path) -> list[str]:
    result=[]
    for path in sorted(root.rglob('*.md')):
        if set(path.relative_to(root).parts)&EXCLUDE:continue
        result.extend(check_file(path,root))
    return result


def main() -> int:
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('root',type=Path,nargs='?',default=ROOT)
    args=parser.parse_args();root=args.root.resolve();issues=scan(root)
    if issues:print('\n'.join(issues));return 1
    print('PASS: local Markdown links, code fences, and display-math structure')
    print('Rendering and scientific correctness are separate checks.')
    return 0

if __name__=='__main__':raise SystemExit(main())
