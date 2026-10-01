#!/usr/bin/env python3
"""Retrieve and safely extract hash-pinned primary coefficient sources."""
import argparse
import gzip
import hashlib
import json
from pathlib import Path
import subprocess
import tarfile

ROOT=Path(__file__).resolve().parents[1]
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--out',type=Path,required=True)
args=parser.parse_args()
out=args.out.resolve()
if out.is_relative_to(ROOT) or ROOT.is_relative_to(out):
    parser.error('Use an external, non-ancestor source directory')
out.mkdir(parents=True,exist_ok=True)
names={'collins_rogers_1705.07167v3.tar':'collins_rogers_1705.07167v3',
       'hkm_cusp.tar':'hkm_cusp','mps_cusp.tar':'mps_cusp',
       'li_zhu.tar':'li_zhu','beta4.tar':'beta4'}
manifest=json.loads((ROOT/'manifests/css2_sources.json').read_text())
for row in manifest['records']:
    if 'archive' not in row:
        continue
    target=out/row['archive']
    if not target.exists():
        subprocess.run(['curl','--fail','--location','--retry','2',row['url'],'--output',str(target)],check=True)
    if hashlib.sha256(target.read_bytes()).hexdigest()!=row['sha256']:
        raise SystemExit('Source archive differs from the audited bytes: '+str(target))
    folder=out/names[row['archive']]
    folder.mkdir(exist_ok=True)
    if row['archive'] == 'beta4.tar':
        # This hash-pinned arXiv object is a single gzipped TeX file, despite
        # the historical .tar name in the source manifest.
        (folder/'source.tex').write_bytes(gzip.decompress(target.read_bytes()))
    else:
        with tarfile.open(target) as archive:
            archive.extractall(folder,filter='data')
for row in manifest['records']:
    if 'extracted_file' in row:
        target = out/row['extracted_file']
        if hashlib.sha256(target.read_bytes()).hexdigest() != row['sha256']:
            raise SystemExit('Extracted source differs from audited bytes: '+str(target))
print('Pinned primary sources available at',out)
