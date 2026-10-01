#!/usr/bin/env python3
"""Build the pinned external LHAPDF runtime on macOS or Linux, on explicit request."""
import argparse,hashlib,json,os,platform,shutil,subprocess,tarfile
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
p=argparse.ArgumentParser();p.add_argument('--work',type=Path,required=True);a=p.parse_args()
work=a.work.resolve()
if work.is_relative_to(ROOT) or ROOT.is_relative_to(work):raise SystemExit('Use an external, non-ancestor work directory')
system=platform.system()
if system not in ('Darwin','Linux'):raise SystemExit('Supported build platforms: macOS and Linux')
dl=work/'downloads';runtime=work/'runtime';dl.mkdir(parents=True,exist_ok=True);runtime.mkdir(exist_ok=True)
m=json.loads((ROOT/'manifests/css2_pdf.json').read_text())
for row in m['files']:
 dest=dl/row['name']
 if not dest.exists():subprocess.run(['curl','-fL','--retry','2',row['url'],'-o',str(dest)],check=True)
 if hashlib.sha256(dest.read_bytes()).hexdigest()!=row['sha256']:raise SystemExit('Source hash mismatch: '+str(dest))
source=runtime/'LHAPDF-6.5.6';prefix=runtime/'lhapdf'
if not source.exists():
 with tarfile.open(dl/'LHAPDF-6.5.6.tar.gz') as t:t.extractall(runtime,filter='data')
pdfdir=prefix/'share/LHAPDF/NNPDF40_nnlo_as_01180';pdfdir.mkdir(parents=True,exist_ok=True)
for name in ('NNPDF40_nnlo_as_01180.info','NNPDF40_nnlo_as_01180_0000.dat'):shutil.copy2(dl/name,pdfdir/name)
subprocess.run(['./configure','--prefix='+str(prefix),'--disable-python'],cwd=source,check=True)
subprocess.run(['make','-j4'],cwd=source,check=True);subprocess.run(['make','install'],cwd=source,check=True)
extension='.dylib' if system=='Darwin' else '.so'
flags=['-dynamiclib'] if system=='Darwin' else ['-shared','-fPIC']
bridge=runtime/('libkt_lhapdf'+extension)
subprocess.run([os.environ.get('CXX','c++'),'-std=c++11','-O2',*flags,str(ROOT/'src/lhapdf_bridge.cc'),
                '-I'+str(prefix/'include'),'-L'+str(prefix/'lib'),'-lLHAPDF',
                '-Wl,-rpath,'+str(prefix/'lib'),'-o',str(bridge)],check=True)
(work/'resources.json').write_text(json.dumps({'library':str(bridge),'pdf_data':str(prefix/'share/LHAPDF'),
    'platform':system,'LHAPDF_version':m['LHAPDF_version']},indent=2)+'\n')
print('Built resources:',work/'resources.json')
