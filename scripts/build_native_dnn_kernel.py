from pathlib import Path
import sys,json,time,hashlib,argparse
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from native_mixture import build_native
from lhapdf_provider import LHAPDFProvider
from native_mixture_reference import build_reference
p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True);p.add_argument('--nQ',type=int,default=4);p.add_argument('--nx',type=int,default=8);args=p.parse_args();args.out.mkdir(parents=True,exist_ok=True)
repo=Path(__file__).resolve().parents[1];runs=repo.parent.parent;start=time.monotonic()
rows=[r for r in json.loads((repo/'manifests/e866_pp_bins.json').read_text())['rows'] if r['selected_core']]
pdf=LHAPDFProvider(runs/'m2_prescription_2026-09-30/runtime/libkt_lhapdf.dylib',runs/'m1_resolution_2026-09-30/runtime/lhapdf/share/LHAPDF')
native=build_native(pdf,rows,nQ=args.nQ,nx=args.nx,nmu=96,order=240,progress=lambda Q,n:print('NATIVE',Q,n,time.monotonic()-start,flush=True));native.save(args.out/'native_kernel.npz')
ref=build_reference(pdf,rows,nQ=args.nQ,nx=args.nx,nmu=96,order=240);ref.save(args.out/'reference_kernel.npz')
fine=build_reference(pdf,rows,nQ=args.nQ*3//2,nx=args.nx*3//2,nmu=128,order=320,nb=48,nq=24);fine.save(args.out/'refined_reference_kernel.npz')
mask=np.abs(ref.kernel)>1e-16
err=float(np.max(abs(native.kernel-ref.kernel)[mask]/abs(ref.kernel[mask])))
cerr=float(np.max(np.where(mask[...,None],abs(native.components-ref.components)/np.maximum(abs(ref.kernel[...,None]),1e-16),0)))
# Test x-correlated non-Gaussian inputs and their coefficient derivatives later;
# these basis-level checks cover every admitted pair and all six insertions.
result=dict(native_pair_kernel_relative=err,separate_component_over_total=cerr,thresholds={'native_pair_kernel_relative':3e-4,'separate_component_over_total':3e-4},passed=err<3e-4 and cerr<3e-4,elapsed_seconds=time.monotonic()-start,settings=native.settings,script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),source_sha256={str(x.relative_to(repo)):hashlib.sha256(x.read_bytes()).hexdigest() for x in (repo/'src').glob('*.py') if x.name!='mixture_dnn.py'})
(args.out/'kernel_validation.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2),flush=True)
if not result['passed']:raise SystemExit(1)
