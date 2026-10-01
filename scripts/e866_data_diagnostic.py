#!/usr/bin/env python3
"""Conditional Gaussian-model comparison to E866; not a calibrated extraction.

The qT averaging convention, diagonal point-systematic treatment and restricted
x/flavor model are explicit. Inconsistent published E-mean values are not used.
"""
from pathlib import Path
import argparse,json,sys,hashlib,time
import numpy as np
from scipy.optimize import least_squares
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from e866_forward import NativeE866,GaussianCSInput
from lhapdf_provider import LHAPDFProvider
from benchmark_e866 import Reference,errors

def main():
 p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True);args=p.parse_args();args.out.mkdir(exist_ok=True,parents=True)
 root=Path(__file__).resolve().parents[1];runs=root.parent.parent
 # Require the completed forward/closure validation before using data.
 validation=json.loads((root.parent/'evidence/final/e866_forward.json').read_text())
 if not validation['passed']:raise RuntimeError('Forward and closure validation must pass first')
 for name,h in validation['code_sha256'].items():
  if hashlib.sha256((root/name).read_bytes()).hexdigest()!=h:raise RuntimeError('Untested source change: '+name)
 pdf=LHAPDFProvider(runs/'m2_prescription_2026-09-30/runtime/libkt_lhapdf.dylib',runs/'m1_resolution_2026-09-30/runtime/lhapdf/share/LHAPDF')
 manifest=json.loads((root/'manifests/e866_pp_bins.json').read_text());rows=[r for r in manifest['rows'] if r['selected_core']]
 engine=NativeE866(pdf,rows);reference=Reference(engine,n=48,nqt=24);data=np.array([r['value'] for r in rows]);sigma=np.sqrt([r['stat']**2+r['syst']**2 for r in rows]);reports=[]
 for mode in ('flat_pt','annular'):
  last={};start=time.monotonic()
  def predict(par):
   key=tuple(par)
   if key not in last:
    val=engine.evaluate(GaussianCSInput(par[0],par[1]),mode);pred=val[:,:,0].sum(axis=1);grad=val[:,:,1:].sum(axis=1)
    residual=np.r_[(par[2]*pred-data)/sigma,(par[2]-1)/.065]
    jac=np.vstack([np.column_stack([par[2]*grad,pred])/sigma[:,None],[0,0,1/.065]])
    last.clear();last[key]=(residual,jac,pred,val)
   return last[key]
  fit=least_squares(lambda p:predict(p)[0],[.7,.05,1.],jac=lambda p:predict(p)[1],bounds=([.05,0,.5],[3.,.5,1.5]),ftol=2e-8,xtol=2e-8,gtol=2e-8,max_nfev=35)
  residual,jac,pred,comp=predict(fit.x);sv=np.linalg.svd(jac,compute_uv=False)
  report=dict(mode=mode,success=bool(fit.success),parameters=dict(lambda0=float(fit.x[0]),g=float(fit.x[1]),normalization=float(fit.x[2])),active_bounds=fit.active_mask.tolist(),evaluations=fit.nfev,chi2_points=float(residual[:-1]@residual[:-1]),chi2_normalization=float(residual[-1]**2),chi2_total=float(residual@residual),points=len(rows),fit_parameters=3,prior_residuals=1,nominal_residual_dof=len(rows)+1-3,jacobian_singular_values=sv.tolist(),positive_core_prediction=bool((pred>0).all()),rows=[dict(id=r['id'],data=r['value'],point_sigma=float(sigma[i]),theory=float(pred[i]),theory_times_nuisance=float(fit.x[2]*pred[i]),pull=float(residual[i])) for i,r in enumerate(rows)],elapsed_seconds=time.monotonic()-start)
  report['best_point_reference_check']=errors(comp,reference.evaluate(GaussianCSInput(fit.x[0],fit.x[1]),mode))
  report['best_point_numerically_validated']=report['best_point_reference_check']['total_relative']<3e-4
  reports.append(report);print(mode,{k:v for k,v in report.items() if k!='rows'},flush=True)
 result=dict(schema='E866-conditional-data-comparison/1',scope='Restricted Gaussian x/flavor-independent model, massless nf5 photon W-only, 41 declared bins. This is not a general proton-TMD or CS-kernel extraction.',statistical_assumptions='Quoted stat and point-to-point syst combined in quadrature and treated as diagonal; one common 6.5% Gaussian normalization nuisance. No full covariance, PDF uncertainty or coverage claim.',data_metadata_issue='Published mean-energy column is inconsistent with stated bins. It is not used or repaired. Dependence of reported invariant on finite-bin conversion remains unconfirmed.',bounds={'lambda0':[.05,3.],'g':[0,.5],'normalization':[.5,1.5]},fits=reports,successful=all(r['success'] and r['best_point_numerically_validated'] for r in reports),code_sha256=validation['code_sha256'],script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),data_manifest_sha256=hashlib.sha256((root/'manifests/e866_pp_bins.json').read_bytes()).hexdigest())
 (args.out/'data_diagnostic.json').write_text(json.dumps(result,indent=2)+'\n')
 return 0 if result['successful'] else 1
if __name__=='__main__':raise SystemExit(main())
