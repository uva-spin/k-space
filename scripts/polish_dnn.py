"""Continue iteration-limited, numerically validated neural fits from checkpoints."""
from pathlib import Path
import sys,json,hashlib,argparse,time
import numpy as np
from scipy.optimize import minimize
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from native_mixture import MixtureKernel
from mixture_dnn import MixtureDNN,Objective,function_sensitivity
p=argparse.ArgumentParser();p.add_argument('--input',type=Path,required=True);p.add_argument('--out',type=Path,required=True);p.add_argument('--evidence-root',type=Path);a=p.parse_args();a.out.mkdir(parents=True,exist_ok=True)
repo=Path(__file__).resolve().parents[1];ev=a.evidence_root or repo.parent/'evidence';old=json.loads((a.input/'data.json').read_text());assert old['stage']=='data'
assert old['network_sha256']==hashlib.sha256((repo/'src/mixture_dnn.py').read_bytes()).hexdigest()
assert all(f['quadrature_relative']<3e-4 and f['directional_gradient_refinement']<5e-4 for f in old['fits'])
net=MixtureDNN();kernel=MixtureKernel.load(ev/'kernel_dense/native_kernel.npz');reference=MixtureKernel.load(ev/'kernel_dense/refined_reference_kernel.npz')
obj=Objective(net,kernel,old['data'],old['sigma'],regularization=old['regularization']);refobj=Objective(net,reference,old['data'],old['sigma'],regularization=old['regularization'])
fits=[];pars=[];start=time.monotonic();rng=np.random.default_rng(830)
for previous in old['fits']:
    seed=previous['seed'];initial=np.load(a.input/f'parameters_seed_{seed}.npz')['parameters']
    opt=minimize(obj,initial,jac=True,method='L-BFGS-B',bounds=[(None,None)]*net.size+[(np.log(.5),np.log(1.5))],options=dict(maxiter=30000,maxfun=120000,maxls=40,maxcor=30,ftol=1e-12,gtol=2e-6))
    loss,grad=obj(opt.x);pred=obj.predict(opt.x);fine=refobj.predict(opt.x);N=float(np.exp(opt.x[-1]));res=(N*pred-np.array(old['data']))/old['sigma'];norm=((N-1)/.065)**2
    err=float(np.max(abs(pred/fine-1)));der=[]
    for _ in range(3):
        direction=rng.normal(size=len(opt.x));direction/=np.linalg.norm(direction);h=1e-5
        dn=(obj.predict(opt.x+h*direction)-obj.predict(opt.x-h*direction))/(2*h);dr=(refobj.predict(opt.x+h*direction)-refobj.predict(opt.x-h*direction))/(2*h)
        der.append(float(np.max(abs(dn-dr))/max(np.max(abs(dr)),1e-12)))
    f=dict(previous);f.update(optimizer=dict(success=bool(opt.success),message=str(opt.message),iterations=opt.nit,evaluations=opt.nfev,loss=loss,gradient_inf=float(np.max(abs(grad))),stage='continued from exact saved parameters with unchanged objective and tolerances'),previous_optimizer=previous['optimizer'],normalization=N,chi2_points=float(res@res),chi2_normalization=norm,chi2_total=float(res@res+norm),predictions=(N*pred).tolist(),bare_predictions=pred.tolist(),reference_predictions=(N*fine).tolist(),quadrature_relative=err,directional_gradient_refinement=max(der),training_max_relative=float(np.max(abs(N*pred/old['data']-1))),function_sensitivity=function_sensitivity(net,opt.x,kernel,np.array(old['sigma'])),polish_max_prediction_change=float(np.max(abs(N*pred/previous['predictions']-1))))
    fits.append(f);pars.append(opt.x);np.savez_compressed(a.out/f'parameters_seed_{seed}.npz',parameters=opt.x)
    print('POLISHED',seed,opt.success,f['chi2_total'],err,max(der),f['polish_max_prediction_change'],flush=True)
best=min(range(len(fits)),key=lambda i:fits[i]['optimizer']['loss']);result=dict(old)
result.update(fits=fits,best_index=best,best_seed=fits[best]['seed'],best_input_profiles=net.density(pars[best][:net.size],np.array(old['probe_x']),np.array(old['probe_k'])).tolist(),passed=all(f['optimizer']['success'] and f['quadrature_relative']<3e-4 and f['directional_gradient_refinement']<5e-4 for f in fits),polish_script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),polish_input_sha256=hashlib.sha256((a.input/'data.json').read_bytes()).hexdigest(),polish_elapsed_seconds=time.monotonic()-start)
(a.out/'data.json').write_text(json.dumps(result,indent=2)+'\n');print('PASSED',result['passed'],flush=True)
raise SystemExit(0 if result['passed'] else 1)
