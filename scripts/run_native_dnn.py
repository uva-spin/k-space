"""Native-only neural training, independent closure and conditional data pilot."""
from pathlib import Path
import sys,json,hashlib,time,argparse
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from native_mixture import MixtureKernel,VARIANCES
from mixture_dnn import MixtureDNN,Objective,truth_weights,train,function_sensitivity
p=argparse.ArgumentParser();p.add_argument('--stage',choices=['closure','data'],required=True);p.add_argument('--out',type=Path,required=True);p.add_argument('--seeds',default='11,29,47');p.add_argument('--maxiter',type=int,default=2500);p.add_argument('--evidence-root',type=Path);args=p.parse_args();args.out.mkdir(parents=True,exist_ok=True)
repo=Path(__file__).resolve().parents[1];ev=args.evidence_root or repo.parent/'evidence';start=time.monotonic()
validation=json.loads((ev/'kernel_dense/kernel_validation.json').read_text());assert validation['passed']
for name,h in validation['source_sha256'].items():assert hashlib.sha256((repo/name).read_bytes()).hexdigest()==h,name
kernel=MixtureKernel.load(ev/'kernel_dense/native_kernel.npz');reference=MixtureKernel.load(ev/'kernel_dense/refined_reference_kernel.npz')
net=MixtureDNN();rows=kernel.rows;actual=np.array([r['value'] for r in rows]);sigma_actual=np.sqrt([r['stat']**2+r['syst']**2 for r in rows])
truth=reference.evaluate(truth_weights(reference.xa),truth_weights(reference.xb))
if args.stage=='closure':
    data=truth;sigma=.03*truth;mask=np.array([r['Qlow']!=7.2 for r in rows]);normalization=False;regularization=1e-5
else:
    accepted=json.loads((ev/'closure/closure.json').read_text());assert accepted['passed'],'Independent closure must pass first'
    data=actual;sigma=sigma_actual;mask=np.ones(len(rows),bool);normalization=True;regularization=1e-4
objective=Objective(net,kernel,data,sigma,mask,normalization=normalization,regularization=regularization)
refobj=Objective(net,reference,data,sigma,mask,normalization=normalization,regularization=regularization)
# Check end-to-end analytic gradients before training with directional differences.
test=net.initialize(100);test=np.r_[test,0.] if normalization else test
loss,grad=objective(test);rng=np.random.default_rng(103);checks=[]
for _ in range(5):
    d=rng.normal(size=len(test));d/=np.linalg.norm(d);h=2e-5
    fd=(objective(test+h*d)[0]-objective(test-h*d)[0])/(2*h)
    checks.append(abs(fd-grad@d)/max(1,abs(fd),abs(grad@d)))
assert max(checks)<2e-6
results=[];parameters=[]
for seed in map(int,args.seeds.split(',')):
    print('TRAIN',args.stage,seed,flush=True)
    par,opt=train(objective,seed,maxiter=args.maxiter,progress=lambda s,i,l:print('OPT',seed,s,i,l,flush=True))
    pred=objective.predict(par);fine=refobj.predict(par);N=float(np.exp(par[-1])) if normalization else 1.
    refine=float(np.max(abs(pred/fine-1)))
    directional=[]
    for _ in range(3):
        d=rng.normal(size=len(par));d/=np.linalg.norm(d);h=1e-5
        dn=(objective.predict(par+h*d)-objective.predict(par-h*d))/(2*h)
        dr=(refobj.predict(par+h*d)-refobj.predict(par-h*d))/(2*h)
        directional.append(float(np.max(abs(dn-dr))/max(np.max(abs(dr)),1e-12)))
    points=((N*pred-data)/sigma)**2;norm=((N-1)/.065)**2 if normalization else 0.
    report=dict(seed=seed,optimizer=opt,normalization=N,chi2_points=float(points[mask].sum()),chi2_normalization=norm,chi2_total=float(points[mask].sum()+norm),training_rows=int(mask.sum()),heldout_rows=int((~mask).sum()),training_max_relative=float(np.max(abs(N*pred[mask]/data[mask]-1))),heldout_max_relative=float(np.max(abs(N*pred[~mask]/data[~mask]-1))) if (~mask).any() else None,predictions=(N*pred).tolist(),bare_predictions=pred.tolist(),reference_predictions=(N*fine).tolist(),quadrature_relative=refine,directional_gradient_refinement=max(directional),function_sensitivity=function_sensitivity(net,par,kernel,sigma),parameter_count=net.size+int(normalization),effective_dof_not_assigned=True)
    parameters.append(par);results.append(report)
    np.savez_compressed(args.out/f'parameters_seed_{seed}.npz',parameters=par)
    print('RESULT',{k:v for k,v in report.items() if k not in ('predictions','bare_predictions','reference_predictions','function_sensitivity')},flush=True)
# Choose the fit by training objective only; never use held-out truth to select.
best=min(range(len(results)),key=lambda i:results[i]['optimizer']['loss']);r=results[best]
passed=all(x['quadrature_relative']<3e-4 and x['directional_gradient_refinement']<5e-4 for x in results)
if args.stage=='closure':passed=passed and r['optimizer']['success'] and r['training_max_relative']<.005 and r['heldout_max_relative']<.02
else:passed=passed and r['optimizer']['success']
xs=np.array([.06,.1,.2,.35,.5]);ks=np.linspace(0,5,151)
profile=net.density(parameters[best][:net.size],xs,ks)
truth_profile=truth_weights(xs)@(np.exp(-ks[None,:]**2/VARIANCES[:,None])/(np.pi*VARIANCES[:,None]))
result=dict(stage=args.stage,passed=bool(passed),scope='Positive Gaussian-basis mixture DNN pilot; shared flavor/legs, fixed CS and perturbative prescription. No full functional identification or calibrated uncertainty claim.',architecture=list(net.widths),network_parameters=net.size,regularization=regularization,smoothness=dict(coefficient=.1,grid_points=257,x_domain=[.04,.70],definition='0.5*coefficient*mean(sum_i (d2 w_i / dlogx2)^2), centered finite differences'),regularization_definition='Objective includes 0.5*coefficient*sum(theta^2), including biases, excluding normalization nuisance',normalization_prior=.065 if normalization else None,closure_generation='Independent refined Bessel predictions from analytic x-correlated mixture weights; generator formula not used by optimizer',selection='Lowest training objective; held-out points do not select the seed',gradient_directional_errors=checks,rows=rows,data=data.tolist(),sigma=sigma.tolist(),training_mask=mask.tolist(),fits=results,best_index=best,best_seed=r['seed'],probe_x=xs.tolist(),probe_k=ks.tolist(),best_input_profiles=profile.tolist(),truth_input_profiles=truth_profile.tolist() if args.stage=='closure' else None,profile_recovery_scope='Full-profile differences are retained as identifiability/model diagnostics, not a claimed recovery guarantee',elapsed_seconds=time.monotonic()-start,script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),network_sha256=hashlib.sha256((repo/'src/mixture_dnn.py').read_bytes()).hexdigest(),kernel_validation_sha256=hashlib.sha256((ev/'kernel_dense/kernel_validation.json').read_bytes()).hexdigest())
(args.out/(args.stage+'.json')).write_text(json.dumps(result,indent=2)+'\n');print('DONE',args.stage,'passed',passed,'best',r['seed'],'elapsed',result['elapsed_seconds'],flush=True)
raise SystemExit(0 if passed else 1)
