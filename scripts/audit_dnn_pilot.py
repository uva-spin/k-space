"""Post-fit native-pilot diagnostics; no refitting or interpretation as coverage."""
from pathlib import Path
import sys,json,hashlib
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from native_mixture import MixtureKernel,VARIANCES
from mixture_dnn import MixtureDNN,Objective,truth_weights
repo=Path(__file__).resolve().parents[1];ev=repo.parent/'evidence'
d=json.loads((ev/'data/data.json').read_text());c=json.loads((ev/'closure/closure.json').read_text());assert d['passed'] and c['passed']
net=MixtureDNN();kernel=MixtureKernel.load(ev/'kernel_dense/native_kernel.npz');data=np.array(d['data']);sigma=np.array(d['sigma']);checks=[]
for fit in d['fits']:
    p=np.load(ev/'data'/f"parameters_seed_{fit['seed']}.npz")['parameters'];theta=p[:net.size]
    values=[]
    for n in (257,513):
        x=np.geomspace(.04,.70,n);h=np.log(x[1]/x[0]);w=net.forward(theta,x);cur=(w[2:]-2*w[1:-1]+w[:-2])/h**2
        values.append(.5*.1*float(np.mean(np.sum(cur*cur,axis=1))))
    checks.append(dict(seed=fit['seed'],smoothness_penalty_257=values[0],smoothness_penalty_513=values[1],absolute_loss_change=abs(values[1]-values[0])))
best=d['fits'][d['best_index']];p=np.load(ev/'data'/f"parameters_seed_{best['seed']}.npz")['parameters'];theta=p[:net.size]
# Exact reconstruction-only change of E_tab to E(kinematic means): row factor.
factors=np.array([np.sqrt(r['Qmean']**2+r['qtmean']**2+(38.8*r['xfmean']/2)**2)/r['Emean'] for r in d['rows']])
pred=np.array(best['bare_predictions'])*factors;N=float((np.sum(pred*data/sigma**2)+1/.065**2)/(np.sum(pred*pred/sigma**2)+1/.065**2))
chi=float(np.sum(((N*pred-data)/sigma)**2)+((N-1)/.065)**2)
xs=np.geomspace(.05,.65,101);profiles=[];boundary_moments=[]
for fit in d['fits']:
    par=np.load(ev/'data'/f"parameters_seed_{fit['seed']}.npz")['parameters'];w=net.forward(par[:net.size],xs);profiles.append(w);boundary_moments.append((w@VARIANCES).tolist())
cp=np.load(ev/'closure'/f"parameters_seed_{c['best_seed']}.npz")['parameters'];probe=np.array(c['probe_x']);k=np.linspace(0,20,4001)
gen=truth_weights(probe)@(np.exp(-k[None,:]**2/VARIANCES[:,None])/(np.pi*VARIANCES[:,None]));got=net.density(cp[:net.size],probe,k)
l1=np.trapezoid(abs(got-gen)*2*np.pi*k,k,axis=1)
spread=np.ptp(np.array([f['predictions'] for f in d['fits']]),axis=0)/np.maximum(np.abs(best['predictions']),1e-14)
out=dict(scope='Initialization and declared-prior diagnostics, not calibrated uncertainty or a refit ensemble',smoothness_grid_checks=checks,seed_prediction_max_fractional_range=float(spread.max()),boundary_x=xs.tolist(),boundary_weights=np.array(profiles).tolist(),basis_variances_GeV2=VARIANCES.tolist(),boundary_second_moments_GeV2=boundary_moments,boundary_moment_scope='Moment of normalized Gaussian-mixture boundary only; not a finite moment of the full renormalized TMD',closure_profile_probe_x=probe.tolist(),closure_profile_L1_distances=l1.tolist(),closure_scope='Observable closure does not guarantee full-profile recovery; these differences expose weak/model-dependent directions',measurement_variation=dict(kind='Reconstruction-only kinematic-mean energy alternative; neural input fixed, only normalization profiled',normalization=N,chi2_total=chi,full_neural_refit=False),script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
(ev/'pilot_diagnostics.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({k:v for k,v in out.items() if k not in ('boundary_x','boundary_second_moments_GeV2')},indent=2))
