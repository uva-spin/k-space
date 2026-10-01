"""Single-TMD snapshots from the accepted pilot; hard factor is excluded.

Native values and separate strict matching components are saved alongside an
independent analytic-Gaussian/Bessel check. No output positivity clipping occurs.
"""
from pathlib import Path
import sys,json,hashlib
import numpy as np
from scipy.interpolate import CubicSpline
from scipy.special import j0
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from e866_forward import HighScaleProfile,canonical_blocks,FLAVORS,nodes_weights
from css2_noncanonical import FixedFlavorWindow,make_fixed_rg
from css2_reference import profile_leg
from native_mixture import BatchRational,VARIANCES,FIXED_G
from mixture_dnn import MixtureDNN
from rational_radial import RadialGrid
from lhapdf_provider import LHAPDFProvider
repo=Path(__file__).resolve().parents[1];runs=repo.parent.parent;ev=repo.parent/'evidence'
result=json.loads((ev/'data/data.json').read_text());assert result['passed']
seed=result['best_seed'];net=MixtureDNN();theta=np.load(ev/'data'/f'parameters_seed_{seed}.npz')['parameters'][:net.size]
assert result['network_sha256']==hashlib.sha256((repo/'src/mixture_dnn.py').read_bytes()).hexdigest()
pdf=LHAPDFProvider(runs/'m2_prescription_2026-09-30/runtime/libkt_lhapdf.dylib',runs/'m1_resolution_2026-09-30/runtime/lhapdf/share/LHAPDF')
xs=np.array([.06,.1,.2,.35,.5]);weights=net.forward(theta,xs);grid=RadialGrid(640.,102400);selected=np.arange(0,np.searchsorted(grid.k,5.),8);k=grid.k[selected]
b=[];bw=[]
for lo,hi in zip([0,.08,.16,.25,.5,1,2,4,8],[.08,.16,.25,.5,1,2,4,8,16]):
    z,w=nodes_weights(lo,hi,48);b.extend(z);bw.extend(w)
b=np.array(b);bw=np.array(bw)*b/(2*np.pi);reports=[]
for Q in (7.,12.):
    profile=HighScaleProfile(Q);low=float(profile.mu(np.inf));rg=make_fixed_rg(FixedFlavorWindow(5,4.92,172.5),low,91.1876,.1180024,Q)
    t=np.linspace(-16*np.finfo(float).eps,np.log(Q/rg.mu0)+16*np.finfo(float).eps,96)
    cache=CubicSpline(t,np.array([canonical_blocks(pdf,rg,xs,rg.mu0*np.exp(tt),240) for tt in t]),axis=0,extrapolate=False)
    def symbol(s):return cache(np.log(profile.mu(s)/rg.mu0))*np.asarray(profile_leg(rg,profile,s))[...,None,None,None]
    rational=BatchRational(rg.s0,64);inf,res=rational.fit(symbol)
    variance=VARIANCES+4*FIXED_G*np.log(Q/5.)
    inputs=np.exp(-grid.k[:,None]**2/variance)/(np.pi*variance);rootweight=np.sqrt(grid.weights);rhs=inputs*rootweight[:,None]
    fields=[inputs[selected]]
    for shift in rational.shifts:fields.append(grid.solve(shift,rhs)[selected]/rootweight[selected,None])
    coeff=np.concatenate([inf[None],res],axis=0)
    native=np.einsum('skv,sxfc,xv->kxfc',np.array(fields),coeff,weights)
    boundary=np.einsum('xv,bv->bx',weights,np.exp(-b[:,None]**2*variance/4))
    reference=np.einsum('kb,b,bx,bxfc->kxfc',j0(k[:,None]*b[None,:]),bw,boundary,symbol(b*b))
    norm=np.maximum(np.max(abs(reference.sum(axis=-1)),axis=0),1e-14)
    err=float(np.max(abs(native.sum(axis=-1)-reference.sum(axis=-1))/norm[None,:,:]))
    cerr=float(np.max(abs(native-reference)/norm[None,:,:,None]));assert err<3e-4 and cerr<3e-4
    reports.append(dict(Q=Q,mu=Q,zeta=Q*Q,x=xs.tolist(),k=k.tolist(),flavors=list(FLAVORS),components=['born','qq_one_loop','qg_one_loop'],native_components=native.tolist(),reference_components=reference.tolist(),max_scaled_error=err,max_component_scaled_error=cerr))
    print('TMD',Q,err,cerr,flush=True)
out=dict(scope='Conditional Collins/CSS2 single-TMD snapshots at mu=Q,zeta=Q^2; prior massless profile, normalized neural boundary, fixed g; no hard factor, no positivity clipping, no calibrated uncertainty',units='GeV^-2, collinear convention f rather than x*f',seed=seed,reports=reports,passed=True,script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
(ev/'single_tmd.json').write_text(json.dumps(out,indent=2)+'\n')
