#!/usr/bin/env python3
"""E866 pp forward validation. No fit to experimental values is performed."""
from pathlib import Path
import sys,json,argparse,time,hashlib
import numpy as np
from scipy.integrate import solve_ivp,quad
from scipy.special import j0
from scipy.optimize import least_squares
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from e866_forward import *
from lhapdf_provider import LHAPDFProvider
from _runtime_paths import add_pdf_arguments,pdf_paths

class Reference:
    def __init__(self,engine,n=32,nqt=16):
        self.engine=engine;self.c=engine.rg.c;self.n=n;self.nqt=nqt
        c=self.c
        sol=solve_ivp(lambda t,y:[-2*sum(b*y[0]**(i+2) for i,b in enumerate(c.beta))],(np.log(91.1876),np.log(engine.rg.mu0)),[.1180024/(4*np.pi)],method='DOP853',rtol=3e-13,atol=3e-15,dense_output=True)
        assert sol.success;self.alpha=sol.sol
        b=[];w=[]
        for lo,hi in zip([0,.08,.16,.25,.5,1,2,4,8],[.08,.16,.25,.5,1,2,4,8,16]):
            bn,wn=nodes_weights(lo,hi,n);b.extend(bn);w.extend(wn)
        self.b=np.array(b);self.bw=np.array(w);self.symbols=[]
        for Q,pref,ids,symbol,fits in engine.nodes:self.symbols.append(self.components(symbol,self.b**2))
    def components(self,symbol,s):
        c=self.c;Q=symbol.profile.Q;mu=symbol.profile.mu(s);low=np.log(mu);hi=np.log(Q)
        z,w=legendre(48);t=low[...,None]+(hi-low[...,None])*(z+1)/2
        a=self.alpha(t.reshape(-1))[0].reshape(t.shape)
        gamma=sum(v*a**(n+1) for n,v in enumerate(c.noncusp));cusp=sum(v*a**(n+1) for n,v in enumerate(c.cusp))
        ai=self.alpha(low.ravel())[0].reshape(low.shape);d=sum(v*ai**(n+1) for n,v in enumerate(c.d))
        logleg=-d*(hi-low)+(hi-low)/2*np.sum(w*(gamma-2*cusp*(hi-t)),axis=-1)
        return symbol.cache(np.log(mu/symbol.rg.mu0))*np.exp(2*logleg)[...,None]
    def evaluate(self,model=GaussianCSInput(),mode='flat_pt'):
        output=np.zeros((len(self.engine.rows),6,3));b=self.b
        for item,components in zip(self.engine.nodes,self.symbols):
            Q,pref,ids,symbol,fits=item;v=model.variance(Q)
            f=np.exp(-v*b*b/4);d=-b*b/4*f
            cols=np.stack([f,2*d,8*np.log(Q/model.Q0)*d],axis=-1)
            for i in ids:
                r=self.engine.rows[i];q,w=nodes_weights(r['qtlow'],r['qthigh'],self.nqt)
                if mode=='flat_pt':w/=r['qthigh']-r['qtlow']
                elif mode=='annular':w*=2*q/(r['qthigh']**2-r['qtlow']**2)
                else:raise ValueError(mode)
                kernel=w@j0(q[:,None]*b[None,:])
                output[i]+=pref*np.einsum('b,bc,bp->cp',self.bw*b/(2*np.pi)*kernel,components,cols)
        return output

def errors(got,ref):
    norm=np.maximum(abs(ref[:,:,0].sum(axis=1)),1e-14)
    return dict(component_over_total=float(np.max(abs(got[:,:,0]-ref[:,:,0])/norm[:,None])),total_relative=float(np.max(abs(got[:,:,0].sum(axis=1)-ref[:,:,0].sum(axis=1))/norm)),gradient_scaled=float(np.max(abs(got[:,:,1:].sum(axis=1)-ref[:,:,1:].sum(axis=1))/np.maximum(np.max(abs(ref[:,:,1:].sum(axis=1)),axis=1)[:,None],1e-14))))

def main():
    p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True);p.add_argument('--quick',action='store_true');add_pdf_arguments(p);args=p.parse_args();args.out.mkdir(parents=True,exist_ok=True)
    root=Path(__file__).resolve().parents[1];runs=root.parent.parent;start=time.monotonic()
    bridge,data=pdf_paths(args,root)
    pdf=LHAPDFProvider(bridge,data)
    pin=json.loads((root/'manifests/css2_pdf.json').read_text())
    for r in pin['files']:
        if not r['name'].endswith('.tar.gz'):assert hashlib.sha256((data/'NNPDF40_nnlo_as_01180'/r['name']).read_bytes()).hexdigest()==r['sha256']
    manifest=json.loads((root/'manifests/e866_pp_bins.json').read_text());rows=[r for r in manifest['rows'] if r['selected_core']]
    tests=[]
    def check(name,ok,detail):tests.append(dict(name=name,passed=bool(ok),detail=detail));print('CHECK',name,ok,detail,flush=True)
    progress=lambda group,n:print('BUILD',group,n,flush=True)
    engine=NativeE866(pdf,rows,nQ=4,nx=8,nmu=64,order=160,degree=64,radius=640.,spacing=.00625,progress=progress)
    print('EVALUATE',time.monotonic()-start,flush=True)
    native=engine.evaluate();ref_engine=Reference(engine);ref=ref_engine.evaluate();err=errors(native,ref)
    print('INITIAL ERRORS',err,flush=True)
    check('native bins',err['total_relative']<3e-4,err)
    check('separate insertions',err['component_over_total']<3e-4,err['component_over_total'])
    check('input and CS gradients',err['gradient_scaled']<5e-4,err['gradient_scaled'])
    ref_fine=Reference(engine,n=48,nqt=24).evaluate();re=errors(ref,ref_fine)
    check('independent Bessel quadrature',re['total_relative']<2e-8,re)
    cache_err=0.;quad_err=0.;point_err=0.
    for _,_,_,sym,_ in engine.nodes[::max(1,len(engine.nodes)//5)]:
        for mu in np.geomspace(float(sym.profile.mu(np.inf))*1.00001,sym.profile.Q*.99999,5):
            direct=sym.direct_luminosity(mu,320);low=sym.direct_luminosity(mu,160)
            norm=max(abs(direct.sum()),1e-12)
            cache_err=max(cache_err,float(np.max(abs(sym.cache(np.log(mu/sym.rg.mu0))-direct))/norm));quad_err=max(quad_err,float(np.max(abs(low-direct))/norm))
        s=np.r_[0.,np.geomspace(1e-6,100,21),np.inf]
        rv=ref_engine.components(sym,s);point_err=max(point_err,float(np.max(abs(sym.components(s)-rv))/max(np.max(abs(rv)),1e-12)))
    # Independently coded adaptive canonical convolutions at new E866 x values.
    from benchmark_css2_point import independent_block
    dat=data/'NNPDF40_nnlo_as_01180/NNPDF40_nnlo_as_01180_0000.dat'
    pdf.validation_xknots=np.array([float(v) for v in dat.read_text().split('---')[1].splitlines()[1].split()])
    ae=0.
    for mu in (5.1,12.5):
        for x in (.09,.22,.55):
            blocks=canonical_blocks(pdf,engine.rg,np.array([x]),mu,240)[0]
            a=float(engine.rg.values(np.log(mu/engine.rg.mu0))[0])
            for i,pid in enumerate(FLAVORS):
                expected=independent_block(pdf,pid,x,mu,a)
                ae=max(ae,float(np.max(abs(blocks[i]-expected)/np.maximum(1.,abs(expected)))))
    check('independent adaptive PDF convolution',ae<3e-6,ae)
    check('PDF matching cache',cache_err<3e-6,cache_err);check('PDF convolution order',quad_err<3e-6,quad_err);check('independent RG',point_err<2e-9,point_err)
    outputs={'schema':'E866-pp-forward-validation/1','scope':'Massless nf5, leading-power photon W-only, explicit shared Gaussian boundary and quadratic NP-CS. Not a complete massive prediction or experimental fit.','settings':engine.settings,'parameters':dict(lambda0=.7,g=.05,Q0=5.),'rows':[], 'tests':tests}
    for i,r in enumerate(rows):outputs['rows'].append(dict(id=r['id'],Qlow=r['Qlow'],Qhigh=r['Qhigh'],xflow=r['xflow'],xfhigh=r['xfhigh'],qtlow=r['qtlow'],qthigh=r['qthigh'],data=r['value'],stat=r['stat'],syst=r['syst'],native_components=dict(zip(COMPONENTS,native[i,:,0].tolist())),reference_components=dict(zip(COMPONENTS,ref[i,:,0].tolist())),native_total=float(native[i,:,0].sum()),reference_total=float(ref[i,:,0].sum()),gradients=native[i,:,1:].sum(axis=0).tolist()))
    outputs['errors']=err
    outputs['observable_sensitivity']={'annular_vs_flat_max_relative':float(np.max(abs(ref_engine.evaluate(mode='annular')[:,:,0].sum(axis=1)/ref[:,:,0].sum(axis=1)-1))),'tabulated_mean_energy_discrepancy_unresolved':True}
    if not args.quick:
        # All tolerances fixed before running the refinements.
        refined=NativeE866(pdf,rows,nQ=6,nx=12,nmu=96,order=240,degree=96,radius=1280.,spacing=.003125,progress=progress)
        rf=Reference(refined,n=48,nqt=24).evaluate();fine=refined.evaluate();fe=errors(fine,rf);long=errors(ref,rf)
        check('mass and xF quadrature refinement',long['total_relative']<3e-4,long)
        check('refined native bins',fe['total_relative']<8e-5,fe)
        check('refined gradients',fe['gradient_scaled']<1.5e-4,fe['gradient_scaled'])
        outputs['refined_errors']=fe;outputs['longitudinal_refinement']=long;outputs['refined_settings']=refined.settings
        # Isolate radius/degree with unchanged physical quadrature and cache.
        from rational_radial import RadialGrid
        oldgrid=engine.grid
        engine.grid=RadialGrid(1280.,round(1280./.00625));box=engine.evaluate();engine.grid=oldgrid
        check('box refinement',errors(box,native)['total_relative']<1e-7,errors(box,native))
        oldnodes=engine.nodes;newnodes=[]
        for Q,pref,ids,sym,fits in oldnodes:
            cache={}
            def comp(s,c):
                key=(np.asarray(s).shape,np.asarray(s).tobytes())
                if key not in cache:cache[key]=sym.components(s)
                return cache[key][...,c]
            fs=[fit_capped_symbol(lambda s,c=c:comp(s,c),engine.rg.s0,96) for c in range(6)]
            newnodes.append((Q,pref,ids,sym,fs))
        engine.nodes=newnodes;degree=engine.evaluate();engine.nodes=oldnodes
        check('degree refinement',errors(degree,native)['total_relative']<1e-5,errors(degree,native))
        truth=ref[:,:,0].sum(axis=1);scale=np.sqrt(np.array([r['stat']**2+r['syst']**2 for r in rows]));train=np.array([r['xflow']<0 for r in rows]);test=~train
        last={}
        def fun_and_jac(par):
            key=tuple(par)
            if key not in last:
                ans=engine.evaluate(GaussianCSInput(*par));last.clear();last[key]=(ans[:,:,0].sum(axis=1),ans[:,:,1:].sum(axis=1))
            return last[key]
        fit=least_squares(lambda p:(fun_and_jac(p)[0][train]-truth[train])/scale[train],[.55,.08],jac=lambda p:fun_and_jac(p)[1][train]/scale[train,None],bounds=([.2,0],[1.5,.2]),xtol=1e-10,ftol=1e-10,gtol=1e-10,max_nfev=20)
        pred,jac=fun_and_jac(fit.x);sv=np.linalg.svd(jac[train]/scale[train,None],compute_uv=False)
        closure=dict(parameters=fit.x.tolist(),truth=[.7,.05],bias=(fit.x-[.7,.05]).tolist(),success=bool(fit.success),evaluations=fit.nfev,whitened_singular_values=sv.tolist(),condition=float(sv[0]/sv[-1]),heldout_max_relative=float(np.max(abs(pred[test]/truth[test]-1))),training_rows=int(sum(train)),heldout_rows=int(sum(test)),scope='Two-parameter closure conditional on x/flavor-independent Gaussian input, fixed perturbative prescription, fixed normalization. No claim of general TMD identifiability.')
        check('conditional multiscale closure',fit.success and abs(fit.x[0]-.7)<1e-3 and abs(fit.x[1]-.05)<5e-4,closure)
        check('heldout xF predictions',closure['heldout_max_relative']<3e-4,closure['heldout_max_relative'])
        outputs['closure']=closure
    outputs['passed']=all(t['passed'] for t in tests);outputs['elapsed_seconds']=time.monotonic()-start
    outputs['code_sha256']={str(p.relative_to(root)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted((root/'src').glob('*')) if p.is_file()}
    outputs['script_sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest();outputs['data_manifest_sha256']=hashlib.sha256((root/'manifests/e866_pp_bins.json').read_bytes()).hexdigest()
    (args.out/'e866_forward.json').write_text(json.dumps(outputs,indent=2)+'\n');print('PASSED',outputs['passed'],time.monotonic()-start,flush=True)
    return 0 if outputs['passed'] else 1
if __name__=='__main__':raise SystemExit(main())
