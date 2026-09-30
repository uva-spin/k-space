#!/usr/bin/env python3
"""Physical-PDF fixed-point validation; independent b-space benchmark only.

No experimental normalization, Y term, nuclear target, fit to data or universal
identity for the moving capped profile is asserted by this diagnostic.
"""
import os
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):os.environ.setdefault(key,'1')
from pathlib import Path
import sys,json,argparse,time,hashlib
import numpy as np
from scipy.integrate import solve_ivp,quad,quad_vec
from scipy.special import roots_legendre,j1
from scipy.optimize import minimize_scalar
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from css2_coefficients import CF
from css2_reference import make_rg,DiagnosticProfile,PhysicalPointSymbol
from lhapdf_provider import LHAPDFProvider
from rational_radial import fit_capped_symbol,RadialGrid


def independent_block(pdf,pid,x,mu,a):
    """Adaptive integration, separately written from production GL matching."""
    f=pdf.f(pid,x,mu)
    knots=pdf.validation_xknots
    points=sorted(x/y for y in knots if x<y<1.)
    qq=CF*(-np.pi**2/6*f+quad(lambda z:2*(1-z)/z*pdf.f(pid,x/z,mu),x,1,epsabs=1e-10,epsrel=1e-9,points=points,limit=500)[0])
    qg=quad(lambda z:2*(1-z)*pdf.f(21,x/z,mu),x,1,epsabs=1e-10,epsrel=1e-9,points=points,limit=500)[0]
    return np.array([f,a*qq,a*qg])

class IndependentBspace:
    def __init__(self,symbol):
        self.symbol=symbol;c=symbol.rg.c;self.c=c
        self.coupling=solve_ivp(lambda l,y:[-2*sum(b*y[0]**(n+2) for n,b in enumerate(c.beta))],
                     (np.log(91.1876),np.log(5.)),[.1180024/(4*np.pi)],method='DOP853',
                     rtol=3e-13,atol=3e-15,dense_output=True).sol
        self.nodes,self.weights=roots_legendre(96)
    def __call__(self,s):
        mu=self.symbol.profile.mu(s);Q=self.symbol.profile.Q
        lo,hi=np.log(mu),np.log(Q)
        ts=(hi+lo)/2+(hi-lo)/2*self.nodes
        aa=self.coupling(ts)[0]
        gamma=sum(v*aa**(n+1) for n,v in enumerate(self.c.noncusp))
        cusp=sum(v*aa**(n+1) for n,v in enumerate(self.c.cusp))
        ai=self.coupling(lo)[0]
        di=sum(v*ai**(n+1) for n,v in enumerate(self.c.d))
        # -S/2 = -d_i ln(Q/mu_i) + integral[gamma -2 Gamma ln(Q/mu)].
        log_leg=-di*(hi-lo)+(hi-lo)/2*np.dot(self.weights,gamma-2*cusp*(hi-ts))
        blocks=self.symbol.cache(np.log(mu/self.symbol.rg.mu0))
        aq=self.coupling(hi)[0]
        h=aq*CF*(-16+7*np.pi**2/3)
        ans=0.
        for i,pid in enumerate(self.symbol.flavors):
            a=blocks[0,i];b=blocks[1,9-i];charge=4/9 if abs(pid) in (2,4) else 1/9
            ans+=charge*((1+h)*a[0]*b[0]+(a[1]+a[2])*b[0]+a[0]*(b[1]+b[2]))
        return float(np.exp(2*log_leg)*ans)


def continuum_bins(target,variance,edges):
    contact=target(np.inf);a,b=edges[:-1],edges[1:]
    exact_contact=contact*(np.exp(-a*a/variance)-np.exp(-b*b/variance))
    # Gaussian continuum input transform is analytic. Only this validation path
    # uses a Bessel function; the production path uses radial elliptic solves.
    integrand=lambda u:(target(u*u)-contact)*np.exp(-variance*u*u/4)*(b*j1(b*u)-a*j1(a*u))
    value,error=quad_vec(integrand,0,16.,epsabs=3e-10,epsrel=3e-10,points=[.04,.055,.1,.2,.3,.5,1.],limit=400)
    return exact_contact+value,float(error)


def main():
    p=argparse.ArgumentParser();p.add_argument('--runtime',type=Path,required=True);p.add_argument('--out',type=Path,required=True);args=p.parse_args()
    args.out.mkdir(exist_ok=True,parents=True);tests=[]
    manifest=json.loads((Path(__file__).resolve().parents[1]/'manifests/css2_pdf.json').read_text())
    verified_files={}
    for item in manifest['files']:
        if item['name'].endswith('.tar.gz'): continue
        path=args.runtime/'lhapdf/share/LHAPDF/NNPDF40_nnlo_as_01180'/item['name']
        actual=hashlib.sha256(path.read_bytes()).hexdigest()
        if actual!=item['sha256']: raise RuntimeError('Pinned PDF hash mismatch: '+str(path))
        verified_files[item['name']]=actual
    def check(name,ok,detail):tests.append({'name':name,'passed':bool(ok),'detail':detail});print(name,ok,flush=True)
    pdf=LHAPDFProvider(args.runtime/'libkt_lhapdf.dylib',args.runtime/'lhapdf/share/LHAPDF')
    datafile=args.runtime/'lhapdf/share/LHAPDF/NNPDF40_nnlo_as_01180/NNPDF40_nnlo_as_01180_0000.dat'
    # Official x-grid breakpoints guide adaptive quadrature only.
    pdf.validation_xknots=np.array([float(v) for v in datafile.read_text().split('---')[1].splitlines()[1].split()])
    rg=make_rg();profile=DiagnosticProfile();start=time.monotonic()
    symbol=PhysicalPointSymbol(pdf,rg,profile)
    print('PDF matching cache ready',time.monotonic()-start,flush=True)
    # Validate interpolation at unseen nodes and refine PDF quadrature independently.
    cache_err=quad_err=0.
    for mu in np.geomspace(5.003,19.991,83):
        direct=symbol.direct_blocks(mu,480)
        scale=np.maximum(1.,abs(direct))
        cache_err=max(cache_err,float(np.max(abs(symbol.cache(np.log(mu/5.))-direct)/scale)))
        quad_err=max(quad_err,float(np.max(abs(symbol.direct_blocks(mu,240)-direct)/scale)))
    check('PDF matching interpolation refinement',cache_err<3e-6,cache_err)
    check('PDF quadrature 240 vs 480',quad_err<3e-6,quad_err)
    adaptive_err=0.
    for mu in (5.,9.7,20.):
        a=float(rg.values(np.log(mu/5.))[0]);direct=symbol.direct_blocks(mu,480)
        for leg,x in enumerate((.1,.2)):
            for i,pid in enumerate(symbol.flavors):
                ref=independent_block(pdf,pid,x,mu,a)
                adaptive_err=max(adaptive_err,float(np.max(abs(ref-direct[leg,i])/np.maximum(1.,abs(ref)))))
    check('Independent adaptive PDF matching',adaptive_err<3e-6,adaptive_err)
    failures=[]
    for pid,x,mu in ((2,.1,1.3),(2,0.,5.),(99,.1,5.)):
        try:pdf.f(pid,x,mu)
        except ValueError:failures.append(True)
        else:failures.append(False)
    check('PDF failures are not replaced by zero',all(failures),failures)
    bspace=IndependentBspace(symbol)
    points=np.r_[0.,np.geomspace(1e-6,1e4,151),np.inf]
    refvals=np.array([bspace(s) for s in points]);got=symbol(points)
    point_err=float(np.max(abs(got-refvals)/np.maximum(1.,abs(refvals))))
    check('Independent evolution and strict assembly',point_err<1e-9,point_err)
    # Physical W structure function at fixed x,Q. Bin measure 2*pi*q*dq.
    edges=np.array([0.,.5,1.,2.,3.,4.,6.,8.]);variances=[1.,1.9,2.8]
    refs=[];qerrors=[]
    for v in variances:
        val,err=continuum_bins(bspace,v,edges);refs.append(val);qerrors.append(err)
    refs=np.stack(refs,axis=1)
    check('Independent b-space quadrature',max(qerrors)<2e-8,qerrors)
    theta=.37;weights=lambda th:np.array([th*th,2*th*(1-th),(1-th)**2])
    derivative=lambda th:np.array([2*th,2-4*th,-2+2*th])
    truth=refs@weights(theta);refgrad=refs@derivative(theta)
    rows=[];fits={};refnorm=float(np.max(abs(truth)))
    settings=[(32,40.,.1),(64,40.,.1),(96,40.,.1),(96,40.,.05),(96,40.,.025),(96,80.,.025),(96,160.,.025),(96,320.,.025),(96,640.,.025),(96,1280.,.025),(128,1280.,.025)]
    for degree,R,dk in settings:
        if degree not in fits:fits[degree]=fit_capped_symbol(symbol,rg.s0,degree)
        fit=fits[degree];grid=RadialGrid(R,round(R/dk))
        inputs=np.stack([np.exp(-grid.k**2/v)/(np.pi*v) for v in variances],axis=1)
        basis=grid.bins(grid.apply(fit,inputs),edges)
        pred=basis@weights(theta);grad=basis@derivative(theta)
        best=minimize_scalar(lambda th:float(np.sum(((basis@weights(th)-truth)/refnorm)**2)),bounds=(0.,1.),method='bounded',options={'xatol':1e-12})
        err=float(np.max(abs(pred-truth)));gerr=float(np.max(abs(grad-refgrad)))
        row={'degree':degree,'radius':R,'spacing':dk,'max_absolute_bin_error':err,
             'max_relative_bin_error':float(np.max(abs(pred-truth)/np.maximum(abs(truth),1e-10))),
             'max_bin_error_over_peak_bin':err/refnorm,'max_gradient_error':gerr,
             'theta_fit':float(best.x),'theta_bias':float(best.x-theta),
             'bins':pred.tolist(),'gradient_bins':grad.tolist(),'fit_diagnostics':fit.diagnostics}
        rows.append(row);print('Native refinement',degree,R,dk,err,gerr,best.x,flush=True)
    final=rows[-1]
    check('Native bin action relative to peak bin',final['max_bin_error_over_peak_bin']<3e-4,final['max_bin_error_over_peak_bin'])
    check('Native closure bias',abs(final['theta_bias'])<5e-4,final['theta_bias'])
    check('Native gradient relative to peak bin',final['max_gradient_error']/refnorm<3e-4,final['max_gradient_error']/refnorm)
    boxchange=float(np.max(abs(np.array(rows[-2]['bins'])-rows[-3]['bins']))/refnorm)
    degreechange=float(np.max(abs(np.array(rows[-1]['bins'])-rows[-2]['bins']))/refnorm)
    check('Box refinement',boxchange<1e-7,boxchange)
    check('Rational degree refinement',degreechange<1e-5,degreechange)
    # Analytic theta derivative against a finite difference of native predictions.
    h=1e-5;fd=basis@(weights(theta+h)-weights(theta-h))/(2*h)
    check('Native theta finite-difference gradient',np.max(abs(fd-grad))<1e-8,float(np.max(abs(fd-grad))))
    result={'scope':'Fixed pp diagnostic with real NNPDF4.0, full N3LL evolution and strict NLO insertions; no experimental cross section',
       'config':{'Q':20.,'xA':.1,'xB':.2,'pdf':'NNPDF40_nnlo_as_01180','member':0,'nf':5,'alpha_mz':.1180024,'mz':91.1876,
                 'coupling_running_loops':4,'pdf_order':'published NNLO','profile_floor':5.,'profile_bmax':1.5,'profile_width':.12,
                 'matching_coordinate':'b_eval=b0/mu_profile(s)','matching_Lb':0.,'initial_zeta':'mu_profile(s)^2',
                 'final_mu':20.,'final_zeta':400.,'NP_CS':0.,'input':'each leg theta*G(0.5)+(1-theta)*G(1.4), G(v)=exp(-k^2/v)/(pi*v)',
                 'theta':theta,'flavors':list(symbol.flavors),'electroweak':'photon charges only; no process prefactor'},
       'verified_pdf_sha256':verified_files,
       'coefficients':{key:list(getattr(rg.c,key)) for key in ('beta','cusp','noncusp','d')},
       'endpoints':{'M0':float(symbol(0.)),'Minfinity':float(symbol(np.inf)),'contact_retained_exactly':True,'uncapped_zero_constraint_used':False},
       'b_point_components':{str(s):{k:float(v) for k,v in symbol.components(s).items()} for s in (0.,.01,.04,1.)},
       'bin_edges':edges.tolist(),'reference_bins':truth.tolist(),'reference_gradients':refgrad.tolist(),
       'refinements':rows,'tests':tests,'passed':all(t['passed'] for t in tests),'elapsed_seconds':time.monotonic()-start}
    repo=Path(__file__).resolve().parents[1]
    result['code_sha256']={str(p.relative_to(repo)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted((repo/'src').glob('*')) if p.is_file()}
    result['benchmark_script_sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    result['runtime_bridge_sha256']=hashlib.sha256((args.runtime/'libkt_lhapdf.dylib').read_bytes()).hexdigest()
    result['configuration_sha256']=hashlib.sha256(json.dumps({'config':result['config'],'coefficients':result['coefficients'],'PDF_files':verified_files},sort_keys=True).encode()).hexdigest()
    (args.out/'physical_point.json').write_text(json.dumps(result,indent=2)+'\n')
    print('PASSED',result['passed'],flush=True)
    return 0 if result['passed'] else 1
if __name__=='__main__':raise SystemExit(main())
