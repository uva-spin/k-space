#!/usr/bin/env python3
"""Fixed-nf noncanonical operator checks; no threshold or data-fit claim."""
from pathlib import Path
import argparse,sys,json,hashlib
import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp
from scipy.optimize import minimize_scalar
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from css2_noncanonical import NoncanonicalProfile,NoncanonicalPointSymbol,FixedFlavorWindow,make_fixed_rg
from css2_reference import make_rg
from css2_coefficients import CF
from lhapdf_provider import LHAPDFProvider
from rational_radial import RadialGrid,fit_capped_symbol
from benchmark_css2_point import continuum_bins


def recurrence_D(c):
    # Independent coefficient recurrence, not the production closed polynomial.
    L=sp.Symbol('L');polys=[sp.Float(c.d[0])+sp.Float(c.cusp[0])*L]
    for n in (2,3):
        rhs=sp.Float(c.cusp[n-1])
        for i in range(n-1):rhs+=c.beta[i]*(n-i-1)*polys[n-i-2]
        polys.append(sp.Float(c.d[n-1])+sp.integrate(rhs,(L,0,L)))
    fs=[sp.lambdify(L,p,'numpy') for p in polys]
    return lambda a,l:sum(a**(i+1)*f(l) for i,f in enumerate(fs))


class ODEReference:
    def __init__(self,symbol,anchor_mu,alpha):
        self.symbol=symbol;self.c=symbol.rg.c;self.db=recurrence_D(self.c)
        self.alpha=solve_ivp(lambda t,y:[-2*self.c.series(y[0],self.c.beta,start=2)],
           (np.log(anchor_mu),np.log(symbol.rg.mu0)),[alpha/(4*np.pi)],method='DOP853',
           rtol=3e-13,atol=3e-15,dense_output=True).sol
        # Both use the same defined input data; evolution is integrated independently.
    def __call__(self,s):
        sym=self.symbol;mu=float(sym.profile.mu(s));Q=sym.profile.Q
        a0=float(self.alpha(np.log(mu))[0]);L=float(sym.profile.Lb(s));d0=self.db(a0,L)
        def rhs(t,y):
            a=y[0]
            return [-2*self.c.series(a,self.c.beta,start=2),
                    self.c.series(a,self.c.noncusp)-2*self.c.series(a,self.c.cusp)*(np.log(Q)-t)]
        sol=solve_ivp(rhs,(np.log(mu),np.log(Q)),[a0,-d0*np.log(Q/mu)],
                      method='RK45',rtol=2e-10,atol=2e-12)
        if not sol.success:raise RuntimeError(sol.message)
        aq=sol.y[0,-1];v2=np.exp(2*sol.y[1,-1]);h=aq*CF*(-16+7*np.pi**2/3)
        kernels=sym.cache(np.log(mu/sym.rg.mu0));parts=[]
        for leg in kernels:
            parts.append([(b,a0*(cq-2*L*pq+CF*(3*L-L*L)*b),a0*(cg-2*L*pg)) for b,pq,pg,cq,cg in leg])
        total=0
        for i,pid in enumerate(sym.flavors):
            a=parts[0][i];b=parts[1][-i-1];charge=4/9 if abs(pid) in (2,4) else 1/9
            total+=charge*((1+h)*a[0]*b[0]+(a[1]+a[2])*b[0]+a[0]*(b[1]+b[2]))
        return v2*total


def main():
    p=argparse.ArgumentParser();p.add_argument('--previous-run',type=Path,required=True);p.add_argument('--runtime',type=Path,required=True);p.add_argument('--out',type=Path,required=True);args=p.parse_args()
    args.out.mkdir(parents=True,exist_ok=True)
    data=args.previous_run/'runtime/lhapdf/share/LHAPDF';pdf=LHAPDFProvider(args.runtime/'libkt_lhapdf.dylib',data)
    manifest=json.loads((Path(__file__).resolve().parents[1]/'manifests/css2_pdf.json').read_text())
    for row in manifest['files']:
        if not row['name'].endswith('.tar.gz'):
            assert hashlib.sha256((data/'NNPDF40_nnlo_as_01180'/row['name']).read_bytes()).hexdigest()==row['sha256']
    cases=[]
    for nf,Q,floor,anchor,alpha in [(4,4.5,1.65,4.5,pdf.alpha(4.5)),(5,7.5,5.,91.1876,.1180024)]:
        print('START',nf,Q,floor,flush=True)
        window=FixedFlavorWindow(nf,1.65 if nf==4 else 4.92,4.92 if nf==4 else 172.5)
        profile=NoncanonicalProfile(Q=Q,floor=floor)
        actual_lower=float(profile.mu(np.inf))
        rg=make_fixed_rg(window,actual_lower,anchor,alpha,Q)
        sym=NoncanonicalPointSymbol(pdf,rg,profile,window,nodes=256,order=320)
        ref=ODEReference(sym,anchor,alpha)
        points=np.r_[0.,np.geomspace(1e-5,1e4,61),np.inf]
        got=sym(points);expected=np.array([ref(x) for x in points]);point_error=float(np.max(abs(got-expected)/np.maximum(1.,abs(expected))))
        interpolation=0.
        for mu in np.geomspace(floor*1.0001,Q*.9999,47):
            direct=sym.direct_kernels(mu)
            interpolation=max(interpolation,float(np.max(abs(sym.cache(np.log(mu/rg.mu0))-direct)/np.maximum(1.,abs(direct)))))
        edges=np.array([0.,.5,1.,2.,3.,4.]);vs=[1.,1.9,2.8]
        references=[continuum_bins(ref,v,edges) for v in vs];basisref=np.stack([r[0] for r in references],axis=1)
        weight=lambda t:np.array([t*t,2*t*(1-t),(1-t)**2]);derivative=lambda t:np.array([2*t,2-4*t,-2+2*t])
        theta=.37;truth=basisref@weight(theta);gradref=basisref@derivative(theta);norm=max(abs(truth))
        rows=[];fits={}
        for degree,R,dk in [(64,320.,.05),(96,320.,.05),(96,320.,.025),(96,640.,.025),(128,640.,.025)]:
            if degree not in fits:fits[degree]=fit_capped_symbol(sym,rg.s0,degree)
            grid=RadialGrid(R,round(R/dk));fit=fits[degree]
            inputs=np.stack([np.exp(-grid.k**2/v)/(np.pi*v) for v in vs],axis=1)
            basis=grid.bins(grid.apply(fit,inputs),edges);pred=basis@weight(theta);gradient=basis@derivative(theta)
            best=minimize_scalar(lambda t:np.sum((basis@weight(t)-truth)**2),bounds=(0,1),method='bounded',options={'xatol':1e-12})
            row={'degree':degree,'radius':R,'spacing':dk,'bins':pred.tolist(),'bin_error_over_peak':float(max(abs(pred-truth))/norm),
                 'gradient_error_over_peak':float(max(abs(gradient-gradref))/norm),'theta_bias':float(best.x-theta),'fit_diagnostics':fit.diagnostics}
            rows.append(row);print('REFINE',nf,degree,R,dk,row['bin_error_over_peak'],row['theta_bias'],flush=True)
        box=float(max(abs(np.array(rows[-2]['bins'])-rows[-3]['bins']))/norm)
        degree=float(max(abs(np.array(rows[-1]['bins'])-rows[-2]['bins']))/norm)
        tests={'independent_ODE':point_error<2e-8,'matching_cache':interpolation<3e-6,
               'b_quadrature':max(r[1] for r in references)<2e-8,'bins':rows[-1]['bin_error_over_peak']<3e-4,
               'gradient':rows[-1]['gradient_error_over_peak']<3e-4,'closure':abs(rows[-1]['theta_bias'])<5e-4,
               'box':box<1e-7,'degree':degree<1e-5}
        cases.append({'nf':nf,'Q':Q,'floor':floor,'actual_minimum_mu':actual_lower,'alpha_anchor_mu':anchor,'alpha_anchor':alpha,'scope':'fixed-nf massless diagnostic with physical PDFs, not a massive-quark prediction',
          'matching_coordinate':'power16 bOPE','CS_boundary':'general Lb through a^3','NP_CS':0.,'point_error':point_error,'interpolation_error':interpolation,
          'endpoints':[float(sym(0.)),float(sym(np.inf))],'Lb_endpoints':[float(profile.Lb(0.)),float(profile.Lb(np.inf))],
          'reference_bins':truth.tolist(),'box_change':box,'degree_change':degree,'refinements':rows,'tests':tests,'passed':all(tests.values())})
    result={'cases':cases,'passed':all(c['passed'] for c in cases)}
    repo=Path(__file__).resolve().parents[1];result['code_sha256']={str(f.relative_to(repo)):hashlib.sha256(f.read_bytes()).hexdigest() for f in sorted((repo/'src').glob('*')) if f.is_file()}
    (args.out/'noncanonical_benchmark.json').write_text(json.dumps(result,indent=2)+'\n')
    print('PASSED',result['passed'],flush=True)
    return 0 if result['passed'] else 1
if __name__=='__main__':raise SystemExit(main())
