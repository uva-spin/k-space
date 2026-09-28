#!/usr/bin/env python3
"""Reproducible structural checks of the v1.5 unexpanded construction.

These checks do not certify a QCD N3LL-prime spectrum. Formal higher-order
coefficients are marked fixtures. The numerical engine uses the explicitly
listed lower-order QCD coefficients from coefficients_demo.json.
"""
from __future__ import annotations
import argparse, csv, hashlib, json, platform, sys
from pathlib import Path
import numpy as np
import scipy
from scipy.integrate import quad
from scipy.linalg import eigh
from scipy.sparse import diags, eye, kron
from scipy.special import k0
import sympy as sp
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from resummed_operator import (Coefficients, CanonicalRG, B0,
                              apply_resolvent_polynomial, bernstein_scalar)

TESTS=[]
def check(name, passed, category, residual=None, tolerance=None, details=''):
    TESTS.append(dict(name=name,passed=bool(passed),category=category,
                      residual=None if residual is None else float(residual),
                      tolerance=tolerance,details=details))
    if not passed: print('FAILED:',name,details,flush=True)

def mul(A,B,N):
    out=[sp.S(0)]*(N+1)
    for n in range(N+1):
        out[n]=sp.expand(sum(A[j]*B[n-j] for j in range(n+1)))
    return out

def power(A,k,N):
    out=[sp.S(0)]*(N+1);out[0]=sp.S(1)
    for _ in range(k):out=mul(out,A,N)
    return out

def compose(A,coeff,N,start=1):
    out=[sp.S(0)]*(N+1)
    for j,c in enumerate(coeff):
        P=power(A,start+j,N)
        out=[sp.expand(x+c*y) for x,y in zip(out,P)]
    return out

def running(beta,N,t):
    A=[sp.S(0)]*(N+1);A[1]=sp.S(1)
    for n in range(2,N+1):
        val=compose(A,beta,n,start=2)[n]
        A[n]=sp.integrate(2*val,t)
    return A

def exponent(A,G,g,d,N,t):
    cu=compose(A,G,N);nc=compose(A,g,N);bd=compose(A,d,N)
    return [sp.expand(sp.integrate(nc[n]-2*t*cu[n],t)-t*bd[n]) for n in range(N+1)]

def symbolic_checks(out):
    t,a=sp.symbols('t a')
    b0,b1,G0,G1,G2,d2,d3=sp.symbols('b0 b1 G0 G1 G2 d2 d3')
    A=running([b0,b1],3,t)
    cu=compose(A,[G0,G1,G2],3);bd=compose(A,[0,d2,d3],3)
    D=[sp.expand(2*sp.integrate(cu[n],t)+bd[n]) for n in range(4)]
    expected=[0,2*G0*t,2*b0*G0*t**2+2*G1*t+d2,
       sp.Rational(8,3)*b0**2*G0*t**3+(4*b0*G1+2*b1*G0)*t**2
       +2*(G2+2*b0*d2)*t+d3]
    for n in range(1,4):
        check(f'canonical_D_expansion_a{n}',sp.expand(D[n]-expected[n])==0,'symbolic')
    check('D3_missing_running_boundary_term_rejected',
          sp.expand(D[3]-(expected[3]-4*b0*d2*t))!=0,'negative_control')
    # Higher orders below are distinct rational fixtures, not physical QCD data.
    N=9
    beta=[sp.S(x) for x in [7,11,13,17,19]]
    G=[sp.S(x) for x in [3,5,7,11,13]]
    g=[sp.S(x) for x in [2,3,5,7]]
    d=[sp.S(x) for x in [0,2,3,11]]
    Ar=running(beta[:4],N,t);Af=running(beta,N,t)
    Er=exponent(Ar,G[:4],g[:3],d[:3],N,t)
    variants={
      'omitted_cusp5':exponent(Ar,G,g[:3],d[:3],N,t),
      'omitted_noncusp4':exponent(Ar,G[:4],g,d[:3],N,t),
      'omitted_CS_boundary4':exponent(Ar,G[:4],g[:3],d,N,t),
      'omitted_beta5':exponent(Af,G[:4],g[:3],d[:3],N,t),
      'all_omissions':exponent(Af,G,g,d,N,t)}
    rows=[]
    for label,E in variants.items():
        for n in range(1,N+1):
            dif=sp.expand(E[n]-Er[n])
            degree=-1 if dif==0 else int(sp.degree(dif,t))
            passed=dif==0 if n<=3 else degree<=n-3
            check(f'{label}_a{n}_tower_bound',passed,'formal_fixture',
                  details=f'log degree={degree}, allowed <= {n-3}; rational fixture')
            rows.append(dict(variation=label,order=n,degree=degree,
                             difference=str(dif),allowed_degree=n-3))
    bf=[sp.S(x) for x in [1,2,3,5]]
    br=compose(Ar,bf[:3],N);bt=compose(Ar,bf,N)
    for n in range(1,N+1):
        dif=sp.expand(bt[n]-br[n]);degree=-1 if dif==0 else int(sp.degree(dif,t))
        check(f'first_omitted_prime_boundary_a{n}',dif==0 if n<=3 else degree<=n-4,
              'formal_fixture',details=f'log degree {degree}; fixture boundary')
    all1=running([b0],7,t)
    ell=exponent(all1,[G0],[],[],7,t)
    for n in range(1,8):
        exact=-G0*2**n*b0**(n-1)*t**(n+1)/sp.Integer(n+1)
        check(f'LL_infinite_tower_a{n}',sp.expand(ell[n]-exact)==0,'symbolic')
    # A coefficient-through-third-order test alone misses an LL change.
    perturb=sp.expand(a**4*t**5)
    check('finite_order_only_certificate_rejected',
          sp.series(perturb,a,0,4).removeO()==0 and sp.degree(perturb.coeff(a,4),t)==5,
          'negative_control')
    with (out/'log_tower_tests.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    (out/'symbolic_D_expansion.txt').write_text('\n'.join(f'a^{n}: {D[n]}' for n in range(1,4))+'\n')


def numerical_checks(out):
    coeff=Coefficients.load(ROOT/'data/coefficients_demo.json')
    rg=CanonicalRG(coeff)
    check('AF_beta_branch',rg.beta_min>0,'lower_order_QCD_demo',details=str(rg.beta_min))
    check('demo_does_not_claim_N3LLprime',not coeff.physical_n3llprime_payload,'scope')
    s=np.geomspace(1e-14,3.,45)
    pa=(3.,25.); pb=(8.,121.); pc=(35.,2500.)
    Lba=rg.log_u(s,pb,pa);Lcb=rg.log_u(s,pc,pb);Lca=rg.log_u(s,pc,pa)
    residual=np.max(abs(Lba+Lcb-Lca))
    check('endpoint_semigroup_log',residual<2e-10,'lower_order_QCD_demo',residual,2e-10)
    tc=rg.tc(s)
    ratio=rg.log_v_tc(tc,*pc)-rg.log_v_tc(tc,*pa)
    residual=np.max(abs(ratio-Lca))
    check('canonical_ratio_endpoint',residual<2e-10,'lower_order_QCD_demo',residual,2e-10)
    residual=np.max(abs(rg.log_u(s,pa,pa)))
    check('endpoint_identity',residual<2e-12,'lower_order_QCD_demo',residual,2e-12)
    residual=np.max(abs(rg.log_u(s,pa,pc)+Lca))
    check('endpoint_inverse',residual<2e-10,'lower_order_QCD_demo',residual,2e-10)
    # Independent scalar quadrature for the rectangle's mu leg and D difference.
    ta,wa=rg.endpoint(*pa);tb,wb=rg.endpoint(*pc)
    def fun(t):
        av=rg.values(t)[0]
        return coeff.series(av,coeff.noncusp)-2*(wa-t)*coeff.series(av,coeff.cusp)
    direct=quad(fun,ta,tb,epsabs=1e-11,epsrel=1e-11)[0]
    ldirect=direct-(wb-wa)*rg.D(s,pc[0])
    residual=np.max(abs(ldirect-Lca))
    check('independent_mu_leg_quadrature',residual<1e-9,'lower_order_QCD_demo',residual,1e-9)
    eps=1e-5
    dg=(rg.D(s,pc[0]*np.exp(eps))-rg.D(s,pc[0]*np.exp(-eps)))/(2*eps)
    value=2*coeff.series(rg.values(np.log(pc[0]/rg.mu0))[0],coeff.cusp)
    residual=np.max(abs(dg-value))
    check('resummed_D_mu_equation',residual<2e-8,'lower_order_QCD_demo',residual,2e-8)
    upr=rg.log_u(s,(pc[0],pc[1]*np.exp(2*eps)),pa)
    umr=rg.log_u(s,(pc[0],pc[1]*np.exp(-2*eps)),pa)
    residual=np.max(abs((upr-umr)/(2*eps)+rg.D(s,pc[0])))
    check('rapidity_endpoint_equation',residual<2e-8,'lower_order_QCD_demo',residual,2e-8)
    upl=rg.log_u(s,(pc[0]*np.exp(eps),pc[1]),pa)
    uml=rg.log_u(s,(pc[0]*np.exp(-eps),pc[1]),pa)
    av=rg.values(tb)[0]
    gamma=coeff.series(av,coeff.noncusp)-coeff.series(av,coeff.cusp)*np.log(pc[1]/pc[0]**2)
    residual=np.max(abs((upl-uml)/(2*eps)-gamma))
    check('mu_endpoint_equation',residual<2e-8,'lower_order_QCD_demo',residual,2e-8)
    check('spectral_canonical_region',np.allclose(rg.tc(np.array([rg.s0/8,rg.s0/4])),
          -0.5*np.log(np.array([1/8,1/4])),rtol=0,atol=1e-14),'numerical')
    check('IR_completion_is_constant',np.all(rg.tc(np.array([rg.s0,2*rg.s0,20*rg.s0]))==0),'numerical')
    uv_t=np.array([100.,300.,1000.])
    slope=rg.log_v_tc(uv_t,20.,400.)/uv_t
    predicted=-coeff.cusp[0]/coeff.beta[0]
    check('canonical_UV_log_slope_smoke',abs(slope[-1]-predicted)<0.035,
          'lower_order_QCD_demo',abs(slope[-1]-predicted),0.035,
          'Large finite t diagnostic; analytic asymptotic proof is in note')
    v=rg.v(np.r_[0.,s],20.,400.)
    check('bounded_canonical_endpoint_sample',v[0]==0 and np.isfinite(v).all(),'numerical')
    # Check closed one-loop exponent against independent scalar quadrature.
    aQ=.012;L=2.3;b0=coeff.beta[0];G0=coeff.cusp[0]
    lam=2*b0*aQ*L
    closed=G0/(2*b0*b0*aQ)*(np.log1p(-lam)+lam)
    integ=quad(lambda t:-2*t*G0*aQ/(1-2*b0*aQ*t),0,L,epsabs=1e-13)[0]
    check('LL_closed_unexpanded_vs_integral',abs(closed-integ)<2e-13,
          'numerical',abs(closed-integ),2e-13)
    tau=rg.s0
    kerint=quad(lambda k:tau*k*k0(np.sqrt(tau)*k),0,np.inf,epsabs=1e-11)[0]
    check('momentum_resolvent_kernel_normalization',abs(kerint-1)<2e-10,
          'numerical',abs(kerint-1),2e-10)
    pl=quad(lambda z:(tau/(tau+z))**2,0,np.inf,epsabs=1e-11)[0]
    check('continuum_contact_Plancherel_constant',abs(pl-tau)<2e-11,
          'numerical',abs(pl-tau),2e-11)
    # A genuine finite momentum grid; no FFT in construction or evaluator.
    n=17;halfwidth=8.;dk=2*halfwidth/(n+1)
    T1=diags([-np.ones(n-1),2*np.ones(n),-np.ones(n-1)],[-1,0,1],format='csc')/dk**2
    T=kron(T1,eye(n,format='csc'))+kron(eye(n,format='csc'),T1)
    grid=np.linspace(-halfwidth+dk,halfwidth-dk,n)
    X,Y=np.meshgrid(grid,grid,indexing='ij')
    f=np.exp(-(X*X+Y*Y)/2).ravel();f/=np.linalg.norm(f)*dk
    # This independent dense matrix oracle is used only for the test reference.
    eig,V=eigh(T.toarray())
    exact=V@(rg.v(eig,20.,400.)*(V.T@f))
    rows=[]
    for degree in [8,16,32,64,128]:
        minf,h=rg.h_nodes(degree,tau,20.,400.)
        approx=apply_resolvent_polynomial(T,f,tau,minf,h)
        r=tau/(tau+eig)
        ps=minf+r*bernstein_scalar(h,r)
        matrix_poly=V@(ps*(V.T@f))
        alg_err=np.linalg.norm(approx-matrix_poly)*dk
        check(f'resolvent_vs_matrix_polynomial_N{degree}',alg_err<2e-10,
              'same_grid_demo',alg_err,2e-10)
        actual=np.linalg.norm(approx-exact)*dk
        samegrid_bound=float(np.max(abs(ps-rg.v(eig,20.,400.))))
        check(f'same_grid_spectral_norm_bound_N{degree}',actual<=samegrid_bound+2e-10,
              'same_grid_demo',actual,samegrid_bound+2e-10,
              'Bound for this finite matrix, not a sampled continuum certificate')
        rows.append(dict(degree=degree,L2_error=actual,finite_matrix_bound=samegrid_bound,
                         recurrence_oracle_residual=alg_err))
    identity=apply_resolvent_polynomial(T,f,tau,.73,np.zeros(8))
    residual=np.linalg.norm(identity-.73*f)*dk
    check('constant_contact_is_preserved',residual<1e-14,'contact',residual,1e-14)
    check('dropped_constant_contact_is_rejected',np.linalg.norm(.73*f)*dk>.7,
          'negative_control')
    with (out/'resolvent_grid_results.csv').open('w',newline='') as fcsv:
        w=csv.DictWriter(fcsv,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    (out/'numeric_scope.json').write_text(json.dumps({
      'description':coeff.description,'physical_n3llprime_payload':False,
      'grid':{'n_each_axis':n,'halfwidth_GeV':halfwidth,'dk_GeV':dk,'boundary':'Dirichlet'},
      'uv_log_slope':slope.tolist(),'predicted_UV_log_slope':predicted,
      'continuum_grid_error_certified':False,'production_evaluator_uses_FFT':False,
      'reference':'Independent diagonalization of the same finite momentum-grid Laplacian'},indent=2)+'\n')


def main():
    pa=argparse.ArgumentParser();pa.add_argument('--out',type=Path,default=ROOT/'outputs')
    args=pa.parse_args();out=args.out;out.mkdir(parents=True,exist_ok=True)
    symbolic_checks(out);numerical_checks(out)
    counts={}
    for item in TESTS:
        q=counts.setdefault(item['category'],{'passed':0,'total':0})
        q['total']+=1;q['passed']+=item['passed']
    hashes={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()
       for p in [Path(__file__),ROOT/'src/resummed_operator.py',ROOT/'data/coefficients_demo.json',ROOT/'manifests/source_provenance.json']}
    report={'all_checks_passed':all(t['passed'] for t in TESTS),
            'passed':sum(t['passed'] for t in TESTS),'total':len(TESTS),
            'categories':counts,'tests':TESTS,
            'environment':{'python':platform.python_version(),'numpy':np.__version__,
                           'scipy':scipy.__version__,'sympy':sp.__version__},
            'sha256':hashes,'physical_N3LLprime_prediction_validated':False,
            'note':'Finite formal/symbolic and lower-order numerical checks supplement, not replace, the written all-order proofs.'}
    (out/'unexpanded_validation.json').write_text(json.dumps(report,indent=2)+'\n')
    with (out/'unexpanded_tests.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(TESTS[0]));w.writeheader();w.writerows(TESTS)
    text=['v1.5 unexpanded-operator validation',
      f"all_checks_passed: {report['all_checks_passed']}",
      f"checks: {report['passed']}/{report['total']}",
      'physical_N3LLprime_prediction_validated: False',
      'continuum_grid_error_certified: False','']
    text += [f"{k}: {v['passed']}/{v['total']}" for k,v in counts.items()]
    text += ['','Details: unexpanded_validation.json, unexpanded_tests.csv,',
             'log_tower_tests.csv, resolvent_grid_results.csv, numeric_scope.json',
             '', report['note']]
    (out/'run_summary.txt').write_text('\n'.join(text)+'\n')
    print('\n'.join(text))
    return 0 if report['all_checks_passed'] else 1
if __name__=='__main__':raise SystemExit(main())
