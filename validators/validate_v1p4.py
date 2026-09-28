#!/usr/bin/env python3
"""Reproducible structural audit of the KCSS2 note (no physical NNLO/N3LO payload).

The tests distinguish exact symbolic identities, numerical weak/bin actions, and
negative controls. Passing them does NOT certify an all-order N3LL-prime solver.
Run from the bundle root: python validators/validate_v1p4.py
"""
from __future__ import annotations
import argparse, csv, hashlib, importlib.metadata, json, platform, sys
from pathlib import Path
import mpmath as mp
import sympy as s

ROOT = Path(__file__).resolve().parents[1]
mp.mp.dps = 55
TESTS: list[dict] = []
BINS: list[dict] = []

def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()

def exact(name: str, residual, kind: str = 'symbolic') -> None:
    if isinstance(residual, s.MatrixBase):
        rr = residual.applyfunc(s.expand)
        ok = all(x == 0 for x in rr)
        text = 'zero matrix' if ok else str(rr)
    else:
        rr = s.simplify(residual)
        ok, text = rr == 0, str(rr)
    TESTS.append(dict(name=name, kind=kind, passed=bool(ok), residual=text))

def condition(name: str, value: bool, detail: str, kind='counterexample') -> None:
    TESTS.append(dict(name=name, kind=kind, passed=bool(value), detail=detail))

def numeric(name: str, left, right, tolerance=mp.mpf('1e-35'), category='numerical') -> None:
    r = abs(left-right)
    TESTS.append(dict(name=name, kind=category, passed=bool(r < tolerance),
                      residual=mp.nstr(r, 18), tolerance=mp.nstr(tolerance, 8)))

def polyint(matrix, var):
    return matrix.applyfunc(lambda v: s.integrate(s.expand(v), (var, 0, var)))

def coefficient(matrix, a, n):
    return matrix.applyfunc(lambda v: s.expand(v).coeff(a,n))

def main() -> int:
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--out',type=Path,default=ROOT/'outputs')
    p.add_argument('--provenance',type=Path,default=ROOT/'manifests/source_provenance.json')
    args=p.parse_args(); args.out.mkdir(parents=True,exist_ok=True)
    eta,X,a,L,ell = s.symbols('eta X a L ell')
    z3,z5=s.zeta(3),s.zeta(5)
    logh=sum(2*s.zeta(j)*eta**j/s.Integer(j) for j in [3,5,7,9])
    h=s.series(s.exp(logh),eta,0,9).removeO()
    invh=s.series(s.exp(-logh-eta*X),eta,0,9).removeO().expand()
    P=[s.expand((-1)**n*s.factorial(n)*invh.coeff(eta,n+1)) for n in range(8)]
    expected=[-X,-X**2/2,-X**3/3-4*z3/3,-X**4/4-4*z3*X,
              -X**5/5-8*z3*X**2-48*z5/5,
              -X**6/6-s.Rational(40,3)*z3*X**3-48*z5*X-s.Rational(80,3)*z3**2]
    for n,v in enumerate(expected): exact(f'Fourier polynomial P{n}',P[n]-v)
    # A triangular inverse dictionary obtained independently by polynomial elimination.
    def project(poly):
        remainder=s.expand(poly); coeffs={}
        for n in reversed(range(8)):
            q=s.expand(remainder).coeff(X,n+1)/s.expand(P[n]).coeff(X,n+1)
            if q!=0: coeffs[n]=s.simplify(q); remainder=s.expand(remainder-q*P[n])
        coeffs[-1]=s.simplify(remainder)
        return coeffs
    J=[s.Integer(1)]+[project(X**m)[-1] for m in range(1,7)]
    target=[1,0,0,-4*z3,0,-48*z5,160*z3**2]
    for m in range(7):
        exact(f'Canonical disk moment J{m}',J[m]-target[m])
        exact(f'Generator/cumulant dictionary J{m}',J[m]-(-1)**m*s.factorial(m)*h.coeff(eta,m))
    conv={(0,0):{-1:0,1:-2},(0,1):{-1:-2*z3,2:-s.Rational(3,2)},
          (0,2):{-1:0,3:-s.Rational(4,3),0:4*z3},
          (1,1):{-1:0,3:-1,0:4*z3},
          (0,3):{-1:-12*z5,4:-s.Rational(5,4),1:12*z3},
          (1,2):{-1:-8*z5,4:-s.Rational(5,6),1:12*z3}}
    for (m,n),terms in conv.items():
        mapped=sum(v*(1 if j==-1 else P[j]) for j,v in terms.items())
        exact(f'Convolution L{m} L{n}',P[m]*P[n]-mapped)
    # Direct massless two-dimensional convolution via its convergent Feynman parameter integral.
    for e,x in [(mp.mpf('.31'),mp.mpf('.27')),(mp.mpf('.23'),mp.mpf('.41'))]:
        beta_integral=mp.quad(lambda u:u**(x-1)*(1-u)**(e-1),[0,mp.mpf('.5'),1])
        massless=mp.pi*mp.gamma(1-e-x)/(mp.gamma(1-e)*mp.gamma(1-x))*beta_integral
        normal=lambda z:mp.exp(-2*mp.euler*z)*mp.gamma(1-z)/(mp.pi*mp.gamma(z))
        numeric(f'Riesz convolution strip eta={e},xi={x}',normal(e)*normal(x)*massless,
                normal(e+x),mp.mpf('1e-12'))
    condition('Reject globally closed strip-generator algebra',1-2*(s.Rational(3,4)+s.Rational(3,4))<=-1,
              'eta=xi=3/4 are individually locally integrable, their b-space product is not.')

    # Weak profiles use d[theta(t) C(t)]/dt with t=k_T^2, NOT the derivative on t>0 alone.
    c=mp.mpf('1.7')
    numeric('Constant cumulative preserves delta contact',mp.quad(lambda t:c*mp.exp(-t),[0,mp.inf]),c)
    condition('Old origin-subtraction identity rejected on delta', c!=0,
              'Old disk identity predicts 0; correct answer for C(t)=1.7 is 1.7.')
    M=mp.mpf('2.3')
    for n in range(4):
        C=lambda t: -mp.log(M/t)**(n+1)/(n+1)
        weak=mp.quad(lambda t:C(t)*mp.exp(-t),[0,M,mp.inf])
        plus=mp.quad(lambda t:mp.log(M/t)**n*mp.expm1(-t)/t,[0,M]) + mp.quad(lambda t:mp.log(M/t)**n*mp.exp(-t)/t,[M,mp.inf])
        numeric(f'Fixed-scale weak/plus identity n={n}',weak,plus,mp.mpf('1e-40'))
        Cp=lambda t: -mp.log((M+t)/t)**(n+1)/(n+1)
        Cpder=lambda t: mp.log((M+t)/t)**n*M/(t*(M+t))
        weakp=mp.quad(lambda t:Cp(t)*mp.exp(-t),[0,1,mp.inf])
        explicit=Cp(1)+mp.quad(lambda t:Cpder(t)*mp.expm1(-t),[0,1])+mp.quad(lambda t:Cpder(t)*mp.exp(-t),[1,mp.inf])
        numeric(f'Frozen profile with divergent origin cumulative n={n}',weakp,explicit,mp.mpf('1e-39'))
    numeric('Running delta coefficient includes continuous derivative',
            mp.quad(lambda t:mp.log(M+t)*mp.exp(-t),[0,1,mp.inf]),
            mp.log(M)+mp.quad(lambda t:mp.exp(-t)/(M+t),[0,1,mp.inf]))
    K,mu=s.symbols('K mu', positive=True); n=s.symbols('n',integer=True,nonnegative=True)
    mup=s.Function('mup')(K); log=s.log(mup**2/K**2)
    exact('Profile chain-rule density',s.diff(-log**(n+1)/(n+1),K)/(2*s.pi*K)
          -log**n/(s.pi*K**2)*(1-K*s.diff(mup,K)/mup))
    exact('Canonical profile annihilates logarithmic cumulative',(-s.log(mu**2/K**2)**3/3).subs(mu,K))
    exact('Canonical profile retains frak-L2 contact',J[3]+4*z3)
    A=mp.mpf('.73'); Qp=mp.mpf('7')
    edges=list(map(mp.mpf,['0','.17','.53','1.21','2.8','5.6','9']))
    for nval in range(4):
        def Cp(t): return -(2*A*mp.exp(-t/Qp**2))**(nval+1)/(nval+1)
        def deriv(t): return (2*A)**(nval+1)/Qp**2*mp.exp(-(nval+1)*t/Qp**2)
        total=mp.mpf('0')
        for lo,hi in zip(edges[:-1],edges[1:]):
            contact=Cp(0) if lo==0 else mp.mpf('0')
            val=contact+mp.quad(deriv,[lo**2,hi**2])
            expected_val=Cp(hi**2)-(Cp(lo**2) if lo>0 else 0)
            total+=val
            numeric(f'Contact-complete bin n={nval},[{lo},{hi}]',val,expected_val)
            BINS.append(dict(n=nval,lower=str(lo),upper=str(hi),contact=mp.nstr(contact,20),
                value=mp.nstr(val,25),endpoint=mp.nstr(expected_val,25),residual=mp.nstr(val-expected_val,12)))
        numeric(f'Contact-complete bin additivity n={nval}',total,Cp(edges[-1]**2))
    Q,r=s.symbols('Q r',positive=True)
    exact('Final rapidity product',Q**2*s.exp(2*r)*Q**2*s.exp(-2*r)-Q**4)

    # Ordered flavor-matrix recurrence; Mellin-space fixture, not physical splitting functions.
    bet=[s.Rational(7,3),s.Rational(5,4),s.Rational(9,7)]
    Gam=[s.diag(4,7),s.diag(3,5),s.diag(2,6)]
    gam=[s.diag(6,8),s.diag(2,9),s.diag(5,4)]
    finiteD=[s.zeros(2),s.diag(s.Rational(2,3),s.Rational(5,7)),s.diag(s.Rational(3,8),s.Rational(4,9))]
    pmat=[s.Matrix([[1,2],[3,4]]),s.Matrix([[2,-1],[5,3]]),s.Matrix([[4,3],[-2,7]])]
    cmat=[s.eye(2),s.Matrix([[1,2],[-3,1]]),s.Matrix([[2,-2],[1,4]]),s.Matrix([[3,1],[2,-1]])]
    D=[s.zeros(2)]
    for j in range(1,4):
        rhs=Gam[j-1]+sum(((j-r)*bet[r-1]*D[j-r] for r in range(1,j)),s.zeros(2))
        D.append(finiteD[j-1]+polyint(rhs,L))
    bar=[s.eye(2)]
    for j in range(1,4):
        rhs=sum(((gam[r-1]-D[r])*bar[j-r]-bar[j-r]*pmat[r-1] for r in range(1,j+1)),s.zeros(2))
        rhs+=sum((2*(j-r)*bet[r-1]*bar[j-r] for r in range(1,j)),s.zeros(2))
        bar.append(cmat[j]+polyint(rhs/2,L))
    Ds=sum((a**j*D[j] for j in range(1,4)),s.zeros(2))
    rapid=s.eye(2)-ell*Ds/2+ell**2*(Ds*Ds)/8-ell**3*(Ds*Ds*Ds)/48
    Cs0=sum((a**j*bar[j] for j in range(4)),s.zeros(2))
    Cs=sum((a**j*coefficient(rapid*Cs0,a,j) for j in range(4)),s.zeros(2))
    beta=-2*sum(bet[r]*a**(r+2) for r in range(3))
    gamma=sum((a**j*(gam[j-1]-ell*Gam[j-1]) for j in range(1,4)),s.zeros(2))
    ps=sum((a**j*pmat[j-1] for j in range(1,4)),s.zeros(2))
    rg=2*Cs.diff(L)-2*Cs.diff(ell)+beta*Cs.diff(a)-gamma*Cs+Cs*ps
    cs=2*Cs.diff(ell)+Ds*Cs
    for j in range(1,4):
        exact(f'Ordered matrix mu-RG order {j}',coefficient(rg,a,j))
        exact(f'Ordered matrix CS order {j}',coefficient(cs,a,j))
        exact(f'Ordered matrix finite boundary order {j}',coefficient(Cs,a,j).subs({L:0,ell:0})-cmat[j])
        degree=max(s.degree(s.expand(v),L) for v in coefficient(Cs,a,j))
        condition(f'Matching polynomial degree bound order {j}',degree<=2*j,f'maximum L degree={degree}',kind='symbolic')
    condition('Reversed flavor ordering rejected',any(s.expand(v)!=0 for v in cmat[1]*pmat[0]-pmat[0]*cmat[1]),
              'C1 P0 != P0 C1 for the stated noncommuting Mellin fixture.')
    (args.out/'matrix_fixture.json').write_text(json.dumps(dict(
        notice='Nonphysical rational Mellin-moment fixture; does not validate physical finite NNLO/N3LO matching.',
        beta=[str(v) for v in bet],c=[str(v) for v in cmat],Pmu=[str(v) for v in pmat],
        Gamma=[str(v) for v in Gam],gamma=[str(v) for v in gam],
        reconstructed=[str(coefficient(Cs,a,j)) for j in range(4)]),indent=2)+'\n')

    # General scalar D3 polynomial and missing all-order integrability.
    b0,b1,b2,G0,G1,G2,G3,d2,d3=s.symbols('b0 b1 b2 G0 G1 G2 G3 d2 d3')
    Dpoly=a*G0*L+a**2*(b0*G0*L**2/2+G1*L+d2)+a**3*(b0**2*G0*L**3/3+(b0*G1+b1*G0/2)*L**2+(G2+2*b0*d2)*L+d3)
    drg=s.expand(2*s.diff(Dpoly,L)-2*(b0*a**2+b1*a**3+b2*a**4)*s.diff(Dpoly,a)-2*(a*G0+a**2*G1+a**3*G2+a**4*G3))
    for j in range(1,4): exact(f'D3 RG coefficient order {j}',drg.coeff(a,j))
    condition('D3 is not an exact four-loop-cusp RG solution',drg.coeff(a,4)!=0,str(drg.coeff(a,4)))
    t=s.symbols('t',real=True); GG=s.Function('Gamma')(t); integ=s.Function('AGamma')(t)
    exact('Integrated RG completion derivative',s.diff(s.Symbol('Dref')+2*integ,t).subs(s.diff(integ,t),GG)-2*GG)
    # Positivity of leading cubic coefficient is the analytic nonexistence test.
    u,cpos=s.symbols('u c',positive=True)
    condition('Reject raw cubic exponent as regular tempered kernel',s.limit((cpos*u**3-2*u)/u,u,s.oo)==s.oo,
              'At b -> 0 with rho>0, radial logarithm of integrand behaves as c*u^3-2*u.')
    hidden=s.series(s.exp(s.Symbol('c')*a**4*L**5),a,0,4).removeO()
    exact('FO3 cannot see altered all-order logarithm',hidden-1)
    exact('Hidden a4 L5 term survives resummed scaling',(a**4*L**5).subs(L,1/a)-1/a)

    # Complete two-leg hard/matching product through cubic order.
    hs=s.symbols('h0:4'); As=s.symbols('A0:4'); Bs=s.symbols('B0:4')
    series=s.expand(sum(hs[i]*a**i for i in range(4))*sum(As[i]*a**i for i in range(4))*sum(Bs[i]*a**i for i in range(4)))
    for j in range(4):
        finite=sum(hs[r]*As[v]*Bs[j-r-v] for r in range(j+1) for v in range(j-r+1))
        exact(f'All ordered two-leg product partitions order {j}',series.coeff(a,j)-finite)
    Wpert,NP,FO=s.symbols('Wpert NP FO')
    exact('Matching leaves NP unexpanded and avoids subtraction loss',(Wpert+NP)+(FO-Wpert)-(FO+NP))

    # Explicit nonperturbative module: forward heat evolution, no b-space numerics.
    x,y,tau,g=s.symbols('x y tau g',positive=True)
    heat=s.exp(-(x*x+y*y)/(4*g*tau))/(4*s.pi*g*tau)
    exact('NP heat equation',s.diff(heat,tau)-g*(s.diff(heat,x,2)+s.diff(heat,y,2)))
    gv=mp.mpf('.2'); tv=mp.mpf('.7')
    norm=mp.quad(lambda tt:mp.exp(-tt/(4*gv*tv))/(4*gv*tv),[0,1,mp.inf])
    moment=mp.quad(lambda tt:tt*mp.exp(-tt/(4*gv*tv))/(4*gv*tv),[0,1,mp.inf])
    numeric('NP heat kernel normalized',norm,1)
    numeric('NP heat kernel second moment',moment,4*gv*tv)
    s1,s2=mp.mpf('.3'),mp.mpf('.8'); q=mp.mpf('.7')
    # Cartesian separability reduces a direct 2D convolution to two 1D quadratures.
    gaussian1=lambda xx,tt:mp.exp(-xx**2/(4*gv*tt))/mp.sqrt(4*mp.pi*gv*tt)
    convx=mp.quad(lambda v:gaussian1(v,s1)*gaussian1(q-v,s2),[-mp.inf,0,q,mp.inf])
    convy=mp.quad(lambda v:gaussian1(v,s1)*gaussian1(-v,s2),[-mp.inf,0,mp.inf])
    numeric('NP heat semigroup direct momentum convolution',convx*convy,
            mp.exp(-q*q/(4*gv*(s1+s2)))/(4*mp.pi*gv*(s1+s2)))

    good=sum(t['passed'] for t in TESTS)
    report=dict(schema='kcss-structural-audit/1.0',tests=TESTS,
        summary=dict(passed=good==len(TESTS),number_passed=good,number_tests=len(TESTS),
            physical_three_loop_coefficients_tested=False,
            all_order_accuracy_certified_by_this_finite_test_run=False,
            numerical_extraction_performed=False),
        provenance=dict(provenance_manifest=args.provenance.name,provenance_manifest_sha256=sha(args.provenance),script_sha256=sha(Path(__file__)),
            python=sys.version,platform=platform.platform(),
            packages={k:importlib.metadata.version(k) for k in ['sympy','mpmath']},mpmath_decimal_digits=mp.mp.dps),
        canonical_disk_moments=[str(j) for j in J])
    (args.out/'audit_results.json').write_text(json.dumps(report,indent=2)+'\n')
    with (args.out/'audit_tests.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=['name','kind','passed','residual','tolerance','detail']);w.writeheader();w.writerows(TESTS)
    with (args.out/'contact_complete_bins.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(BINS[0]));w.writeheader();w.writerows(BINS)
    binmax=max(abs(mp.mpf(row['residual'])) for row in BINS)
    lines=['KCSS2 substantive structural audit','='*40,f'Checks passed: {good}/{len(TESTS)}',
           f'Maximum contact-complete bin residual: {mp.nstr(binmax,16)}',
           'Precision: 55 decimal digits; tolerances are recorded per test.',
           'Physical NNLO/N3LO payload checked: False',
           'All-order accuracy certified by this finite test run: False','',
           *[f"{'PASS' if t['passed'] else 'FAIL'} [{t['kind']}] {t['name']} :: {t.get('residual',t.get('detail',''))}" for t in TESTS]]
    (args.out/'audit_report.txt').write_text('\n'.join(lines)+'\n')
    print('\n'.join(lines[:7]))
    for test in TESTS:
        if not test['passed']: print('FAILED:',test)
    return 0 if good==len(TESTS) else 1

if __name__=='__main__': raise SystemExit(main())
