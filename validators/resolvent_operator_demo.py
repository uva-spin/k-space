#!/usr/bin/env python3
"""Direct momentum-grid resolvent approximation; no Fourier/Hankel transforms.

A bounded analytic TEST multiplier m(r)=exp(c*r), not QCD, is evaluated via
operator de Casteljau recursion. Its same-grid reference is scipy.linalg.expm.
Neither this finite-grid experiment nor its error bound proves continuum-grid
convergence or physical N3LL-prime accuracy.
"""
from __future__ import annotations
import argparse, hashlib, importlib.metadata, json, platform
from pathlib import Path
import numpy as np
from scipy import sparse
from scipy.sparse.linalg import splu
from scipy.linalg import expm
from scipy.integrate import quad
from scipy.special import k0

ROOT=Path(__file__).resolve().parents[1]

def main()->int:
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--out',type=Path,default=ROOT/'outputs');args=p.parse_args();args.out.mkdir(parents=True,exist_ok=True)
    count=9; length=4.0; h=2*length/(count+1);tau=.6;c=.7
    k=np.linspace(-length+h,length-h,count)
    T1=sparse.diags([-np.ones(count-1),2*np.ones(count),-np.ones(count-1)],[-1,0,1],format='csc')/h**2
    T=sparse.kronsum(T1,T1,format='csc');I=sparse.eye(count**2,format='csc')
    solver=splu(tau*I+T)
    def apply_R(v): return tau*solver.solve(np.asarray(v))
    xx,yy=np.meshgrid(k,k,indexing='ij');f=np.exp(-.7*(xx**2+yy**2)).ravel()
    R=apply_R(np.eye(f.size))
    ref=expm(c*R)@f
    norm=lambda v: float(np.linalg.norm(v)*h)
    rows=[]
    # Operator-valued de Casteljau: (I-R)v_i + R v_(i+1).
    for N in [4,8,16,32,64]:
        work=f[:,None]*np.exp(c*np.arange(N+1)/N)[None,:]
        for level in range(N): work=work[:,:-1]+apply_R(work[:,1:]-work[:,:-1])
        approximate=work[:,0];error=norm(approximate-ref)
        bound=c*c*np.exp(c)/(8*N)*norm(f)
        rows.append(dict(N=N,L2_same_grid_error=error,theorem_bound=bound,passed=bool(error<=bound*(1+1e-12))))
    normalization,estimate=quad(lambda r:tau*r*k0(np.sqrt(tau)*r),0,np.inf,epsabs=2e-12,epsrel=2e-12)
    tests=[dict(name=f'Bernstein operator error bound N={r["N"]}',passed=r['passed']) for r in rows]
    tests += [dict(name='same-grid errors decrease under degree refinement',passed=all(rows[j+1]['L2_same_grid_error']<rows[j]['L2_same_grid_error'] for j in range(4))),
              dict(name='resolvent kernel normalization',passed=abs(normalization-1)<2e-11),
              dict(name='resolvent positive spectrum and contraction',passed=bool(np.all(np.linalg.eigvalsh(R)>0) and np.max(np.linalg.eigvalsh(R))<1))]
    d=dict(summary=dict(passed=all(t['passed'] for t in tests),checks=len(tests),passed_count=sum(t['passed'] for t in tests),
          physical_QCD_test=False,continuum_grid_error_established=False),tests=tests,refinement=rows,
          kernel_normalization=normalization,quadrature_error_estimate=estimate,
          setup=dict(interior_grid=[count,count],half_extent=length,spacing=h,tau=tau,c=c,multiplier='m(r)=exp(c*r)',
                     boundary='Dirichlet finite grid',reference='independent dense matrix exponential on same grid'),
          provenance=dict(script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),python=platform.python_version(),
                          packages={p:importlib.metadata.version(p) for p in ['numpy','scipy']}))
    (args.out/'resolvent_demo.json').write_text(json.dumps(d,indent=2)+'\n')
    lines=['Direct-k resolvent demonstration (nonphysical test multiplier)',f"Checks: {d['summary']['passed_count']}/{len(tests)}",
          'No Fourier or Hankel transforms are used. This is not a continuum/QCD validation.',*map(str,rows)]
    (args.out/'resolvent_demo.txt').write_text('\n'.join(lines)+'\n');print('\n'.join(lines))
    return 0 if d['summary']['passed'] else 1
if __name__=='__main__':raise SystemExit(main())
