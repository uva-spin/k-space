"""Collins/CSS2 one-loop OPE and insertion-aware DY assembly.

PDF callback returns f(x,mu), not x*f(x,mu). Couplings are explicit inputs.
Pqq = CF [(1+z^2)/(1-z)]_+ (whole-function plus); df/dln(mu)=4a P*f.
Finite canonical constants: Collins--Rogers (74a,b); logs follow (75),(76).
"""
from dataclasses import dataclass
import numpy as np
from scipy.special import roots_legendre
from css2_coefficients import CF, TF

@dataclass(frozen=True)
class MatchingBlock:
    born: float
    qq: float
    qg: float


def _nodes(x, order):
    if not np.isfinite(x) or not 0 < x < 1 or order < 16:
        raise ValueError('Require 0<x<1 and quadrature order>=16')
    t,w=roots_legendre(order)
    return x+(1-x)*(t+1)/2, w*(1-x)/2


def convolutions(x, quark, gluon, order=160):
    z,w=_nodes(x,order)
    fx=float(quark(x))
    q=np.asarray(quark(x/z),dtype=float)/z
    g=np.asarray(gluon(x/z),dtype=float)/z
    if not (np.isfinite(fx) and np.isfinite(q).all() and np.isfinite(g).all()):
        raise ValueError('Nonfinite PDF data')
    # Endpoint subtraction is applied to the complete test function f(x/z)/z.
    primitive=-x*x/2-x-2*np.log1p(-x)
    pqq=CF*(np.dot(w,(1+z*z)/(1-z)*(q-fx))-fx*primitive)
    pqg=TF*np.dot(w,(z*z+(1-z)**2)*g)
    cqq=CF*(-np.pi**2*fx/6+2*np.dot(w,(1-z)*q))
    cqg=4*TF*np.dot(w,z*(1-z)*g)
    return fx,pqq,pqg,cqq,cqg


def match_nlo(x,quark,gluon,a,Lb=0.,ell=0.,order=160):
    """Return LO and actual O(a) insertions. ell=ln(zeta/mu^2)."""
    if not np.isfinite([a,Lb,ell]).all() or a<0:
        raise ValueError('Invalid matching parameters')
    fx,pqq,pqg,cqq,cqg=convolutions(x,quark,gluon,order)
    qq=cqq-2*Lb*pqq+CF*((3-2*ell)*Lb-Lb*Lb)*fx
    qg=cqg-2*Lb*pqg
    return MatchingBlock(fx,a*qq,a*qg)


def hard_nlo(a,Q,mu):
    if not np.isfinite([a,Q,mu]).all() or a<0 or min(Q,mu)<=0:
        raise ValueError('Invalid hard-factor parameters')
    T=2*np.log(Q/mu)
    return a*CF*(-16+7*np.pi**2/3+6*T-2*T*T)


def strict_nlo(A,B,delta_h,evolution=1.):
    """Each return component has epsilon degree <=1; evolution is unexpanded."""
    terms={'born':A.born*B.born,'hard_1':delta_h*A.born*B.born,
           'ope_qq_leg_a':A.qq*B.born,'ope_qg_leg_a':A.qg*B.born,
           'ope_qq_leg_b':A.born*B.qq,'ope_qg_leg_b':A.born*B.qg}
    terms={k:evolution*v for k,v in terms.items()}
    terms['strict_total']=sum(terms.values())
    return terms
