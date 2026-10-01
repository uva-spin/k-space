"""General-scale CS boundary and explicitly scoped noncanonical profiles.

Derived from dD/dln(mu)=2 Gamma and da/dln(mu)=-2 sum beta_n a^(n+2).
This is a massless fixed-flavor implementation. A threshold crossing is rejected;
changing nf alone is not massive-TMD matching.
"""
from dataclasses import dataclass
import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import CubicSpline
from css2_coefficients import qcd_coefficients,CF
from css2_matching import convolutions,hard_nlo
from css2_reference import DiagnosticProfile
from resummed_operator import B0,CanonicalRG


def cs_boundary(a,L,c):
    if len(c.cusp)<3 or len(c.beta)<2 or len(c.d)<3 or c.d[0]!=0:
        raise ValueError('Need three-loop CS boundary and canonical d1=0')
    a,L=np.broadcast_arrays(np.asarray(a,float),np.asarray(L,float))
    if not np.isfinite(a).all() or not np.isfinite(L).all() or (a<0).any():raise ValueError('Invalid CS arguments')
    b0,b1=c.beta[:2];g0,g1,g2=c.cusp[:3];d2,d3=c.d[1:3]
    return (a*g0*L+a*a*(b0*g0*L*L/2+g1*L+d2)
            +a**3*(b0*b0*g0*L**3/3+(b1*g0/2+b0*g1)*L*L+(g2+2*b0*d2)*L+d3))


@dataclass(frozen=True)
class FixedFlavorWindow:
    nf: int
    lower: float
    upper: float
    pdf_min: float=1.65
    def validate(self,lo,hi):
        if not np.isfinite([lo,hi,self.lower,self.upper,self.pdf_min]).all() or self.nf not in (4,5):
            raise ValueError('Invalid supported flavor window')
        if lo<self.pdf_min-1e-12:raise ValueError('Profile below the pinned PDF domain')
        if lo<self.lower-1e-12 or hi>=self.upper or hi<lo:
            raise ValueError('Threshold crossing requires massive-TMD matching; fixed-flavor mode refuses it')


def make_fixed_rg(window,mu0,anchor_mu,anchor_alpha,Q):
    window.validate(min(mu0,anchor_mu,Q),max(mu0,anchor_mu,Q))
    if not np.isfinite(anchor_alpha) or min(mu0,anchor_mu,anchor_alpha)<=0:raise ValueError('Positive coupling and scales required')
    c=qcd_coefficients(window.nf)
    sol=solve_ivp(lambda t,y:[-2*c.series(y[0],c.beta,start=2)],(np.log(anchor_mu),np.log(mu0)),
                  [anchor_alpha/(4*np.pi)],method='DOP853',rtol=2e-13,atol=2e-15)
    if not sol.success or sol.y[0,-1]<=0:raise RuntimeError('Coupling integration failed')
    return CanonicalRG(c,mu0,float(sol.y[0,-1]),tmax=np.log(window.upper/mu0))


@dataclass(frozen=True)
class NoncanonicalProfile(DiagnosticProfile):
    power: float=16.
    C5: float=1.
    def __post_init__(self):
        super().__post_init__()
        if not np.isfinite([self.power,self.C5]).all() or self.power<=1 or self.C5<=0:
            raise ValueError('Invalid matching coordinate profile')
    def b_eval(self,s):
        s=np.asarray(s,float)
        if np.isnan(s).any() or (s<0).any():raise ValueError('Require s>=0')
        bs=np.zeros(s.shape);sel=s>0
        bs[sel]=1/np.sqrt(1/s[sel]+1/self.bmax**2)
        bmin=B0/(self.C5*self.Q)
        largest=np.maximum(bs,bmin)
        return largest*((bs/largest)**self.power+(bmin/largest)**self.power)**(1/self.power)
    def Lb(self,s):return 2*np.log(self.mu(s)*self.b_eval(s)/B0)


class BoundaryEvolution:
    """Universal two-scale evolution at a FIXED b and fixed boundary anchor.

    Re-anchoring by re-expanding the boundary has a finite-order remainder.
    Transporting this same anchor gives exact numerical path independence.
    """
    def __init__(self,rg,b,mu_anchor):
        self.rg=rg;self.b=np.asarray(b,float);self.t0=np.log(np.asarray(mu_anchor,float)/rg.mu0)
        a,g,_,_=rg.values(self.t0)
        self.g0=g
        self.D0=cs_boundary(a,2*np.log(np.asarray(mu_anchor)*self.b/B0),rg.c)
    def D(self,mu):
        g=self.rg.values(np.log(np.asarray(mu)/self.rg.mu0))[1]
        return self.D0+2*(g-self.g0)
    def log_u(self,target,initial):
        rg=self.rg
        tf,wf=rg.endpoint(*target);ti,wi=rg.endpoint(*initial)
        _,gf,nf,jf=rg.values(tf);_,gi,ni,ji=rg.values(ti)
        return nf-ni-2*wi*(gf-gi)+2*(jf-ji)-(wf-wi)*self.D(target[0])
    def log_leg(self,mu_initial,Q):
        rg=self.rg;ti=np.log(np.asarray(mu_initial)/rg.mu0);tf=np.log(Q/rg.mu0)
        _,gi,ni,ji=rg.values(ti);_,gf,nf,jf=rg.values(tf)
        return -self.D(mu_initial)*(tf-ti)+nf-ni-2*tf*(gf-gi)+2*(jf-ji)


class NoncanonicalPointSymbol:
    def __init__(self,pdf,rg,profile,window,xA=.1,xB=.2,nodes=192,order=240):
        window.validate(float(profile.mu(np.inf)),profile.Q)
        # Distinct PDF/evolution flavor counts cannot be silently combined.
        expected=qcd_coefficients(window.nf)
        for key in ('beta','cusp','noncusp','d'):
            actual,wanted=getattr(rg.c,key),getattr(expected,key)
            if len(actual)!=len(wanted) or not np.allclose(actual,wanted,rtol=1e-12,atol=1e-12):
                raise ValueError('RG coefficients do not match declared flavor window')
        if float(profile.mu(np.inf))<rg.mu0-1e-12:
            raise ValueError('RG interval must include the actual composed profile endpoint')
        self.pdf,self.rg,self.profile=pdf,rg,profile
        self.window,self.xA,self.xB,self.order=window,xA,xB,order
        self.flavors=tuple(range(-window.nf,0))+tuple(range(1,window.nf+1))
        self.ts=np.linspace(0,np.log(profile.Q/rg.mu0),nodes)
        self.blocks=np.array([self.direct_kernels(rg.mu0*np.exp(t)) for t in self.ts])
        self.cache=CubicSpline(self.ts,self.blocks,axis=0,extrapolate=False)
    def direct_kernels(self,mu):
        return np.array([[convolutions(x,lambda y:self.pdf.f(pid,y,mu),lambda y:self.pdf.f(21,y,mu),self.order)
                          for pid in self.flavors] for x in (self.xA,self.xB)])
    def components(self,s):
        mu=self.profile.mu(s);t=np.log(mu/self.rg.mu0);L=self.profile.Lb(s)
        kernels=self.cache(t);a=self.rg.values(t)[0]
        born,pqq,pqg,cqq,cqg=np.moveaxis(kernels,-1,0)
        qq=a[...,None,None]*(cqq-2*L[...,None,None]*pqq+CF*(3*L-L*L)[...,None,None]*born)
        qg=a[...,None,None]*(cqg-2*L[...,None,None]*pqg)
        blocks=np.stack([born,qq,qg],axis=-1)
        A=blocks[...,0,:,:];B=blocks[...,1,::-1,:]
        charges=np.array([4/9 if abs(pid) in (2,4) else 1/9 for pid in self.flavors])
        lum=np.sum(charges*A[...,0]*B[...,0],axis=-1)
        aq=float(self.rg.values(np.log(self.profile.Q/self.rg.mu0))[0])
        evo=BoundaryEvolution(self.rg,self.profile.b_eval(s),mu)
        ev=np.exp(2*evo.log_leg(mu,self.profile.Q))
        terms={'born':lum,'hard_1':hard_nlo(aq,self.profile.Q,self.profile.Q)*lum,
               'ope_qq_leg_a':np.sum(charges*A[...,1]*B[...,0],axis=-1),
               'ope_qg_leg_a':np.sum(charges*A[...,2]*B[...,0],axis=-1),
               'ope_qq_leg_b':np.sum(charges*A[...,0]*B[...,1],axis=-1),
               'ope_qg_leg_b':np.sum(charges*A[...,0]*B[...,2],axis=-1)}
        terms={k:ev*v for k,v in terms.items()};terms['strict_total']=sum(terms.values())
        return terms
    def __call__(self,s):return self.components(s)['strict_total']
