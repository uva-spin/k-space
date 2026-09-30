"""A declared capped diagnostic prescription composed with universal CSS2 RG.

The Q-dependent b remapping is a process prescription. RG/CS derivatives are
with b_eval (and its boundary) held fixed. No universal-TMD claim is made for
its derivative when Q is also used to move the remapping.
"""
from dataclasses import dataclass
import numpy as np
from scipy.integrate import solve_ivp
from resummed_operator import B0,CanonicalRG
from css2_coefficients import qcd_coefficients

@dataclass(frozen=True)
class DiagnosticProfile:
    Q: float=20.
    floor: float=5.
    bmax: float=1.5
    width: float=.12
    def __post_init__(self):
        if not np.isfinite([self.Q,self.floor,self.bmax,self.width]).all() or min(self.floor,self.bmax,self.width)<=0 or self.Q<=self.floor:
            raise ValueError('Require finite Q>floor>0, bmax>0, width>0')
    def mu(self,s):
        s=np.asarray(s,float)
        if np.isnan(s).any() or (s<0).any(): raise ValueError('Require s>=0')
        raw=np.full(s.shape,np.inf);sel=s>0
        raw[sel]=B0*np.sqrt(1/s[sel]+1/self.bmax**2)
        lower=self.floor+self.width*np.logaddexp(0.,(raw-self.floor)/self.width)
        capped=self.Q-self.width*np.logaddexp(0.,(self.Q-lower)/self.width)
        return capped
    def b_eval(self,s): return B0/self.mu(s)


def make_rg(mu0=5.,alpha_mz=.1180024,mz=91.1876):
    c=qcd_coefficients(5)
    # Explicit four-loop nf=5 alpha_s anchored at the PDF set's stated MZ.
    # PDF interpolation itself remains the published NNLO set. Its tabulated
    # three-loop alpha_s is NOT relabeled as four-loop evolution.
    sol=solve_ivp(lambda t,y: [-2*c.series(y[0],c.beta,start=2)],
                  (np.log(mz),np.log(mu0)),[alpha_mz/(4*np.pi)],
                  method='DOP853',rtol=2e-13,atol=2e-15)
    if not sol.success or sol.y[0,-1]<=0: raise RuntimeError('Coupling integration failed')
    return CanonicalRG(c,mu0,float(sol.y[0,-1]),tmax=8.)


def profile_leg(rg,profile,s):
    tc=np.log(profile.mu(s)/rg.mu0)
    return np.exp(rg.log_v_tc(tc,profile.Q,profile.Q**2))

class PhysicalPointSymbol:
    """Fixed pp luminosity, all five quark flavors, flavor-independent NP input.

    Caches only the longitudinal PDF/matching blocks versus ln(mu). Fitted
    transverse inputs never undergo a transform. The cache can be independently
    refined and checked against direct official-LHAPDF evaluations.
    """
    flavors=(-5,-4,-3,-2,-1,1,2,3,4,5)
    def __init__(self,pdf,rg,profile,xA=.1,xB=.2,nodes=192,order=240):
        from scipy.interpolate import CubicSpline
        self.pdf,self.rg,self.profile=pdf,rg,profile
        self.xA,self.xB,self.order=xA,xB,order
        self.ts=np.linspace(0,np.log(profile.Q/rg.mu0),nodes)
        self.blocks=np.array([self.direct_blocks(rg.mu0*np.exp(t)) for t in self.ts])
        self.cache=CubicSpline(self.ts,self.blocks,axis=0,extrapolate=False)
    def direct_blocks(self,mu,order=None):
        from css2_matching import match_nlo
        a=float(self.rg.values(np.log(mu/self.rg.mu0))[0])
        return np.array([[[m.born,m.qq,m.qg] for pid in self.flavors
                          for m in [match_nlo(x,lambda y:self.pdf.f(pid,y,mu),
                                              lambda y:self.pdf.f(21,y,mu),a,
                                              order=order or self.order)]] for x in (self.xA,self.xB)])
    def components(self,s):
        from css2_matching import hard_nlo
        mu=self.profile.mu(s);t=np.log(mu/self.rg.mu0)
        blocks=self.cache(t)
        a=blocks[...,0,:,:];b=blocks[...,1,::-1,:] # reversal maps pid -> -pid
        charges=np.array([4/9 if abs(pid) in (2,4) else 1/9 for pid in self.flavors])
        born=np.sum(charges*a[...,0]*b[...,0],axis=-1)
        aq=float(self.rg.values(np.log(self.profile.Q/self.rg.mu0))[0])
        ev=profile_leg(self.rg,self.profile,s)**2
        vals={'born':born,'hard_1':hard_nlo(aq,self.profile.Q,self.profile.Q)*born,
              'ope_qq_leg_a':np.sum(charges*a[...,1]*b[...,0],axis=-1),
              'ope_qg_leg_a':np.sum(charges*a[...,2]*b[...,0],axis=-1),
              'ope_qq_leg_b':np.sum(charges*a[...,0]*b[...,1],axis=-1),
              'ope_qg_leg_b':np.sum(charges*a[...,0]*b[...,2],axis=-1)}
        vals={k:ev*v for k,v in vals.items()};vals['strict_total']=sum(vals.values())
        return vals
    def __call__(self,s): return self.components(s)['strict_total']
