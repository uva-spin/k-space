"""Declared leading-power, massless nf=5 pp baseline for E866 bin studies.

This is a new prescription, not a reproduction of the earlier p-Cu fit.
PDF/perturbative matching scales stay above 4.92 GeV. Long-distance input and
CS kernel are explicit model parameters, not claimed to be perturbative.
"""
from dataclasses import dataclass
from functools import lru_cache
import numpy as np
from scipy.special import roots_legendre
from scipy.interpolate import CubicSpline
from css2_coefficients import CF
from css2_matching import hard_nlo
from css2_noncanonical import FixedFlavorWindow,make_fixed_rg
from css2_reference import profile_leg
from resummed_operator import B0
from rational_radial import fit_capped_symbol,RadialGrid
from dy_measurement import HBARC2_PB

COMPONENTS=('born','hard_1','ope_qq_leg_a','ope_qg_leg_a','ope_qq_leg_b','ope_qg_leg_b')
FLAVORS=(-5,-4,-3,-2,-1,1,2,3,4,5)
CHARGES=np.array([4/9 if abs(p) in (2,4) else 1/9 for p in FLAVORS])

@lru_cache(maxsize=32)
def legendre(n):return roots_legendre(n)

def nodes_weights(lo,hi,n):
    x,w=legendre(n);return (hi+lo)/2+(hi-lo)/2*x,(hi-lo)/2*w

def parton_x(Q,xf,sqrts=38.8):
    root=np.sqrt(np.asarray(xf)**2+4*(Q/sqrts)**2)
    return (root+xf)/2,(root-xf)/2

@dataclass(frozen=True)
class HighScaleProfile:
    Q:float
    match_scale:float=5.
    power:float=16.
    def __post_init__(self):
        if not np.isfinite([self.Q,self.match_scale,self.power]).all() or self.Q<=self.match_scale or self.match_scale<=0 or self.power<=1:raise ValueError('Invalid profile')
        FixedFlavorWindow(5,4.92,172.5).validate(float(self.mu(np.inf)),self.Q)
    def b_eval(self,s):
        s=np.asarray(s,float)
        if np.isnan(s).any() or (s<0).any():raise ValueError('Require nonnegative spectral s=b^2')
        bs=np.zeros(s.shape);sel=s>0
        bs[sel]=1/np.sqrt(1/s[sel]+(self.match_scale/B0)**2)
        low=B0/self.Q;large=np.maximum(bs,low)
        return large*((bs/large)**self.power+(low/large)**self.power)**(1/self.power)
    def mu(self,s):return B0/self.b_eval(s)


def canonical_blocks(pdf,rg,x,mu,order=160):
    """Vectorized independent-x canonical finite OPE kernels, f rather than xf."""
    x=np.asarray(x,float)
    if not np.isfinite(x).all() or (x<=0).any() or (x>=1).any():raise ValueError('PDF fraction outside (0,1)')
    t,w=legendre(order);z=x[...,None]+(1-x[...,None])*(t+1)/2
    weights=w*(1-x[...,None])/2
    g=pdf.f(21,x[...,None]/z,mu)
    cqg=2*np.sum(weights*(1-z)*g,axis=-1)
    coupling=float(rg.values(np.log(mu/rg.mu0))[0]);out=[]
    for pid in FLAVORS:
        f=pdf.f(pid,x,mu);q=pdf.f(pid,x[...,None]/z,mu)
        cqq=CF*(-np.pi**2*f/6+2*np.sum(weights*(1-z)/z*q,axis=-1))
        out.append(np.stack([f,coupling*cqq,coupling*cqg],axis=-1))
    return np.stack(out,axis=-2)

class XFAveragedSymbol:
    """Averages luminosity over xF before the linear transverse operator.

    This commutation requires the declared x/flavor-independent transverse input.
    It must not be reused unchanged for an x- or flavor-dependent input model.
    """
    def __init__(self,pdf,rg,Q,xflo,xfhi,*,nx=8,nmu=64,order=160,match_scale=5.):
        if xfhi<=xflo:raise ValueError('Invalid xF interval')
        self.pdf,self.rg=pdf,rg;self.profile=HighScaleProfile(Q,match_scale)
        self.xf,self.xw=nodes_weights(xflo,xfhi,nx);self.xw/=xfhi-xflo
        self.xa,self.xb=parton_x(Q,self.xf);self.order=order
        # Include floating-point evaluation of both exact profile endpoints.
        # This pads only the interpolation grid; no scale/profile value is clipped.
        pad=16*np.finfo(float).eps
        self.t=np.linspace(np.log(float(self.profile.mu(np.inf))/rg.mu0)-pad,np.log(Q/rg.mu0)+pad,nmu)
        values=np.array([self.direct_luminosity(rg.mu0*np.exp(t)) for t in self.t])
        self.cache=CubicSpline(self.t,values,axis=0,extrapolate=False)
    def direct_luminosity(self,mu,order=None):
        a=canonical_blocks(self.pdf,self.rg,self.xa,mu,order or self.order)
        b=canonical_blocks(self.pdf,self.rg,self.xb,mu,order or self.order)[...,::-1,:]
        h=hard_nlo(float(self.rg.values(np.log(self.profile.Q/self.rg.mu0))[0]),self.profile.Q,self.profile.Q)
        base=a[...,0]*b[...,0]
        terms=np.stack([base,h*base,a[...,1]*b[...,0],a[...,2]*b[...,0],a[...,0]*b[...,1],a[...,0]*b[...,2]],axis=-1)
        return np.einsum('x,f,xfc->c',self.xw,CHARGES,terms)
    def components(self,s):
        t=np.log(self.profile.mu(s)/self.rg.mu0)
        return self.cache(t)*np.asarray(profile_leg(self.rg,self.profile,s)**2)[...,None]
    def __call__(self,s):return self.components(s).sum(axis=-1)

@dataclass(frozen=True)
class GaussianCSInput:
    """Each leg has G(k;lambda0) at sqrt(zeta)=Q0 and D_NP=g*b^2.

    Two legs give v_pair=2*lambda0+8*g*ln(Q/Q0). This is a declared shared,
    x/flavor-independent model. g is held universal; b denotes the physical
    coordinate, not the perturbative b_eval. Parameters are in GeV^2.
    """
    lambda0:float=.7
    g:float=.05
    Q0:float=5.
    def variance(self,Q):
        v=2*self.lambda0+8*self.g*np.log(Q/self.Q0)
        if not np.isfinite([v,self.lambda0,self.g,self.Q0]).all() or self.lambda0<=0 or self.g<0 or self.Q0<=0 or v<=0:raise ValueError('Invalid NP parameters')
        return v
    def native_columns(self,k,Q):
        v=self.variance(Q);f=np.exp(-k*k/v)/(np.pi*v);dv=f*(k*k/v-1)/v
        return np.stack([f,2*dv,8*np.log(Q/self.Q0)*dv],axis=-1)


def pt_average(grid,values,lo,hi,mode='flat_pt'):
    """Average a differential invariant; not an integrated annular cross section."""
    if hi<=lo or lo<0 or hi>grid.radius:raise ValueError('Invalid qT interval')
    a=np.maximum(lo,grid.edges[:-1]);b=np.minimum(hi,grid.edges[1:]);b=np.maximum(a,b)
    if mode=='flat_pt':w=(b-a)/(hi-lo)
    elif mode=='annular':w=(b*b-a*a)/(hi*hi-lo*lo)
    else:raise ValueError('Explicit flat_pt or annular average required')
    return np.einsum('k,k...->...',w,values)

class NativeE866:
    def __init__(self,pdf,rows,*,nQ=4,nx=8,nmu=64,order=160,degree=64,radius=640.,spacing=.00625,match_scale=5.,progress=None):
        self.rows=rows;self.grid=RadialGrid(radius,round(radius/spacing))
        low=min(float(HighScaleProfile(r['Qlow'],match_scale).mu(np.inf)) for r in rows)
        self.rg=make_fixed_rg(FixedFlavorWindow(5,4.92,172.5),low,91.1876,.1180024,max(r['Qhigh'] for r in rows))
        self.groups={};self.nodes=[];self.settings=dict(nQ=nQ,nx=nx,nmu=nmu,order=order,degree=degree,radius=radius,spacing=spacing,match_scale=match_scale)
        for i,r in enumerate(rows):self.groups.setdefault((r['Qlow'],r['Qhigh'],r['xflow'],r['xfhigh']),[]).append(i)
        for group,indices in self.groups.items():
            q,w=nodes_weights(group[0],group[1],nQ)
            for Q,Qw in zip(q,w):
                symbol=XFAveragedSymbol(pdf,self.rg,float(Q),group[2],group[3],nx=nx,nmu=nmu,order=order,match_scale=match_scale)
                # Every insertion gets its own endpoint-preserving rational action.
                cache={}
                def component(s,c):
                    key=(np.asarray(s).shape,np.asarray(s).tobytes())
                    if key not in cache:cache[key]=symbol.components(s)
                    return cache[key][...,c]
                fits=[fit_capped_symbol(lambda s,c=c:component(s,c),self.rg.s0,degree) for c in range(6)]
                # dQ^2=2QdQ; invariant d sigma/(dy d^2qT), in pb/GeV^2.
                pref=Qw*2*Q*4*np.pi*(1/137.035999084)**2/(9*38.8**2*Q**2)*HBARC2_PB
                self.nodes.append((float(Q),float(pref),indices,symbol,fits))
            if progress:progress(group,len(self.nodes))
    def evaluate(self,model=GaussianCSInput(),mode='flat_pt'):
        result=np.zeros((len(self.rows),6,3));grid=self.grid
        for Q,pref,indices,symbol,fits in self.nodes:
            inputs=model.native_columns(grid.k,Q)
            # Streaming shared shifted solves keeps memory independent of pole count.
            weight=np.sqrt(grid.weights);rhs=inputs*weight[:,None]
            values=np.array([f.infinity*inputs for f in fits])
            for k,shift in enumerate(fits[0].shifts):
                response=grid.solve(shift,rhs)/weight[:,None]
                values+=np.array([f.residues[k] for f in fits])[:,None,None]*response[None,:,:]
            for i in indices:
                r=self.rows[i];result[i]+=pref*pt_average(grid,values.transpose(1,0,2),r['qtlow'],r['qthigh'],mode)
        return result
