"""Native-kT Gaussian-basis pair operator with x-resolved measurement weights.

A mixture-density network learns boundary weights, not perturbative functions.
Gaussian pair convolution is performed analytically in kT; no production Fourier
transform is used. This finite positive basis is an explicit model restriction.
"""
from pathlib import Path
import json
import numpy as np
from scipy.linalg import null_space,lstsq
from scipy.interpolate import CubicSpline
from e866_forward import (HighScaleProfile,canonical_blocks,parton_x,CHARGES,
                          nodes_weights,GaussianCSInput)
from css2_matching import hard_nlo
from css2_noncanonical import FixedFlavorWindow,make_fixed_rg
from css2_reference import profile_leg
from rational_radial import RadialGrid
from dy_measurement import HBARC2_PB

VARIANCES=np.array([.08,.20,.50,1.25,3.125,7.8125])
FIXED_G=.15790301183804753

def pair_indices(size):return np.array([(i,j) for i in range(size) for j in range(i,size)])

def pair_weights(a,b,pairs):
    i,j=pairs.T
    return a[...,i]*b[...,j]+(i!=j)*a[...,j]*b[...,i]

def pair_vjp(a,b,pairs,gradient):
    da=np.zeros_like(a);db=np.zeros_like(b)
    for p,(i,j) in enumerate(pairs):
        da[:,i]+=gradient[:,p]*b[:,j];db[:,j]+=gradient[:,p]*a[:,i]
        if i!=j:da[:,j]+=gradient[:,p]*b[:,i];db[:,i]+=gradient[:,p]*a[:,j]
    return da,db

class BatchRational:
    """Same endpoint projector and SVD cutoff as fit_capped_symbol, batch RHS."""
    def __init__(self,scale,degree=64):
        self.shifts=scale*np.geomspace(1e-2,1e2,degree);self.s=scale*np.geomspace(1e-8,1e6,2400)
        self.r=scale/(scale+self.s);self.phi=(self.shifts/(self.s[:,None]+self.shifts))/self.r[:,None]
        self.Z=null_space(np.ones((1,degree)));self.A=self.phi@self.Z
    def fit(self,target):
        inf=target(np.inf);zero=target(0.);shape=inf.shape
        values=target(self.s).reshape(len(self.s),-1);inff=inf.ravel();d=(zero-inf).ravel()
        c0=np.broadcast_to(d[None,:]/len(self.shifts),(len(self.shifts),len(d))).copy()
        h=(values-inff)/self.r[:,None]
        z,_,_,_=lstsq(self.A,h-self.phi@c0,cond=1e-10,lapack_driver='gelsd')
        c=c0+self.Z@z;c+=(d-c.sum(axis=0))/len(self.shifts)
        return inf,c.reshape((len(self.shifts),)+shape)

class ResolvedLuminosity:
    def __init__(self,pdf,rg,Q,xf,nmu=64,order=160):
        self.profile=HighScaleProfile(float(Q));self.rg=rg
        self.xa,self.xb=parton_x(Q,xf)
        pad=16*np.finfo(float).eps
        self.t=np.linspace(np.log(float(self.profile.mu(np.inf))/rg.mu0)-pad,np.log(Q/rg.mu0)+pad,nmu)
        h=hard_nlo(float(rg.values(np.log(Q/rg.mu0))[0]),Q,Q)
        def direct(mu):
            a=canonical_blocks(pdf,rg,self.xa,mu,order);b=canonical_blocks(pdf,rg,self.xb,mu,order)[...,::-1,:]
            base=a[...,0]*b[...,0]
            c=np.stack([base,h*base,a[...,1]*b[...,0],a[...,2]*b[...,0],a[...,0]*b[...,1],a[...,0]*b[...,2]],axis=-1)
            return np.einsum('f,xfc->xc',CHARGES,c)
        self.cache=CubicSpline(self.t,np.array([direct(rg.mu0*np.exp(t)) for t in self.t]),axis=0,extrapolate=False)
    def __call__(self,s):
        mu=self.profile.mu(s);v=self.cache(np.log(mu/self.rg.mu0))
        return v*np.asarray(profile_leg(self.rg,self.profile,s)**2)[...,None,None]

def measurement_weights(grid,row,Q,xf,mode='reported_mean'):
    lo=np.maximum(row['qtlow'],grid.edges[:-1]);hi=np.minimum(row['qthigh'],grid.edges[1:]);hi=np.maximum(lo,hi)
    if mode=='uniform_invariant':return (hi-lo)/(row['qthigh']-row['qtlow'])
    C=Q*Q+(38.8*xf/2)**2
    # Exact integral of q/E within each cell; stable at q=0.
    integral=(hi*hi-lo*lo)/(np.sqrt(C+hi*hi)+np.sqrt(C+lo*lo))
    if mode=='reported_mean':E=row['Emean']
    elif mode=='physical_mean':E=np.sqrt(row['Qmean']**2+row['qtmean']**2+(38.8*row['xfmean']/2)**2)
    else:raise ValueError('Unknown registered measurement functional')
    return E/row['qtmean']/(row['qthigh']-row['qtlow'])*integral

class MixtureKernel:
    def __init__(self,kernel,xa,xb,Q,pairs,rows,settings,components=None):
        self.kernel,self.xa,self.xb,self.Q,self.pairs=kernel,xa,xb,Q,pairs
        self.rows,self.settings,self.components=rows,settings,components
    def evaluate(self,wa,wb):return np.einsum('nrp,np->r',self.kernel,pair_weights(wa,wb,self.pairs))
    def vjp(self,wa,wb,residual_gradient):
        gp=np.einsum('nrp,r->np',self.kernel,residual_gradient)
        return pair_vjp(wa,wb,self.pairs,gp)
    def save(self,path):
        np.savez_compressed(path,kernel=self.kernel,xa=self.xa,xb=self.xb,Q=self.Q,pairs=self.pairs,components=self.components,rows=json.dumps(self.rows),settings=json.dumps(self.settings))
    @classmethod
    def load(cls,path):
        with np.load(path,allow_pickle=False) as z:return cls(z['kernel'],z['xa'],z['xb'],z['Q'],z['pairs'],json.loads(str(z['rows'])),json.loads(str(z['settings'])),z['components'])

def build_native(pdf,rows,*,variances=VARIANCES,g=FIXED_G,nQ=4,nx=8,nmu=64,order=160,degree=64,radius=640.,spacing=.00625,progress=None):
    variances=np.asarray(variances);pairs=pair_indices(len(variances))
    if np.any(variances<=0) or g<0:raise ValueError('Invalid boundary/CS parameters')
    low=min(float(HighScaleProfile(r['Qlow']).mu(np.inf)) for r in rows)
    rg=make_fixed_rg(FixedFlavorWindow(5,4.92,172.5),low,91.1876,.1180024,max(r['Qhigh'] for r in rows))
    grid=RadialGrid(radius,round(radius/spacing));rational=BatchRational(rg.s0,degree)
    cut=min(grid.cells,np.searchsorted(grid.k,max(r['qthigh'] for r in rows))+2)
    comp=[];xa=[];xb=[];allQ=[];weight=np.sqrt(grid.weights)
    windows=sorted(set((r['Qlow'],r['Qhigh']) for r in rows))
    for window in windows:
        Qs,Qws=nodes_weights(*window,nQ)
        for Q,Qw in zip(Qs,Qws):
            variance=variances[pairs[:,0]]+variances[pairs[:,1]]+8*g*np.log(Q/5.)
            inputs=np.exp(-grid.k[:,None]**2/variance)/(np.pi*variance)
            rhs=weight[:,None]*inputs;fields=[inputs[:cut]]
            for shift in rational.shifts:fields.append(grid.solve(shift,rhs)[:cut]/weight[:cut,None])
            fields=np.array(fields)
            pref=Qw*2*Q*4*np.pi*(1/137.035999084)**2/(9*38.8**2*Q**2)*HBARC2_PB
            groups=sorted(set((r['xflow'],r['xfhigh']) for r in rows if (r['Qlow'],r['Qhigh'])==window))
            for group in groups:
                xf,xw=nodes_weights(*group,nx);xw/=group[1]-group[0]
                lum=ResolvedLuminosity(pdf,rg,Q,xf,nmu,order);inf,residues=rational.fit(lum)
                coeff=np.concatenate([inf[None,:,:],residues],axis=0)
                ids=[i for i,r in enumerate(rows) if (r['Qlow'],r['Qhigh'],r['xflow'],r['xfhigh'])==window+group]
                for x,xx in enumerate(xf):
                    out=np.zeros((len(rows),len(pairs),6))
                    for i in ids:
                        mw=measurement_weights(grid,rows[i],Q,xx)[:cut]
                        moments=np.einsum('k,skp->sp',mw,fields)
                        out[i]=pref*xw[x]*np.einsum('sp,sc->pc',moments,coeff[:,x,:])
                    comp.append(out);xa.append(lum.xa[x]);xb.append(lum.xb[x]);allQ.append(Q)
            if progress:progress(float(Q),len(comp))
    comp=np.array(comp)
    settings=dict(method='native elliptic solves and direct kT Gaussian convolution',variances=variances.tolist(),g=g,Q0=5.,nQ=nQ,nx=nx,nmu=nmu,order=order,degree=degree,radius=radius,spacing=spacing,measurement='reported_mean')
    return MixtureKernel(comp.sum(axis=-1),np.array(xa),np.array(xb),np.array(allQ),pairs,rows,settings,comp)
