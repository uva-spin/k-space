"""Shared-shift, endpoint-constrained radial CSS2 evaluator prototype.

The fit diagnostics are sampled errors, NOT certified continuum suprema.
The scalar target is provided externally; no replacement of its RG solution occurs.
"""
from __future__ import annotations
from dataclasses import dataclass
import numpy as np
from scipy.linalg import null_space, lstsq, solve_banded
from scipy.sparse import diags

@dataclass
class RationalSymbol:
    shifts: np.ndarray
    residues: np.ndarray
    infinity: float
    zero: float
    tau_star: float
    diagnostics: dict

    def __call__(self, s):
        s=np.asarray(s,dtype=float)
        return self.infinity+(self.shifts/(s[...,None]+self.shifts))@self.residues

    def weighted_tail_bound(self, lower):
        """For target constant at s>=lower; uses actual floating constraint residual."""
        if self.diagnostics.get('tail_constraint_imposed') is False:
            raise ValueError('Constant-target tail bound does not apply to this capped fit')
        A1=abs(np.dot(self.residues,self.shifts))
        A2=np.dot(abs(self.residues),self.shifts**2)
        return float((1/lower+1/self.tau_star)*(A1+A2/(lower+self.shifts.min())))

def fit_symbol(target,s0,degree=48):
    if s0<=0 or degree<4: raise ValueError('s0>0 and degree>=4 required')
    tau=s0
    shifts=s0*np.geomspace(1e-12,10,degree)
    sin=np.unique(np.r_[s0*np.geomspace(1e-18,1e3,2200),np.linspace(.35*s0,1.3*s0,550)])
    r=tau/(tau+sin)
    mi=float(target(np.array([2*s0]))[0]);mz=float(target(np.array([0.]))[0])
    h=(target(sin)-mi)/r
    # φ_j = R_tau_j / R_tau_star. Two linear endpoint constraints:
    # sum c=M(0)-M∞; sum c*tau=0 makes the large-s regular tail O(s^-2).
    phi=(shifts/(sin[:,None]+shifts))/r[:,None]
    C=np.vstack((np.ones(degree),shifts/shifts.max()))
    d=np.array([mz-mi,0.])
    c0=C.T@np.linalg.solve(C@C.T,d)
    Z=null_space(C)
    z,_,rank,sv=lstsq(phi@Z,h-phi@c0,cond=1e-8,lapack_driver='gelsd')
    c=c0+Z@z
    c+=C.T@np.linalg.solve(C@C.T,d-C@c)
    sout=np.unique(np.r_[s0*np.geomspace(1e-22,1e5,7000),np.linspace(.3*s0,1.5*s0,1301)])
    residual=(mi+(shifts/(sout[:,None]+shifts))@c-target(sout))*(1+sout/tau)
    diag={'fit_degree':degree,'fit_samples':len(sin),'check_samples':len(sout),
          'sampled_weighted_error':float(np.max(abs(residual))),
          'sampled_relative_to_symbol_scale':float(np.max(abs(residual))/max(1.,np.max(abs(target(sout))))),
          'finite_grid_error_is_global_certificate':False,'svd_relative_cutoff':1e-8,'rank':int(rank),
          'effective_condition':float(sv[0]/sv[rank-1]),'sum_abs_residues':float(np.sum(abs(c))),
          'zero_constraint_residual':float(sum(c)-(mz-mi)),
          'tail_constraint_residual':float(np.dot(c,shifts))}
    return RationalSymbol(shifts,c,mi,mz,tau,diag)

class RadialGrid:
    """Cell-centered finite-volume -radial Laplacian in the 2πk dk metric.

    Regular no-flux origin; outer homogeneous Dirichlet boundary.
    Weighted similarity produces a real symmetric positive matrix.
    This finite box does not impose a physical TMD cutoff.
    """
    def __init__(self,radius=40.,cells=800):
        if radius<=0 or cells<4:raise ValueError('Invalid grid')
        self.radius,self.cells=radius,cells
        edges=np.linspace(0,radius,cells+1);h=radius/cells
        self.k=(edges[:-1]+edges[1:])/2
        self.weights=np.pi*(edges[1:]**2-edges[:-1]**2)
        conduct=2*np.pi*edges[1:-1]/h
        stiffness=np.zeros(cells);stiffness[:-1]+=conduct;stiffness[1:]+=conduct
        stiffness[-1]+=4*np.pi*radius/h
        self.diagonal=stiffness/self.weights
        self.off=-conduct/np.sqrt(self.weights[:-1]*self.weights[1:])
        self.edges=edges
    def matrix(self):
        """Return the symmetric weighted finite-volume radial Laplacian."""
        return diags((self.off,self.diagonal,self.off),(-1,0,1),format='csc')
    def solve(self,shift,rhs):
        band=np.zeros((3,self.cells));band[1]=self.diagonal+shift
        band[0,1:]=self.off;band[2,:-1]=self.off
        return solve_banded((1,1),band,shift*rhs,check_finite=False)
    def responses(self,shifts,inputs):
        inputs=np.asarray(inputs); weighted=np.sqrt(self.weights)[:,None]*inputs
        return np.stack([self.solve(t,weighted) for t in shifts],axis=0)
    def apply(self,symbol,inputs,responses=None):
        inputs=np.asarray(inputs)
        if responses is None:responses=self.responses(symbol.shifts,inputs)
        return symbol.infinity*inputs+np.tensordot(symbol.residues,responses,axes=(0,0))/np.sqrt(self.weights)[:,None]
    def bins(self,values,edges):
        return np.array([np.sum(self.weights[(self.k>=a)&(self.k<b),None]*values[(self.k>=a)&(self.k<b)],axis=0)
                         for a,b in zip(edges[:-1],edges[1:])])


def fit_capped_symbol(target,scale,degree=64):
    """Endpoint-constrained fit for a capped target with a possible 1/s tail.

    Both endpoints are evaluated explicitly. No sum(c*tau)=0 constraint is
    imposed. Diagnostics are sampled and do not certify a global supremum.
    """
    if scale<=0 or degree<8: raise ValueError('Invalid fit scale or degree')
    shifts=scale*np.geomspace(1e-2,1e2,degree)
    sin=scale*np.geomspace(1e-8,1e6,2400)
    r=scale/(scale+sin)
    mi=float(target(np.inf));mz=float(target(0.))
    h=(target(sin)-mi)/r
    phi=(shifts/(sin[:,None]+shifts))/r[:,None]
    C=np.ones((1,degree));d=np.array([mz-mi]);c0=C.T[:,0]*d[0]/degree
    Z=null_space(C)
    v,_,rank,sv=lstsq(phi@Z,h-phi@c0,cond=1e-10,lapack_driver='gelsd')
    c=c0+Z@v;c+=(d[0]-sum(c))/degree
    check=scale*np.geomspace(1e-10,1e7,7601)
    delta=(mi+(shifts/(check[:,None]+shifts))@c-target(check))*(1+check/scale)
    return RationalSymbol(shifts,c,mi,mz,scale,{
        'fit_degree':degree,'rank':int(rank),'svd_relative_cutoff':1e-10,
        'sampled_weighted_error':float(np.max(abs(delta))),
        'finite_grid_error_is_global_certificate':False,
        'sum_abs_residues':float(np.sum(abs(c))),
        'zero_constraint_residual':float(sum(c)-(mz-mi)),
        'tail_constraint_imposed':False,'actual_rational_tail_coefficient':float(c@shifts)})
