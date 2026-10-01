"""Explicit invariant DY observable at a published fixed-xF qT point.

W callback: W(Q,y,qT) in GeV^-2, using unweighted PDFs, charge sum, strict
hard/OPE and the chosen resummed/NP prescription. The target mixture belongs
in W. This adapter supplies no guessed nuclear model, qT bins or fit rescaling.
"""
from dataclasses import dataclass
import numpy as np
from scipy.integrate import quad
HBARC2_PB=3.893793656e8


def rapidity_from_xf(Q,qt,xf,sqrts,kinematics):
    if not np.isfinite([Q,qt,xf,sqrts]).all() or Q<=0 or qt<0 or sqrts<=Q:raise ValueError('Invalid kinematics')
    if kinematics=='exact_pair_geometry':mt=np.hypot(Q,qt)
    elif kinematics=='leading_power':mt=Q
    else:raise ValueError('Specify exact_pair_geometry or leading_power')
    return np.arcsinh(xf*sqrts/(2*mt))


@dataclass(frozen=True)
class InvariantPoint:
    Qlow: float
    Qhigh: float
    qt: float
    xf: float
    sqrts: float
    source_id: str
    def evaluate(self,W,*,alpha_em,kinematics):
        if not 0<self.Qlow<self.Qhigh<self.sqrts or not 0<alpha_em<1:
            raise ValueError('Invalid mass interval or alpha_em')
        def integrand(Q):
            y=rapidity_from_xf(Q,self.qt,self.xf,self.sqrts,kinematics)
            value=W(Q,y,self.qt)
            if not np.isfinite(value):raise ValueError('Nonfinite structure function')
            # dQ^2=2Q dQ; invariant = d sigma/(dy d^2qT) at each Q.
            return 2*Q*(4*np.pi*alpha_em**2/(9*self.sqrts**2*Q**2))*value*HBARC2_PB
        value,error=quad(integrand,self.Qlow,self.Qhigh,epsabs=1e-10,epsrel=2e-8,limit=200)
        return {'invariant_pb_GeV_minus2_per_nucleon':value,'quadrature_error_estimate':error,
                'source_id':self.source_id,'kinematics':kinematics,'mass_averaged':False,'qT_binned':False}
