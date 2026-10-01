"""Independent continuum Bessel/bin reference; never used by production training."""
import sys
from pathlib import Path
import numpy as np
from scipy.special import j0
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from native_mixture import MixtureKernel,ResolvedLuminosity,VARIANCES,FIXED_G,pair_indices
from e866_forward import nodes_weights,HighScaleProfile
from css2_noncanonical import FixedFlavorWindow,make_fixed_rg
from dy_measurement import HBARC2_PB

def build_reference(pdf,rows,*,variances=VARIANCES,g=FIXED_G,nQ=4,nx=8,nmu=64,order=160,nb=32,nq=16,measurement='reported_mean'):
    variances=np.asarray(variances);pairs=pair_indices(len(variances))
    low=min(float(HighScaleProfile(r['Qlow']).mu(np.inf)) for r in rows)
    rg=make_fixed_rg(FixedFlavorWindow(5,4.92,172.5),low,91.1876,.1180024,max(r['Qhigh'] for r in rows))
    b=[];bw=[]
    for lo,hi in zip([0,.08,.16,.25,.5,1,2,4,8],[.08,.16,.25,.5,1,2,4,8,16]):
        xx,ww=nodes_weights(lo,hi,nb);b.extend(xx);bw.extend(ww)
    b=np.array(b);bw=np.array(bw)*b/(2*np.pi);out=[];xa=[];xb=[];allQ=[]
    for window in sorted(set((r['Qlow'],r['Qhigh']) for r in rows)):
        qs,qws=nodes_weights(*window,nQ)
        for Q,Qw in zip(qs,qws):
            variance=variances[pairs[:,0]]+variances[pairs[:,1]]+8*g*np.log(Q/5.)
            boundary=np.exp(-b[:,None]**2*variance[None,:]/4)
            pref=Qw*2*Q*4*np.pi*(1/137.035999084)**2/(9*38.8**2*Q**2)*HBARC2_PB
            for group in sorted(set((r['xflow'],r['xfhigh']) for r in rows if (r['Qlow'],r['Qhigh'])==window)):
                xf,xw=nodes_weights(*group,nx);xw/=group[1]-group[0]
                lum=ResolvedLuminosity(pdf,rg,Q,xf,nmu,order);components=lum(b*b)
                ids=[i for i,r in enumerate(rows) if (r['Qlow'],r['Qhigh'],r['xflow'],r['xfhigh'])==window+group]
                for ix,x in enumerate(xf):
                    ans=np.zeros((len(rows),len(pairs),6))
                    for i in ids:
                        r=rows[i];q,qw=nodes_weights(r['qtlow'],r['qthigh'],nq);qw/=r['qthigh']-r['qtlow']
                        if measurement=='reported_mean':qw*=r['Emean']/r['qtmean']*q/np.sqrt(Q*Q+q*q+(38.8*x/2)**2)
                        elif measurement=='uniform_invariant':pass
                        elif measurement=='physical_mean':qw*=np.sqrt(r['Qmean']**2+r['qtmean']**2+(38.8*r['xfmean']/2)**2)/r['qtmean']*q/np.sqrt(Q*Q+q*q+(38.8*x/2)**2)
                        else:raise ValueError(measurement)
                        kernel=qw@j0(q[:,None]*b[None,:])
                        ans[i]=pref*xw[ix]*np.einsum('b,bp,bc->pc',bw*kernel,boundary,components[:,ix,:])
                    out.append(ans);xa.append(lum.xa[ix]);xb.append(lum.xb[ix]);allQ.append(Q)
    out=np.array(out)
    return MixtureKernel(out.sum(axis=-1),np.array(xa),np.array(xb),np.array(allQ),pairs,rows,dict(method='independent Bessel reference',variances=variances.tolist(),g=g,nQ=nQ,nx=nx,nmu=nmu,order=order,nb=nb,nq=nq,measurement=measurement),out)
