"""Small fully connected native-kT mixture-density network and exact reverse mode.

No external ML framework, finite-difference training or transformed trainable
input is used. Tanh hidden layers and softmax coefficients define an explicit
finite positive radial model, not an unrestricted renormalized-TMD prior.
"""
import numpy as np
from scipy.special import softmax
from scipy.optimize import minimize
from native_mixture import VARIANCES

class MixtureDNN:
    def __init__(self,widths=(1,12,12,12,6)):
        if widths[0]!=1 or widths[-1]!=len(VARIANCES):raise ValueError('Unsupported input/basis')
        self.widths=tuple(widths);self.slices=[];offset=0
        for a,b in zip(widths[:-1],widths[1:]):
            self.slices.append((slice(offset,offset+a*b),(a,b),slice(offset+a*b,offset+a*b+b)))
            offset+=a*b+b
        self.size=offset
    def initialize(self,seed):
        rng=np.random.default_rng(seed);p=np.zeros(self.size)
        for j,(ws,shape,bs) in enumerate(self.slices):
            p[ws]=rng.normal(size=np.prod(shape))*np.sqrt(2/sum(shape))*(.15 if j==len(self.slices)-1 else 1)
        p[self.slices[-1][2]]=np.log([.04,.12,.45,.28,.09,.02])
        return p
    def forward(self,p,x,cache=False):
        x=np.asarray(x)
        if np.any(x<=0) or np.any(x>=1):raise ValueError('x outside (0,1)')
        a=(np.log(x/.1)/np.log(4))[:,None];states=[a]
        for j,(ws,shape,bs) in enumerate(self.slices):
            z=a@p[ws].reshape(shape)+p[bs]
            a=softmax(z,axis=1) if j==len(self.slices)-1 else np.tanh(z)
            states.append(a)
        return (a,states) if cache else a
    def reverse(self,p,states,dweights):
        w=states[-1];dz=w*(dweights-np.sum(w*dweights,axis=1,keepdims=True));grad=np.zeros_like(p)
        for j in range(len(self.slices)-1,-1,-1):
            ws,shape,bs=self.slices[j];grad[ws]=(states[j].T@dz).ravel();grad[bs]=dz.sum(axis=0)
            if j:dz=(dz@p[ws].reshape(shape).T)*(1-states[j]**2)
        return grad
    def density(self,p,x,k):
        w=self.forward(p,np.asarray(x));v=VARIANCES
        return w@(np.exp(-np.asarray(k)[None,:]**2/v[:,None])/(np.pi*v[:,None]))

def truth_weights(x):
    """Smooth non-Gaussian closure generator; training receives predictions only."""
    z=np.log(np.asarray(x)/.12)
    logits=np.stack([-2.2+.25*np.sin(z),-.2-.55*z,.9+.35*np.cos(1.3*z),.25+.65*z,-1.3+.25*np.sin(2*z),-3.5+.2*z],axis=-1)
    return softmax(logits,axis=-1)

class Objective:
    def __init__(self,net,kernel,data,sigma,mask=None,*,normalization=True,regularization=1e-4,fixed_leg_b=None,smoothness=.1):
        self.net,self.kernel=net,kernel;self.data=np.asarray(data);self.sigma=np.asarray(sigma)
        self.mask=np.ones(len(data),bool) if mask is None else np.asarray(mask,bool)
        if np.any(self.sigma<=0):raise ValueError('Require positive uncertainties')
        self.normalization,self.regularization=normalization,regularization
        self.x=np.r_[kernel.xa,kernel.xb];self.fixed_leg_b=fixed_leg_b
        self.smoothness=smoothness;self.smooth_x=np.geomspace(.04,.70,257)
        self.smooth_h=np.log(self.smooth_x[1]/self.smooth_x[0])
    def predict(self,p):
        w=self.net.forward(p[:self.net.size],self.x);n=len(self.kernel.xa)
        return self.kernel.evaluate(w[:n],w[n:] if self.fixed_leg_b is None else self.fixed_leg_b)
    def __call__(self,p):
        theta=p[:self.net.size];N=np.exp(p[-1]) if self.normalization else 1.
        weights,states=self.net.forward(theta,self.x,cache=True);n=len(self.kernel.xa)
        wa,wb=weights[:n],weights[n:] if self.fixed_leg_b is None else self.fixed_leg_b
        pred=self.kernel.evaluate(wa,wb);res=(N*pred-self.data)/self.sigma
        res=np.where(self.mask,res,0);loss=.5*np.dot(res,res)+.5*self.regularization*np.dot(theta,theta)
        da,db=self.kernel.vjp(wa,wb,N*res/self.sigma)
        if self.fixed_leg_b is not None:db*=0
        grad=self.net.reverse(theta,states,np.r_[da,db])+self.regularization*theta
        if self.smoothness:
            sw,scache=self.net.forward(theta,self.smooth_x,cache=True)
            curvature=(sw[2:]-2*sw[1:-1]+sw[:-2])/self.smooth_h**2
            loss+=.5*self.smoothness*np.sum(curvature**2)/len(curvature)
            dcur=self.smoothness*curvature/(len(curvature)*self.smooth_h**2)
            ds=np.zeros_like(sw);ds[2:]+=dcur;ds[1:-1]-=2*dcur;ds[:-2]+=dcur
            grad+=self.net.reverse(theta,scache,ds)
        if self.normalization:
            rn=(N-1)/.065;loss+=.5*rn*rn
            gn=np.dot(res,N*pred/self.sigma)+rn*N/.065
            grad=np.r_[grad,gn]
        return float(loss),grad

def train(objective,seed,*,adam_steps=800,maxiter=1500,progress=None):
    p=objective.net.initialize(seed)
    if objective.normalization:p=np.r_[p,0.]
    m=np.zeros_like(p);v=m.copy();history=[];best=(np.inf,p.copy())
    for step in range(1,adam_steps+1):
        loss,grad=objective(p)
        if loss<best[0]:best=(loss,p.copy())
        # Clip only the optimizer's vector norm, never physics values.
        grad*=min(1.,100/max(np.linalg.norm(grad),1e-100))
        m=.9*m+.1*grad;v=.999*v+.001*grad*grad
        lr=.012*(.2+.8*(1-step/adam_steps))
        p-=lr*(m/(1-.9**step))/(np.sqrt(v/(1-.999**step))+1e-8)
        if objective.normalization:p[-1]=np.clip(p[-1],np.log(.5),np.log(1.5))
        if step%200==0:
            history.append(dict(stage='Adam',step=step,loss=loss))
            if progress:progress('Adam',step,loss)
    bounds=[(None,None)]*len(p)
    if objective.normalization:bounds[-1]=(np.log(.5),np.log(1.5))
    opt=minimize(objective,best[1],jac=True,method='L-BFGS-B',bounds=bounds,options=dict(maxiter=maxiter,maxfun=4*maxiter,maxls=40,ftol=1e-12,gtol=2e-6,maxcor=30))
    loss,grad=objective(opt.x)
    return opt.x,dict(success=bool(opt.success),message=str(opt.message),iterations=opt.nit,evaluations=opt.nfev,loss=loss,gradient_inf=float(np.max(abs(grad))),history=history)

def function_sensitivity(net,p,kernel,sigma,knots=8):
    """Whitened Jacobian for local logit functions, not redundant neural weights."""
    x=np.r_[kernel.xa,kernel.xb];lx=np.log(x);grid=np.linspace(lx.min(),lx.max(),knots)
    H=np.array([np.interp(lx,grid,np.eye(knots)[j]) for j in range(knots)]).T
    w=net.forward(p[:net.size],x);n=len(kernel.xa);J=[]
    for r in range(len(kernel.rows)):
        e=np.zeros(len(kernel.rows));e[r]=1
        da,db=kernel.vjp(w[:n],w[n:],e);dw=np.r_[da,db]
        logits=w*(dw-np.sum(w*dw,axis=1,keepdims=True))
        J.append((H.T@logits[:,:-1]).ravel()/sigma[r])
    u,s,vh=np.linalg.svd(J,full_matrices=False)
    return dict(coordinates='8 log-x hat functions times 5 relative mixture logits',jacobian_shape=list(np.asarray(J).shape),singular_values=s.tolist(),relative_rank_1e3=int(np.sum(s>s[0]*1e-3)),relative_rank_1e6=int(np.sum(s>s[0]*1e-6)),scope='Local sensitivity in declared coordinates; not a unique functional reconstruction or neural-parameter covariance')
