import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
import numpy as np
from scipy.integrate import quad
from native_mixture import pair_weights,pair_vjp,pair_indices,measurement_weights,BatchRational
from mixture_dnn import MixtureDNN
from rational_radial import RadialGrid,fit_capped_symbol

def test_pair_convolution_weights_and_adjoint():
    rng=np.random.default_rng(91);a=rng.dirichlet(np.ones(6),4);b=rng.dirichlet(np.ones(6),4);pairs=pair_indices(6)
    assert np.allclose(pair_weights(a,b,pairs).sum(axis=1),1)
    c=rng.normal(size=(4,len(pairs)));da,db=pair_vjp(a,b,pairs,c);d=rng.normal(size=a.shape);h=1e-6
    fd=np.sum(c*(pair_weights(a+h*d,b,pairs)-pair_weights(a-h*d,b,pairs)))/(2*h)
    assert abs(fd-np.sum(d*da))<1e-8

def test_mixture_density_normalization_and_reverse_mode():
    net=MixtureDNN();p=net.initialize(13);x=np.array([.04,.12,.4]);rng=np.random.default_rng(8)
    w,states=net.forward(p,x,True);assert np.all(w>0) and np.allclose(w.sum(axis=1),1)
    for xx in x:
        norm=quad(lambda k:2*np.pi*k*net.density(p,np.array([xx]),np.array([k]))[0,0],0,np.inf)[0]
        assert abs(norm-1)<1e-10
    dw=rng.normal(size=w.shape);gradient=net.reverse(p,states,dw);direction=rng.normal(size=p.shape);h=1e-6
    fd=np.sum(dw*(net.forward(p+h*direction,x)-net.forward(p-h*direction,x)))/(2*h)
    assert abs(fd-gradient@direction)<1e-8

def test_conditional_measurement_constant_field_and_origin():
    grid=RadialGrid(8,1600);row=dict(qtlow=0.,qthigh=.253,qtmean=.17,Emean=10.)
    Q=7.;xf=.22;w=measurement_weights(grid,row,Q,xf)
    exact=10/.17/.253*(np.sqrt(Q**2+(38.8*xf/2)**2+.253**2)-np.sqrt(Q**2+(38.8*xf/2)**2))
    assert np.all(np.isfinite(w)) and np.all(w>=0) and abs(w.sum()-exact)<1e-12

def test_batch_rational_keeps_endpoint_and_original_action():
    def scalar(s):return 2+3/(1+np.asarray(s)/.2)-.4/(1+np.asarray(s)/2)
    fit=fit_capped_symbol(scalar,.1,32);batch=BatchRational(.1,32)
    inf,res=batch.fit(lambda s:np.stack([scalar(s),2*scalar(s)],axis=-1))
    assert np.allclose(inf,[2,4])
    s=np.geomspace(1e-7,1e5,100)
    got=inf+(batch.shifts/(s[:,None]+batch.shifts))@res
    assert np.max(abs(got[:,0]-fit(s)))<1e-8
    assert np.allclose(inf+res.sum(axis=0),[scalar(0),2*scalar(0)])

def test_gaussian_pair_input_matches_direct_two_dimensional_convolution():
    from scipy.special import roots_legendre
    q=.8;va=.2;vb=1.25
    z,w=roots_legendre(100);k=4*(z+1);kw=4*w;phi=np.pi*(z+1);pw=np.pi*w
    integral=np.sum(kw[:,None]*pw[None,:]*k[:,None]*np.exp(-k[:,None]**2/va)/(np.pi*va)*np.exp(-(q*q+k[:,None]**2-2*q*k[:,None]*np.cos(phi))/vb)/(np.pi*vb))
    expected=np.exp(-q*q/(va+vb))/(np.pi*(va+vb))
    assert abs(integral/expected-1)<1e-11

def test_full_neural_pair_objective_gradient_including_normalization():
    from native_mixture import MixtureKernel
    from mixture_dnn import Objective
    rng=np.random.default_rng(19);pairs=pair_indices(6);kernel=MixtureKernel(rng.uniform(.01,.2,(5,3,len(pairs))),np.linspace(.06,.4,5),np.linspace(.2,.05,5),np.full(5,8.),pairs,[{}, {}, {}],{})
    net=MixtureDNN();par=np.r_[net.initialize(3),.02];obj=Objective(net,kernel,[.5,.8,.6],[.03,.04,.05]);loss,grad=obj(par)
    d=rng.normal(size=par.shape);d/=np.linalg.norm(d);h=1e-5
    fd=(obj(par+h*d)[0]-obj(par-h*d)[0])/(2*h)
    assert abs(fd-grad@d)/max(1,abs(fd))<1e-7
