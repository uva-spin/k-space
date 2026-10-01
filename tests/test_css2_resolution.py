"""Scientific regression tests for the corrected CSS2 prescription."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
import numpy as np
import pytest
from scipy.integrate import quad
from css2_coefficients import qcd_coefficients,CF
from css2_matching import match_nlo,convolutions,MatchingBlock,strict_nlo,hard_nlo
from css2_reference import make_rg,DiagnosticProfile,profile_leg
from resummed_operator import Coefficients,CanonicalRG,B0,apply_chebyshev_resolvent

@pytest.fixture(scope='module')
def rg(): return make_rg()


def test_full_evolution_orders_and_demo_overlap():
    c=qcd_coefficients(5)
    assert tuple(map(len,(c.beta,c.cusp,c.noncusp,c.d)))==(4,4,3,3)
    assert c.evolution_accuracy=='N3LL'
    assert c.hard_and_ope_order=='strict NLO'
    assert c.nf==5
    assert c.source_manifest=='manifests/css2_sources.json'
    demo=Coefficients.load(Path(__file__).resolve().parents[1]/'data/coefficients_demo.json')
    for key in ('beta','cusp','noncusp','d'):
        ref=getattr(demo,key)
        assert np.allclose(getattr(c,key)[:len(ref)],ref,rtol=2e-14,atol=1e-12)
    assert not c.physical_n3llprime_payload


def test_physical_coefficient_file_matches_formula_implementation():
    path=Path(__file__).resolve().parents[1]/'data/coefficients_css2_nf5_n3ll.json'
    loaded=Coefficients.load(path);generated=qcd_coefficients(5)
    for key in ('beta','cusp','noncusp','d'):
        assert np.allclose(getattr(loaded,key),getattr(generated,key),rtol=0,atol=2e-12)
    assert loaded.evolution_accuracy=='N3LL'
    assert loaded.hard_and_ope_order=='strict NLO'
    assert loaded.nf==5


def test_n3ll_metadata_requires_complete_evolution_vectors(tmp_path):
    import json
    path=tmp_path/'bad.json'
    path.write_text(json.dumps({'beta':[1.], 'cusp':[1.], 'noncusp':[1.], 'd':[0.],
                                'description':'bad', 'evolution_accuracy':'N3LL'}))
    with pytest.raises(ValueError,match='4/4/3/3'):
        Coefficients.load(path)


def test_chebyshev_resolvent_matches_dense_function_of_matrix(rg):
    from scipy.sparse import diags
    n=15
    matrix=diags([-np.ones(n-1),2*np.ones(n),-np.ones(n-1)],[-1,0,1],format='csc')
    values,vectors=np.linalg.eigh(matrix.toarray())
    initial=(1.0+np.arange(n)**2)**-2
    coefficients=rg.chebyshev_coefficients(96,rg.s0,20.,400.)
    actual=apply_chebyshev_resolvent(matrix,initial,rg.s0,coefficients)
    scalar=np.polynomial.chebyshev.chebval(2*rg.s0/(rg.s0+values)-1,coefficients)
    expected=vectors@(scalar*(vectors.T@initial))
    assert np.linalg.norm(actual-expected)<2e-11
    assert np.polynomial.chebyshev.chebval(1.,coefficients)==pytest.approx(0.,abs=2e-13)
    assert np.polynomial.chebyshev.chebval(-1.,coefficients)==pytest.approx(
        rg.v(np.array([2*rg.s0]),20.,400.)[0],abs=2e-13)


def test_method_paper_accuracy_manifest_keeps_branches_and_orders_separate():
    import json
    root=Path(__file__).resolve().parents[1]
    manifest=json.loads((root/'manifests/method_paper_accuracy.json').read_text())
    assert manifest['central_accuracy_statement']=='N3LL evolution with strict NLO hard and OPE boundaries'
    assert manifest['canonical_branch']['coefficient_orders']=={
        'beta':4,'cusp':4,'noncusp':3,'canonical_CS_boundary':3}
    assert manifest['capped_noncanonical_branch']['inherits_canonical_N3LL_tower_theorem'] is False
    assert 'physical N3LL-prime spectrum' in manifest['not_claimed']

@pytest.mark.parametrize('s',[.002,.017,.04])
def test_two_scale_rg_cs_and_path_independence(rg,s):
    mu,zeta,h=17.,31.**2,1e-5
    a=rg.values(np.log(mu/rg.mu0))[0]
    dmu=(rg.log_v_tc(rg.tc(s),mu*np.exp(h),zeta)-rg.log_v_tc(rg.tc(s),mu*np.exp(-h),zeta))/(2*h)
    expected=rg.c.series(a,rg.c.noncusp)-rg.c.series(a,rg.c.cusp)*np.log(zeta/mu**2)
    assert dmu==pytest.approx(expected,abs=2e-9)
    dz=(rg.log_v_tc(rg.tc(s),mu,zeta*np.exp(2*h))-rg.log_v_tc(rg.tc(s),mu,zeta*np.exp(-2*h)))/(2*h)
    assert dz==pytest.approx(-rg.D(s,mu),abs=2e-9)
    dD=(rg.D(s,mu*np.exp(h))-rg.D(s,mu*np.exp(-h)))/(2*h)
    assert dD==pytest.approx(2*rg.c.series(a,rg.c.cusp),abs=2e-9)
    e0=(6.,49.);e1=(12.,100.);e2=(23.,700.)
    assert rg.log_u(s,e2,e0)==pytest.approx(rg.log_u(s,e2,e1)+rg.log_u(s,e1,e0),abs=2e-11)
    # Independently integrate each side of a rectangular scale path.
    mu_leg=quad(lambda t: rg.c.series(rg.values(t)[0],rg.c.noncusp)
                -rg.c.series(rg.values(t)[0],rg.c.cusp)*np.log(e0[1]/(rg.mu0*np.exp(t))**2),
                np.log(e0[0]/rg.mu0),np.log(e2[0]/rg.mu0),epsabs=1e-12)[0]
    rect=mu_leg-rg.D(s,e2[0])*np.log(np.sqrt(e2[1]/e0[1]))
    assert rect==pytest.approx(rg.log_u(s,e2,e0),abs=3e-11)


def test_fixed_coupling_sudakov_normalization():
    a,L=.013,1.7
    single=quad(lambda t:a*(6*CF-8*CF*(L-t)),0,L)[0]
    correct_S=-2*single
    assert correct_S==pytest.approx(a*(8*CF*L*L-12*CF*L))
    historical_S=a*(4*CF*L*L-6*CF*L)
    assert historical_S==pytest.approx(correct_S/2)

@pytest.mark.parametrize('x',[.03,.2,.7])
def test_general_matching_virtuality_and_rapidity_equations(x):
    q=lambda y:(1-y)**3
    g=lambda y:2*(1-y)**5
    L,ell,h=.7,-.4,1e-5
    # a=1 extracts the formal one-loop coefficient (not a physical coupling).
    def coef(L,e):
        m=match_nlo(x,q,g,1.,L,e,order=240)
        return m.qq+m.qg
    _,pqq,pqg,_,_=convolutions(x,q,g,240)
    dmu=(coef(L+2*h,ell-2*h)-coef(L-2*h,ell+2*h))/(2*h)
    assert dmu+4*(pqq+pqg)==pytest.approx((6*CF-4*CF*ell)*q(x),abs=2e-8)
    dz=(coef(L,ell+2*h)-coef(L,ell-2*h))/(2*h)
    assert dz==pytest.approx(-4*CF*L*q(x),abs=2e-8)


def test_plus_prescription_against_independent_adaptive_integral():
    x=.2;q=lambda t:(1-t)**3;g=lambda t:(1-t)**5
    fx,p,_,_,_=convolutions(x,q,g,200)
    expected=CF*(quad(lambda z:(1+z*z)/(1-z)*(q(x/z)/z-fx),x,1,epsabs=1e-12)[0]
                 -fx*quad(lambda z:(1+z*z)/(1-z),0,x,epsabs=1e-12)[0])
    assert p==pytest.approx(expected,abs=1e-12)


def test_strict_insertion_projection():
    A=MatchingBlock(2.,.3,-.2);B=MatchingBlock(3.,-.1,.4);dh=.17;ev=.83
    result=strict_nlo(A,B,dh,ev)
    product=lambda e: ev*(1+e*dh)*(A.born+e*(A.qq+A.qg))*(B.born+e*(B.qq+B.qg))
    # A cubic polynomial has a five-point derivative exact up to rounding.
    h=.1;deriv=(product(-2*h)-8*product(-h)+8*product(h)-product(2*h))/(12*h)
    assert result['strict_total']==pytest.approx(product(0)+deriv,abs=1e-13)
    assert abs(result['strict_total']-product(1))>.01


def test_hard_virtuality_cancels_two_legs_at_one_loop():
    a,Q,mu,h=.01,20.,17.,1e-5
    der=(hard_nlo(a,Q,mu*np.exp(h))-hard_nlo(a,Q,mu*np.exp(-h)))/(2*h)
    assert der+2*a*(6*CF-4*CF*np.log(Q**2/mu**2))==pytest.approx(0,abs=1e-10)


def test_capped_endpoints_do_not_inherit_uncapped_zero(rg):
    profile=DiagnosticProfile()
    assert profile.mu(0)==20.
    assert profile.mu(np.inf)==pytest.approx(5.)
    assert profile_leg(rg,profile,0)==1.
    assert np.isfinite(profile_leg(rg,profile,np.inf))
    assert profile_leg(rg,profile,np.inf)>0
    assert rg.v(0.,20.,400.)==0.
    assert np.array_equal(rg.v(np.zeros(3),20.,400.),np.zeros(3))

@pytest.mark.parametrize('s',[-1.,np.nan])
def test_profile_rejects_invalid_s(s):
    with pytest.raises(ValueError): DiagnosticProfile().mu(s)


def test_cusp_against_independent_published_values():
    import json
    p=Path(__file__).resolve().parents[1]/'data/css2_cusp_independent.json'
    expected=json.loads(p.read_text())['expected']
    for nf,values in expected.items():
        assert np.allclose(qcd_coefficients(int(nf)).cusp,values,rtol=1e-12,atol=2e-9)


def test_capped_fit_keeps_nonzero_origin_and_nonzero_tail():
    from rational_radial import fit_capped_symbol
    symbol=fit_capped_symbol(lambda s:2+3/(1+np.asarray(s)),1.,48)
    assert symbol(0)==pytest.approx(5.,abs=1e-11)
    assert symbol(np.inf)==2.
    s=np.geomspace(1e-4,1e4,301)
    assert np.max(abs((symbol(s)-(2+3/(1+s)))*(1+s)))<2e-6
    with pytest.raises(ValueError):symbol.weighted_tail_bound(10.)
