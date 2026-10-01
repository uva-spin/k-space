from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
import numpy as np
import pytest
import sympy as sp
from css2_noncanonical import cs_boundary,NoncanonicalProfile,BoundaryEvolution,FixedFlavorWindow
from css2_coefficients import qcd_coefficients
from css2_reference import make_rg
from dy_measurement import InvariantPoint,rapidity_from_xf,HBARC2_PB
from resummed_operator import B0


def test_general_CS_polynomial_satisfies_RG_through_three_loops():
    a,L=sp.symbols('a L');b0,b1,g0,g1,g2,d2,d3=sp.symbols('b0 b1 g0 g1 g2 d2 d3')
    D=a*g0*L+a**2*(b0*g0*L**2/2+g1*L+d2)+a**3*(b0**2*g0*L**3/3+(b1*g0/2+b0*g1)*L**2+(g2+2*b0*d2)*L+d3)
    residual=2*sp.diff(D,L)-2*(b0*a**2+b1*a**3)*sp.diff(D,a)-2*(g0*a+g1*a**2+g2*a**3)
    assert sp.series(residual,a,0,4).removeO().expand()==0
    # Independent integration of the coefficient recurrence verifies the numeric API.
    c=qcd_coefficients(4)
    substitutions={b0:c.beta[0],b1:c.beta[1],g0:c.cusp[0],g1:c.cusp[1],g2:c.cusp[2],d2:c.d[1],d3:c.d[2]}
    for av,lv in ((.015,-.7),(.026,1.6)):
        assert cs_boundary(av,lv,c)==pytest.approx(float(D.subs(substitutions|{a:av,L:lv})),abs=2e-15)


def test_general_boundary_canonical_limit():
    c=qcd_coefficients(5);a=np.array([.01,.02,.03])
    assert np.allclose(cs_boundary(a,0,c),c.series(a,c.d),atol=1e-16)


def test_noncanonical_fixed_anchor_two_scale_paths():
    rg=make_rg();ev=BoundaryEvolution(rg,.4,7.)
    mu,zeta,h=13.,400.,1e-5
    a=rg.values(np.log(mu/5))[0]
    assert (ev.D(mu*np.exp(h))-ev.D(mu*np.exp(-h)))/(2*h)==pytest.approx(2*rg.c.series(a,rg.c.cusp),abs=1e-9)
    initial=(7.,49.)
    dmu=(ev.log_u((mu*np.exp(h),zeta),initial)-ev.log_u((mu*np.exp(-h),zeta),initial))/(2*h)
    assert dmu==pytest.approx(rg.c.series(a,rg.c.noncusp)-rg.c.series(a,rg.c.cusp)*np.log(zeta/mu**2),abs=1e-9)
    dz=(ev.log_u((mu,zeta*np.exp(2*h)),initial)-ev.log_u((mu,zeta*np.exp(-2*h)),initial))/(2*h)
    assert dz==pytest.approx(-ev.D(mu),abs=1e-9)
    middle=(9.,81.);target=(20.,400.)
    assert ev.log_u(target,initial)==pytest.approx(ev.log_u(target,middle)+ev.log_u(middle,initial),abs=1e-12)


def test_power16_coordinate_has_real_noncanonical_logs():
    p=NoncanonicalProfile(Q=20.,floor=5.)
    assert p.b_eval(0.)==pytest.approx(B0/20)
    assert p.Lb(0.)==pytest.approx(0.,abs=1e-15)
    assert p.b_eval(np.inf)==pytest.approx(1.5)
    assert p.Lb(np.inf)>3
    assert np.isfinite(p.b_eval(np.array([0,.1,np.inf]))).all()

@pytest.mark.parametrize('nf,lo,hi',[(5,1.3,20),(5,1.65,20),(4,1.65,7.5),(4,1.65,4.92)])
def test_unsupported_thresholds_fail_closed(nf,lo,hi):
    window=FixedFlavorWindow(nf,1.65 if nf==4 else 4.92,4.92 if nf==4 else 172.5)
    with pytest.raises(ValueError):window.validate(lo,hi)


def test_single_flavor_windows_remain_available():
    FixedFlavorWindow(4,1.65,4.92).validate(1.65,4.5)
    FixedFlavorWindow(5,4.92,172.5).validate(5.,20.)


def test_invariant_observable_mass_integral_and_units():
    p=InvariantPoint(7.,8.,.1,-.1,38.8,'analytic-fixture')
    alpha=1/137.;constant=2.3
    expected=4*np.pi*alpha**2/(9*38.8**2)*HBARC2_PB*constant*2*np.log(8/7)
    got=p.evaluate(lambda Q,y,qt:constant,alpha_em=alpha,kinematics='exact_pair_geometry')
    assert got['invariant_pb_GeV_minus2_per_nucleon']==pytest.approx(expected,rel=1e-13)
    assert not got['mass_averaged'] and not got['qT_binned']


def test_fixed_xf_geometry_and_leading_power_limit():
    y=rapidity_from_xf(7.5,.1,-.1,38.8,'exact_pair_geometry')
    assert 2*np.hypot(7.5,.1)*np.sinh(y)/38.8==pytest.approx(-.1)
    y0=rapidity_from_xf(7.5,0,-.1,38.8,'leading_power')
    assert y0==pytest.approx(-.255865709,abs=1e-9)
    assert y0==rapidity_from_xf(7.5,0,-.1,38.8,'exact_pair_geometry')


def test_actual_softcap_endpoint_sets_the_RG_interval():
    from css2_noncanonical import make_fixed_rg
    p=NoncanonicalProfile(Q=7.5,floor=5.)
    actual=float(p.mu(np.inf))
    assert actual<5. and 5.-actual<1e-8
    rg=make_fixed_rg(FixedFlavorWindow(5,4.92,172.5),actual,91.1876,.1180024,7.5)
    assert np.isfinite(rg.values(np.log(p.mu(np.inf)/rg.mu0))).all()


@pytest.mark.parametrize('anchor_mu,anchor_alpha',[(1.3,.3),(4.92,.22),(4.5,float('nan'))])
def test_coupling_anchor_respects_window_and_finite_input(anchor_mu,anchor_alpha):
    from css2_noncanonical import make_fixed_rg
    with pytest.raises(ValueError):
        make_fixed_rg(FixedFlavorWindow(4,1.65,4.92),1.65,anchor_mu,anchor_alpha,4.5)


def test_symbol_rejects_inconsistent_flavor_before_accessing_PDF():
    from css2_noncanonical import NoncanonicalPointSymbol
    # A declared nf=4 window cannot use the nf=5 coefficient vector.
    with pytest.raises(ValueError,match='RG coefficients'):
        NoncanonicalPointSymbol(None,make_rg(),NoncanonicalProfile(Q=4.5,floor=1.65),FixedFlavorWindow(4,1.65,4.92))


def test_symbol_rejects_uncovered_softcap_before_accessing_PDF():
    from css2_noncanonical import NoncanonicalPointSymbol
    with pytest.raises(ValueError,match='actual composed profile endpoint'):
        NoncanonicalPointSymbol(None,make_rg(),NoncanonicalProfile(Q=7.5,floor=5.),FixedFlavorWindow(5,4.92,172.5))
