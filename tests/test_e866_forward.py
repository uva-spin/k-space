from pathlib import Path
import json,sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
import numpy as np
import pytest
from e866_forward import HighScaleProfile,GaussianCSInput,parton_x,pt_average,canonical_blocks,FLAVORS
from css2_noncanonical import FixedFlavorWindow
from css2_reference import make_rg
from css2_matching import match_nlo
from rational_radial import RadialGrid
from resummed_operator import B0

def test_primary_proton_tables_and_actual_bin_edges():
    j=json.loads((Path(__file__).resolve().parents[1]/'manifests/e866_pp_bins.json').read_text());rows=j['rows']
    assert len(rows)==279 and len({r['id'] for r in rows})==279
    assert {r['table'] for r in rows}=={'5.3','5.4','5.5','5.6'}
    assert sum(r['selected_core'] for r in rows)==41
    assert all(r['target']=='proton' and r['qtlow']<r['qthigh'] for r in rows)
    assert rows[0]['value']==13.3 and rows[0]['stat']==1.56 and rows[0]['syst']==.88
    assert rows[0]['Emean']>rows[0]['E_cm_max_in_bin']
    # Preserve inconsistent source metadata rather than silently repairing it.
    assert sum(not r['Emean_within_bin_bound'] for r in rows)==253

@pytest.mark.parametrize('Q',[6.2,7.2,10.85,16.85])
def test_high_scale_matching_is_canonical_and_does_not_cross_bottom(Q):
    p=HighScaleProfile(Q)
    s=np.r_[0.,np.geomspace(1e-8,1e8,80),np.inf]
    assert np.allclose(p.mu(s)*p.b_eval(s),B0,rtol=2e-16)
    assert p.mu(0)==pytest.approx(Q)
    assert float(p.mu(np.inf))>4.92
    FixedFlavorWindow(5,4.92,172.5).validate(float(p.mu(np.inf)),Q)

def test_profile_rejects_threshold_crossing_instead_of_clipping():
    with pytest.raises(ValueError,match='Threshold crossing'):HighScaleProfile(5.01)

def test_np_CS_two_leg_factor_and_input_gradients():
    m=GaussianCSInput(.7,.05);Q=11.;b=.8
    # Pair boundary exp(-lambda0*b^2/2) times two single-leg CS factors.
    expected=np.exp(-m.lambda0*b*b/2)*np.exp(-2*m.g*b*b*np.log(Q/m.Q0))
    assert np.exp(-m.variance(Q)*b*b/4)==pytest.approx(expected)
    k=np.linspace(0,5,50);cols=m.native_columns(k,Q);h=1e-6
    for column,p in [(1,'lambda0'),(2,'g')]:
        pars=dict(lambda0=.7,g=.05);pars[p]+=h;plus=GaussianCSInput(**pars).native_columns(k,Q)[:,0]
        pars[p]-=2*h;minus=GaussianCSInput(**pars).native_columns(k,Q)[:,0]
        assert np.max(abs(cols[:,column]-(plus-minus)/(2*h)))<1e-10

def test_leading_power_fraction_map():
    xf=np.array([-.04,.1,.3]);a,b=parton_x(8.,xf)
    assert np.allclose(a-b,xf)
    assert np.allclose(a*b,(8/38.8)**2)

def test_invariant_average_is_not_integrated_annular_yield():
    grid=RadialGrid(5.,1000);constant=np.ones((1000,2))*[3.,7.]
    assert np.allclose(pt_average(grid,constant,.013,1.246),[3.,7.])
    assert np.allclose(pt_average(grid,constant,.013,1.246,'annular'),[3.,7.])
    values=grid.k[:,None]**2
    assert pt_average(grid,values,0.,1.)[0]==pytest.approx(1/3,abs=3e-6)
    assert pt_average(grid,values,0.,1.,'annular')[0]==pytest.approx(.5,abs=1e-5)

class AnalyticPDF:
    def f(self,pid,x,mu):
        x=np.asarray(x);return (1+.01*pid)*x**.3*(1-x)**3*(1+.1*np.log(mu))

def test_vectorized_canonical_matching_matches_existing_endpoint_API():
    pdf=AnalyticPDF();rg=make_rg();mu=8.;xs=np.array([.08,.24,.6]);blocks=canonical_blocks(pdf,rg,xs,mu,160)
    a=float(rg.values(np.log(mu/rg.mu0))[0])
    for i,x in enumerate(xs):
        for f,pid in enumerate(FLAVORS):
            m=match_nlo(x,lambda y:pdf.f(pid,y,mu),lambda y:pdf.f(21,y,mu),a,order=160)
            assert np.allclose(blocks[i,f],[m.born,m.qq,m.qg],rtol=1e-13,atol=1e-15)

def test_single_energy_log_x_width_and_CS_slope_have_an_exact_null_direction():
    # With lambda(x)=a+b ln(1/x), v=2a+b ln(s/Q^2)+8g ln(Q/Q0).
    Q=np.array([6.5,8.,11.5,15.]);s=38.8**2;Q0=5.
    J=np.stack([np.full_like(Q,2),np.log(s/Q**2),8*np.log(Q/Q0)],axis=1)
    null=np.array([-4*np.log(s/Q0**2),8,2])
    assert np.max(abs(J@null))<1e-13
    assert np.linalg.matrix_rank(J)==2
