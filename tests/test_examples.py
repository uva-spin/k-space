"""Software checks for the documented example clients; separate from physics suites."""
from __future__ import annotations
import importlib.util
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

import numpy as np
import pytest

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from resummed_operator import Coefficients,CanonicalRG

@pytest.fixture
def output_dir():
    parent=ROOT/'outputs';parent.mkdir(exist_ok=True)
    path=Path(tempfile.mkdtemp(prefix='example-test-',dir=parent))
    yield path
    shutil.rmtree(path)


def run_example(name,*options):
    env=os.environ.copy();env.update(OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',PYTHONDONTWRITEBYTECODE='1')
    return subprocess.run([sys.executable,str(ROOT/'examples'/name),*map(str,options)],cwd=ROOT,
                          text=True,capture_output=True,timeout=90,env=env)


def test_canonical_example_matches_operator(output_dir):
    proc=run_example('canonical_evolution.py','--out',output_dir)
    assert proc.returncode==0,proc.stdout+proc.stderr
    report=json.loads((output_dir/'canonical_evolution.json').read_text())
    samples=report['samples'];rg=CanonicalRG(Coefficients.load(ROOT/'data/coefficients_demo.json'))
    expected=rg.v(np.array([r['s_GeV_minus2'] for r in samples]),20.,400.)
    assert np.allclose([r['canonical_leg'] for r in samples],expected,rtol=2e-11,atol=2e-13)
    assert samples[0]['D_c'] is None and samples[0]['mu_c_GeV'] is None
    assert report['scope']['physical_N3LLprime_prediction_validated'] is False


def test_canonical_example_exercises_physical_n3ll_evolution_input(output_dir):
    coefficient=ROOT/'data/coefficients_css2_nf5_n3ll.json'
    proc=run_example('canonical_evolution.py','--coefficients',coefficient,'--out',output_dir)
    assert proc.returncode==0,proc.stdout+proc.stderr
    report=json.loads((output_dir/'canonical_evolution.json').read_text())
    assert report['coefficient_metadata']['evolution_accuracy']=='N3LL'
    assert report['coefficient_metadata']['hard_and_ope_order']=='strict NLO'
    assert report['scope']['physical_N3LL_evolution_input'] is True
    assert report['scope']['physical_N3LLprime_prediction_validated'] is False
    assert 'N3LL; boundaries: strict NLO' in proc.stdout

@pytest.mark.parametrize('arguments',[['--s','-1'],['--s','nan'],['--mu','1'],['--zeta','0'],['--a0','0'],['--mu','inf']])
def test_canonical_rejects_bad_arguments(arguments,output_dir):
    assert run_example('canonical_evolution.py',*arguments,'--out',output_dir).returncode!=0

@pytest.mark.parametrize('reference',[False,True])
def test_grid_example_scopes_and_shapes(reference,output_dir):
    options=['--grid','5','--degrees','4','8','--out',output_dir]
    if reference: options.append('--reference')
    proc=run_example('momentum_grid.py',*options)
    assert proc.returncode==0,proc.stdout+proc.stderr
    report=json.loads((output_dir/'momentum_grid.json').read_text())
    assert report['reference_computed'] is reference
    assert report['degrees']==[4,8]
    assert report['scope']['continuum_grid_error_certified'] is False
    assert report['scope']['data_fit_performed'] is False
    array=np.genfromtxt(output_dir/'momentum_grid.csv',delimiter=',',names=True)
    assert array.shape==(25,)
    for row in report['results']:
        assert abs(row['input_L2']-1)<1e-12
        if reference:
            assert row['L2_same_grid_error']<=row['finite_matrix_bound']+2e-10
            assert row['recurrence_oracle_residual']<2e-10
        else:assert row['L2_same_grid_error'] is None

@pytest.mark.parametrize('arguments',[['--grid','2'],['--grid','26','--reference'],['--halfwidth','0'],['--degrees','0'],['--degrees','513']])
def test_grid_rejects_bad_arguments(arguments,output_dir):
    assert run_example('momentum_grid.py',*arguments,'--out',output_dir).returncode!=0


def test_examples_reject_source_output():
    assert run_example('canonical_evolution.py','--out',ROOT).returncode!=0


def test_documentation_links_and_fences():
    spec=importlib.util.spec_from_file_location('documentation_check',ROOT/'tools/check_documentation.py')
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    assert module.scan(ROOT)==[]
