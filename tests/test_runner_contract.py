"""Software-contract tests; separate from the scientific regression suites."""
import importlib.util
from pathlib import Path
import pytest
import subprocess
ROOT=Path(__file__).resolve().parents[1]
def load(name,path):
    spec=importlib.util.spec_from_file_location(name,path)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    return module
runner=load('verification_runner',ROOT/'run_all_checks.py')
policy=load('verification_policy',ROOT/'tools/check_repository_policy.py')

def test_explicit_true():
    assert runner.validate_report({'passed':True},['passed'])==(True,None,None)
@pytest.mark.parametrize('value',[False,'True','False',1,0,None])
def test_truthy_values_rejected(value):
    assert not runner.validate_report({'passed':value},['passed'])[0]
def test_missing_field_fails():
    with pytest.raises(ValueError):runner.validate_report({},['passed'])
def test_individual_failure_overrides_aggregate():
    assert runner.validate_report({'passed':True,'tests':[{'passed':False}]},['passed'])==(False,0,1)
def test_empty_tests_rejected():
    assert not runner.validate_report({'passed':True,'tests':[]},['passed'])[0]
def test_malformed_tests_rejected():
    assert not runner.validate_report({'passed':True,'tests':[{'passed':'True'}]},['passed'])[0]
def test_cleanup_cannot_target_root():
    with pytest.raises(ValueError):runner.checked_output_dir(str(ROOT))
def test_cleanup_cannot_target_outputs_root():
    with pytest.raises(ValueError):runner.checked_output_dir(str(ROOT/'outputs'))
def test_cleanup_path_is_local():
    assert runner.checked_output_dir(str(ROOT/'outputs/test'))==ROOT/'outputs/test'
def test_manuscript_policy(tmp_path):
    (tmp_path/'draft.tex').write_text('separate writing workspace\n')
    assert any('draft.tex' in s for s in policy.scan(tmp_path))
def test_archive_policy(tmp_path):
    (tmp_path/'notes.zip').write_bytes(b'data')
    assert any('notes.zip' in s for s in policy.scan(tmp_path))

@pytest.mark.parametrize('name',['AGENTS.md','CLAUDE.md','.codex/session.json',
                                '.cursor/rules/task.md','launch_codex.sh','console.log'])
def test_local_work_artifacts_rejected(tmp_path,name):
    p=tmp_path/name;p.parent.mkdir(parents=True,exist_ok=True);p.write_text('{}\n')
    assert any(name in s for s in policy.scan(tmp_path))

def test_dependencies_docs_and_local_outputs_allowed(tmp_path):
    for name in ['requirements.txt','requirements-dev.txt','README.md',
                 'outputs/latest/report.txt','outputs/latest/console.log','data/measurement.txt','docs/scientific_report.txt']:
        p=tmp_path/name;p.parent.mkdir(parents=True,exist_ok=True);p.write_text('content\n')
    assert policy.scan(tmp_path)==[]

@pytest.mark.parametrize('name',['outputs/latest/report.txt','validators/__pycache__/cached.pyc'])
def test_tracked_generated_artifacts_rejected(tmp_path,name):
    subprocess.run(['git','init','-q',str(tmp_path)],check=True)
    p=tmp_path/name;p.parent.mkdir(parents=True,exist_ok=True);p.write_text('content\n')
    subprocess.run(['git','-C',str(tmp_path),'add',name],check=True)
    assert any(name in s for s in policy.scan(tmp_path))

import json
import shutil
import sys

@pytest.fixture
def isolated_runner(tmp_path):
    for name in ['run_all_checks.py','tools/execute_suite.py']:
        target=tmp_path/name;target.parent.mkdir(parents=True,exist_ok=True)
        shutil.copyfile(ROOT/name,target)
    (tmp_path/'validators').mkdir()
    (tmp_path/'manifests').mkdir()
    spec=dict(id='probe',script='validators/probe.py',output='report.json',
              pass_key=['passed'],invocation='environment',category='software_contract',
              description='Runner failure injection')
    (tmp_path/'manifests/verification_suite.json').write_text(json.dumps(
        dict(suites=[spec],required_suite_ids=['probe'])))
    return tmp_path

@pytest.mark.parametrize('mode',['exit_failure','missing','malformed','individual_failure',
                                 'skipped','unsupported_claim','success'])
def test_process_contract_without_manuscripts(isolated_runner,mode):
    root=isolated_runner
    stale=root/'outputs/latest/probe/report.json';stale.parent.mkdir(parents=True)
    stale.write_text('{"passed": true}')
    payload={'passed':True,'tests':[{'passed':True}]}
    if mode=='individual_failure':payload['tests'][0]['passed']=False
    if mode=='skipped':payload['tests'][0]['status']='skipped'
    if mode=='unsupported_claim':payload['physical_N3LLprime_prediction_validated']=True
    body="import os, sys\nfrom pathlib import Path\np=Path(os.environ['KCSS_OUTPUT_DIR'])/'report.json'\n"
    if mode=='exit_failure':body+='sys.exit(7)\n'
    elif mode=='malformed':body+="p.write_text('{broken')\n"
    elif mode!='missing':body+='p.write_text('+repr(json.dumps(payload))+')\n'
    (root/'validators/probe.py').write_text(body)
    result=subprocess.run([sys.executable,str(root/'run_all_checks.py')],
                          cwd=root,capture_output=True,text=True)
    assert not list(root.rglob('*.tex')) and not list(root.rglob('*.bib'))
    assert result.returncode==(0 if mode=='success' else 1),result.stdout+result.stderr
    summary=json.loads((root/'outputs/latest/run_all_checks_summary.json').read_text())
    assert summary['all_suites_passed'] is (mode=='success')
    if mode in {'exit_failure','missing'}:assert not stale.exists()

@pytest.mark.parametrize('field',['suite_id','run_id','report_sha256'])
def test_wrong_execution_identity_rejected(tmp_path,field):
    report=tmp_path/'report.json';report.write_text('{"passed":true}')
    envelope=dict(schema='k-space-suite-execution/1',suite_id='probe',run_id='current',
                  report_sha256=runner.sha(report))
    envelope[field]='wrong'
    (tmp_path/'execution.json').write_text(json.dumps(envelope))
    with pytest.raises(ValueError):runner.validate_execution(tmp_path,report,'probe','current')

@pytest.mark.parametrize('mutation',['removed_suite','unregistered_script','empty','duplicate'])
def test_required_discovery_contract(isolated_runner,mutation):
    root=isolated_runner;(root/'validators/probe.py').write_text('')
    manifest=json.loads((root/'manifests/verification_suite.json').read_text())
    if mutation=='removed_suite':manifest['suites']=[]
    if mutation=='unregistered_script':(root/'validators/extra.py').write_text('')
    if mutation=='empty':manifest=dict(suites=[],required_suite_ids=[])
    if mutation=='duplicate':manifest['suites']*=2
    with pytest.raises(ValueError):runner.validate_manifest(manifest,root)

def test_nested_assertion_rejected():
    assert not runner.validate_report({'passed':True,'nested':{'passed':False}},['passed'])[0]

def test_generated_agent_artifact_rejected(tmp_path):
    p=tmp_path/'outputs/latest/AGENTS.md';p.parent.mkdir(parents=True);p.write_text('instructions\n')
    assert any('agent instructions' in s for s in policy.scan(tmp_path))

def test_scientific_pdf_plot_allowed(tmp_path):
    p=tmp_path/'figures/operator.pdf';p.parent.mkdir();p.write_bytes(b'%PDF scientific figure')
    assert policy.scan(tmp_path)==[]

def test_compound_build_suffix_rejected(tmp_path):
    (tmp_path/'draft.run.xml').write_text('<build/>\n')
    assert any('draft.run.xml' in s for s in policy.scan(tmp_path))

def test_optional_acquisition_ledger_is_repeatable(tmp_path):
    # Pre-stage every requested file, so this workflow regression runs offline.
    root=tmp_path/'sources'
    members={
        '2006.05329v2': ['TMDPDF.m','PTBF.m','PTBF_ZExpansion.m','PTBF_ZbExpansion.m',
                       'PT_SoftFunction.m','PT_SoftFunction_Renormalized.m','info.txt',
                       'e-print-2006.05329v2.tar'],
        '2012.03256v2': ['BeamFunction.m','BeamfunctionN.m','FFMatchingKernels.m',
                       'FFMatchingKernelsN.m','TMDFF.m','TMDFFN.m','TMDPDF.m','TMDPDFN.m',
                       'ancillary_readme.pdf','resummedFFsinglets.m','softfunction.m',
                       'e-print-2012.03256v2.tar']}
    for directory,names in members.items():
        for name in names:
            path=root/directory/name;path.parent.mkdir(parents=True,exist_ok=True)
            path.write_bytes(b'offline workflow fixture')
    script=ROOT/'scripts/stage_public_coefficients.sh'
    subprocess.run(['bash',str(script),str(root)],check=True,capture_output=True)
    first=(root/'SHA256SUMS.txt').read_bytes()
    subprocess.run(['bash',str(script),str(root)],check=True,capture_output=True)
    assert (root/'SHA256SUMS.txt').read_bytes()==first
    assert b'SHA256SUMS.txt' not in first
