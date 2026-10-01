import json,math
from pathlib import Path
import pytest
ROOT=Path(__file__).resolve().parents[1]

def test_primary_e605_rows_are_unique_and_keep_error_components():
    j=json.loads((ROOT/'manifests/e605_primary_points.json').read_text())
    rows=j['rows']
    assert len(rows)==54 and len({r['source_id'] for r in rows})==54
    assert sum(r['selected_small_qt_core'] for r in rows)==52
    assert j['target_proton_fraction'] is None
    assert j['production_prediction_allowed'] is False
    for r in rows:
        assert r['qt_bin_edges'] is None and r['xf_bin_edges'] is None
        assert abs(r['CS_rescaling_residual'])<1e-6
        assert r['quoted_systematics'][0]['symerror']=='15.0%'
        assert r['quoted_systematics'][1]['symerror']=='5.0%'
        assert r['xf']==-.1
    first=rows[0]
    assert first['source_table_doi']=='10.17182/hepdata.22831.v1/t17'
    assert first['invariant_pb_GeV_minus2_per_nucleon']==.785
    assert first['quoted_absolute_error']==.115
