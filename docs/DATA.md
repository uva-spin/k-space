# Preserved Drell–Yan CSV tables

`Data_DY_test/` contains 21 CSV tables retained byte-for-byte from the supplied
repository snapshot. They provide possible inputs for later phenomenological
work; neither the registered scientific verification suites nor the runnable
operator examples load them.

## Inventory

The following counts are file-level properties obtained with Python's CSV reader,
not claims that all rows pass a physics selection. The exact column names and
SHA256 values are recorded in [`manifests/dy_tables.json`](../manifests/dy_tables.json).

| File | Data rows | Columns |
|---|---:|---:|
| [`ATLAS_13.csv`](../Data_DY_test/ATLAS_13.csv) | 78 | 47 |
| [`ATLAS_7.csv`](../Data_DY_test/ATLAS_7.csv) | 78 | 45 |
| [`ATLAS_8.csv`](../Data_DY_test/ATLAS_8.csv) | 160 | 44 |
| [`CDF_RUN_1.csv`](../Data_DY_test/CDF_RUN_1.csv) | 50 | 28 |
| [`CDF_RUN_2.csv`](../Data_DY_test/CDF_RUN_2.csv) | 82 | 28 |
| [`CMS_13.csv`](../Data_DY_test/CMS_13.csv) | 170 | 37 |
| [`CMS_7.csv`](../Data_DY_test/CMS_7.csv) | 18 | 43 |
| [`CMS_8.csv`](../Data_DY_test/CMS_8.csv) | 18 | 46 |
| [`D0_RUN_1.csv`](../Data_DY_test/D0_RUN_1.csv) | 26 | 29 |
| [`D0_RUN_2.csv`](../Data_DY_test/D0_RUN_2.csv) | 24 | 31 |
| [`D0_RUN_2N.csv`](../Data_DY_test/D0_RUN_2N.csv) | 13 | 45 |
| [`E288_200.csv`](../Data_DY_test/E288_200.csv) | 89 | 22 |
| [`E288_300.csv`](../Data_DY_test/E288_300.csv) | 118 | 22 |
| [`E288_400.csv`](../Data_DY_test/E288_400.csv) | 177 | 22 |
| [`E605.csv`](../Data_DY_test/E605.csv) | 93 | 21 |
| [`E772.csv`](../Data_DY_test/E772.csv) | 102 | 21 |
| [`LHCb_13.csv`](../Data_DY_test/LHCb_13.csv) | 14 | 37 |
| [`LHCb_7.csv`](../Data_DY_test/LHCb_7.csv) | 14 | 37 |
| [`LHCb_8.csv`](../Data_DY_test/LHCb_8.csv) | 14 | 37 |
| [`PHENIX200.csv`](../Data_DY_test/PHENIX200.csv) | 12 | 28 |
| [`STAR_510.csv`](../Data_DY_test/STAR_510.csv) | 23 | 30 |

## Interpretation before fitting

The tables preserve their original labels and values. For example, the
`E288_200.csv` header includes `qT`, `A`, `PreFactor`, `CS`, `dA`, `error`,
`sysNorm`, `sysP2P`, rapidity and mass-bin columns, `SqrtS`, `xF`, `x1`, `x2`,
`BeamE`, and `frac_error`. Those column names alone do not establish physical
units, preprocessing conventions, the experimental covariance, or the correct
use of a normalization uncertainty. Other files have different schemas.

Before using any table in an extraction, establish its original publication,
kinematics, units, normalization and acceptance conventions, treatment of pointwise
and correlated uncertainties, and any nuclear or fiducial corrections. Source
attribution and redistribution conditions for these CSV values are not resolved
by the software license. No missing covariance or experimental metadata is filled
in by the verification scripts.

This inventory is intended to preserve useful data without representing it as a
completed or validated data-analysis input. The coefficient fixture used by the
operator examples is a different file, described in
[COEFFICIENTS.md](COEFFICIENTS.md).
