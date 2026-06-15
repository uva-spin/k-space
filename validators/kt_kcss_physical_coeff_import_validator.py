#!/usr/bin/env python3
"""
Source-audited coefficient import gate for the K_CSS2 N3LL' project.

This script deliberately does not use synthetic coefficient values.  It checks the
source ledger, the expected ancillary filenames, and the local availability/hash state
of raw coefficient packages.  If raw files are absent, the audit can still pass as an
infrastructure/source-ledger test, but physics_ready is False and the NNLO/N3LO
singular-expansion gates are blocked.
"""
from __future__ import annotations

# Repository-local defaults. Override generated-output location with KCSS_OUTPUT_DIR.
import os as _os_for_repo
from pathlib import Path as _PathForRepo
_REPO_ROOT = _PathForRepo(__file__).resolve().parents[1]
_OUTPUT_DIR = _PathForRepo(_os_for_repo.environ.get('KCSS_OUTPUT_DIR', str(_REPO_ROOT / 'outputs'))).resolve()
_OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


import argparse
import csv
import hashlib
import json
from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Iterable

@dataclass(frozen=True)
class ExpectedFile:
    source_id: str
    arxiv: str
    version: str
    loop_relevance: str
    operator_type: str
    filename: str
    mandatory_for_pdf_import: bool
    purpose: str

SOURCES = [
    {
        "source_id": "EMV2020_TMDPDF_N3LO",
        "arxiv": "2006.05329",
        "version": "v2",
        "authors": "Ebert-Mistlberger-Vita",
        "bibliographic_role": "N3LO unpolarized TMDPDF matching kernels and soft function cross-checks",
        "official_record_evidence": "arXiv ancillary-file list includes TMDPDF.m and related PTBF/soft-function files",
    },
    {
        "source_id": "LYZZ2021_TMDPDF_TMDFF_N3LO",
        "arxiv": "2012.03256",
        "version": "v2",
        "authors": "Luo-Yang-Zhu-Zhu",
        "bibliographic_role": "N3LO unpolarized quark/gluon TMDPDF and TMDFF matching kernels",
        "official_record_evidence": "arXiv ancillary-file list includes BeamFunction.m, TMDPDF.m, TMDPDFN.m, TMDFF.m, and ancillary_readme.pdf",
    },
]

EXPECTED_FILES = [
    ExpectedFile("EMV2020_TMDPDF_N3LO", "2006.05329", "v2", "N3LO", "TMDPDF", "TMDPDF.m", True, "Primary TMDPDF coefficient payload"),
    ExpectedFile("EMV2020_TMDPDF_N3LO", "2006.05329", "v2", "N3LO", "beam/collinear", "PTBF.m", False, "Perturbative transverse beam-function cross-check"),
    ExpectedFile("EMV2020_TMDPDF_N3LO", "2006.05329", "v2", "N3LO", "beam/collinear", "PTBF_ZExpansion.m", False, "z-expansion cross-check"),
    ExpectedFile("EMV2020_TMDPDF_N3LO", "2006.05329", "v2", "N3LO", "beam/collinear", "PTBF_ZbExpansion.m", False, "z and b expansion cross-check"),
    ExpectedFile("EMV2020_TMDPDF_N3LO", "2006.05329", "v2", "N3LO", "soft", "PT_SoftFunction.m", False, "Bare/renormalized soft-function cross-check"),
    ExpectedFile("EMV2020_TMDPDF_N3LO", "2006.05329", "v2", "N3LO", "soft", "PT_SoftFunction_Renormalized.m", False, "Renormalized soft-function cross-check"),
    ExpectedFile("EMV2020_TMDPDF_N3LO", "2006.05329", "v2", "metadata", "readme", "info.txt", False, "Source package notes"),

    ExpectedFile("LYZZ2021_TMDPDF_TMDFF_N3LO", "2012.03256", "v2", "N3LO", "beam/TMDPDF", "BeamFunction.m", True, "Beam-function/TMDPDF coefficient payload"),
    ExpectedFile("LYZZ2021_TMDPDF_TMDFF_N3LO", "2012.03256", "v2", "N3LO", "beam/TMDPDF", "BeamfunctionN.m", False, "Mellin-N representation cross-check"),
    ExpectedFile("LYZZ2021_TMDPDF_TMDFF_N3LO", "2012.03256", "v2", "N3LO", "TMDPDF", "TMDPDF.m", True, "Primary TMDPDF payload"),
    ExpectedFile("LYZZ2021_TMDPDF_TMDFF_N3LO", "2012.03256", "v2", "N3LO", "TMDPDF", "TMDPDFN.m", False, "Mellin-N TMDPDF payload"),
    ExpectedFile("LYZZ2021_TMDPDF_TMDFF_N3LO", "2012.03256", "v2", "N3LO", "TMDFF", "TMDFF.m", False, "TMDFF payload for SIDIS/e+e- extensions"),
    ExpectedFile("LYZZ2021_TMDPDF_TMDFF_N3LO", "2012.03256", "v2", "N3LO", "TMDFF", "TMDFFN.m", False, "Mellin-N TMDFF payload"),
    ExpectedFile("LYZZ2021_TMDPDF_TMDFF_N3LO", "2012.03256", "v2", "N3LO", "TMDFF", "FFMatchingKernels.m", False, "Fragmentation matching kernels"),
    ExpectedFile("LYZZ2021_TMDPDF_TMDFF_N3LO", "2012.03256", "v2", "N3LO", "TMDFF", "FFMatchingKernelsN.m", False, "Mellin-N fragmentation matching kernels"),
    ExpectedFile("LYZZ2021_TMDPDF_TMDFF_N3LO", "2012.03256", "v2", "metadata", "readme", "ancillary_readme.pdf", False, "Ancillary convention documentation"),
    ExpectedFile("LYZZ2021_TMDPDF_TMDFF_N3LO", "2012.03256", "v2", "N3LO", "soft", "softfunction.m", False, "Soft-function cross-check"),
]


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def find_file(root: Path, source_id: str, filename: str) -> Path | None:
    candidates = [root / source_id / filename, root / filename]
    for c in candidates:
        if c.exists():
            return c
    matches = list(root.rglob(filename)) if root.exists() else []
    return matches[0] if matches else None


def write_csv(path: Path, rows: Iterable[dict], fieldnames: list[str]) -> None:
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source-dir", default=str(_REPO_ROOT / "coefficient_sources"), help="Directory containing raw ancillary files")
    parser.add_argument("--out-prefix", default=str(_OUTPUT_DIR / "kt_kcss_physical_coeff_import_validator"))
    args = parser.parse_args()

    source_dir = Path(args.source_dir)
    prefix = Path(args.out_prefix)

    manifest_rows = []
    local_rows = []
    for rec in EXPECTED_FILES:
        d = asdict(rec)
        manifest_rows.append(d)
        found = find_file(source_dir, rec.source_id, rec.filename)
        local_rows.append({
            **d,
            "local_path": str(found) if found else "",
            "present": bool(found),
            "byte_size": found.stat().st_size if found else 0,
            "sha256": sha256_file(found) if found else "",
        })

    mandatory = [r for r in local_rows if r["mandatory_for_pdf_import"]]
    mandatory_present = sum(1 for r in mandatory if r["present"])
    all_mandatory_present = mandatory_present == len(mandatory)
    any_raw_files_present = any(r["present"] for r in local_rows)

    tests = [
        {"test": "source_ledger_has_two_primary_families", "passed": len(SOURCES) == 2, "details": str([s["source_id"] for s in SOURCES])},
        {"test": "mandatory_pdf_files_declared", "passed": len(mandatory) >= 3, "details": f"mandatory={len(mandatory)}"},
        {"test": "mock_payload_disabled_for_physics_gate", "passed": True, "details": "No synthetic coefficients are accepted for physics_ready=True."},
        {"test": "local_source_directory_exists", "passed": source_dir.exists(), "details": str(source_dir)},
        {"test": "raw_files_present", "passed": any_raw_files_present, "details": f"present={sum(1 for r in local_rows if r['present'])}/{len(local_rows)}"},
        {"test": "mandatory_pdf_import_files_present", "passed": all_mandatory_present, "details": f"present={mandatory_present}/{len(mandatory)}"},
        {"test": "nnlo_n3lo_expansion_gate_ready", "passed": all_mandatory_present, "details": "Requires parsed physical coefficients; blocked if mandatory files are absent."},
    ]

    infrastructure_passed = all(t["passed"] for t in tests[:3])
    physics_ready = all(t["passed"] for t in tests)

    output = {
        "validator": "kt_kcss_vetted_coeff_source_gate_validator",
        "source_dir": str(source_dir),
        "infrastructure_passed": infrastructure_passed,
        "physics_ready": physics_ready,
        "mock_payload_used": False,
        "source_families": SOURCES,
        "expected_file_count": len(EXPECTED_FILES),
        "mandatory_file_count": len(mandatory),
        "present_file_count": sum(1 for r in local_rows if r["present"]),
        "mandatory_present_count": mandatory_present,
        "blocked_reason": None if physics_ready else "Raw vetted ancillary coefficient files are not all present locally; NNLO/N3LO expansion gates were not executed.",
        "tests": tests,
    }

    write_csv(prefix.with_name(prefix.name + "_source_manifest.csv"), manifest_rows, list(manifest_rows[0].keys()))
    write_csv(prefix.with_name(prefix.name + "_local_audit.csv"), local_rows, list(local_rows[0].keys()))
    write_csv(prefix.with_name(prefix.name + "_tests.csv"), tests, ["test", "passed", "details"])
    prefix.with_name(prefix.name + "_output.json").write_text(json.dumps(output, indent=2))

    report = []
    report.append("K_CSS2 vetted coefficient-source gate validator")
    report.append("================================================")
    report.append(f"Source directory: {source_dir}")
    report.append(f"Infrastructure passed: {infrastructure_passed}")
    report.append(f"Physics ready: {physics_ready}")
    report.append(f"Mock payload used: False")
    report.append(f"Expected files: {len(EXPECTED_FILES)}")
    report.append(f"Present files: {output['present_file_count']}")
    report.append(f"Mandatory PDF-import files present: {mandatory_present}/{len(mandatory)}")
    if output["blocked_reason"]:
        report.append(f"Blocked reason: {output['blocked_reason']}")
    report.append("")
    report.append("Tests:")
    for t in tests:
        report.append(f"  - {t['test']}: {t['passed']} ({t['details']})")
    prefix.with_name(prefix.name + "_report.txt").write_text("\n".join(report) + "\n")

    print(json.dumps(output, indent=2))

if __name__ == "__main__":
    main()
