#!/usr/bin/env python3
"""Run the KCSS pure-kT formalism validation suite.

The suite is intentionally split into two categories:
  1. formal/infrastructure checks that must pass for the formalism paper, and
  2. physical N3LL' implementation gates that are expected to remain blocked until
     external NNLO/N3LO coefficient payloads are staged and converted.
"""
from __future__ import annotations

import argparse
import json
import os
import shutil
import subprocess
import sys
from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Callable, Optional

REPO_ROOT = Path(__file__).resolve().parent
VALIDATOR_DIR = REPO_ROOT / "validators"


@dataclass
class CheckResult:
    name: str
    script: str
    returncode: int
    formal_passed: bool
    expected_blocked: bool
    detail: str
    log_file: str


def load_json(path: Path) -> dict:
    if not path.exists():
        return {}
    return json.loads(path.read_text())


def bool_from_json(path: Path, key_path: tuple[str, ...], default: Optional[bool] = None) -> Optional[bool]:
    data = load_json(path)
    cur = data
    for key in key_path:
        if not isinstance(cur, dict) or key not in cur:
            return default
        cur = cur[key]
    return bool(cur)


def summarize_toypdf(path: Path) -> bool:
    data = load_json(path)
    tests = data.get("tests", {})
    if not tests:
        return False
    bools = [v for v in tests.values() if isinstance(v, bool)]
    numerics = [abs(v) for v in tests.values() if isinstance(v, (int, float)) and not isinstance(v, bool)]
    return all(bools) and all(v < 1e-10 for v in numerics)


def summarize_physical_gate(path: Path) -> tuple[bool, bool, str]:
    data = load_json(path)
    infra = bool(data.get("infrastructure_passed", False))
    physics_ready = bool(data.get("physics_ready", False))
    detail = f"infrastructure_passed={infra}; physics_ready={physics_ready}"
    # For this formalism repository, physics_ready=False is an expected blocked gate.
    return infra, not physics_ready, detail


def summarize_label_certificate(path: Path) -> tuple[bool, bool, str]:
    data = load_json(path)
    formal = bool(data.get("formal_prescription_ready", False))
    impl = bool(data.get("implementation_claim_allowed", False))
    detail = f"formal_prescription_ready={formal}; implementation_claim_allowed={impl}"
    # Formalism paper requires formal=True and implementation claim blocked until external payloads are present.
    return formal, not impl, detail


def one_loop_summary(outdir: Path) -> tuple[bool, bool, str]:
    data = load_json(outdir / "kt_kcss_one_loop_validator_report.json")
    tests = data.get("tests", {})
    passed = bool(tests.get("pass", False))
    detail = f"pass={passed}; plus_null={tests.get('plus_null_max_abs')}; additivity={tests.get('additivity_max_abs')}; offzero={tests.get('offzero_diag_max_abs_difference')}"
    return passed, False, detail


VALIDATORS: list[tuple[str, str, Callable[[Path], tuple[bool, bool, str]], str]] = []


def add_passed(name: str, script: str, json_name: str) -> None:
    def summarizer(outdir: Path, json_name: str = json_name) -> tuple[bool, bool, str]:
        data = load_json(outdir / json_name)
        passed = bool(data.get("passed", False))
        return passed, False, f"passed={passed}"
    VALIDATORS.append((name, script, summarizer, json_name))


VALIDATORS.append(("one-loop binned singular spectrum", "kt_kcss_one_loop_validator.py", one_loop_summary, "kt_kcss_one_loop_validator_report.json"))
add_passed("toy-PDF x-convolution", "kt_kcss_toypdf_xconv_validator.py", "kt_kcss_toypdf_xconv_validator_output.json")
add_passed("generator convolution algebra", "kt_kcss_generator_convolution_validator.py", "kt_kcss_generator_convolution_validator_output.json")
add_passed("one-loop evolution Green function", "kt_kcss_one_loop_evolution_validator.py", "kt_kcss_one_loop_evolution_validator_output.json")
add_passed("two-loop ingredient checks", "kt_kcss_two_loop_ingredients_validator.py", "kt_kcss_two_loop_ingredients_validator_output.json")
add_passed("running-coupling evolution", "kt_kcss_running_coupling_evolution_validator.py", "kt_kcss_running_coupling_evolution_validator_output.json")
add_passed("binned W-term prototype", "kt_kcss_binned_wterm_prototype_validator.py", "kt_kcss_binned_wterm_prototype_validator_output.json")
add_passed("accuracy ingredient manifest", "kt_kcss_accuracy_ingredient_validator.py", "kt_kcss_accuracy_ingredient_validator_output.json")
add_passed("N3LO mock matching import", "kt_kcss_n3lo_matching_import_validator.py", "kt_kcss_n3lo_matching_import_output.json")


def formalism_completion_summary(outdir: Path) -> tuple[bool, bool, str]:
    data = load_json(outdir / "kt_kcss_formalism_completion_validator_output.json")
    summary = data.get("summary", {})
    passed = bool(summary.get("passed", False))
    detail = f"passed={passed}; status={summary.get('formalism_note_status')}"
    return passed, False, detail


def toypdf_summary(outdir: Path) -> tuple[bool, bool, str]:
    passed = summarize_toypdf(outdir / "kt_kcss_toypdf_xconv_validator_output.json")
    return passed, False, f"passed={passed}"

# Override toy-PDF summarizer because it does not expose top-level passed.
VALIDATORS[1] = ("toy-PDF x-convolution", "kt_kcss_toypdf_xconv_validator.py", toypdf_summary, "kt_kcss_toypdf_xconv_validator_output.json")
VALIDATORS.append(("physical coefficient source gate", "kt_kcss_physical_coeff_import_validator.py", lambda out: summarize_physical_gate(out / "kt_kcss_physical_coeff_import_validator_output.json"), "kt_kcss_physical_coeff_import_validator_output.json"))
VALIDATORS.append(("N3LL-prime label certificate", "kt_kcss_n3llprime_label_certificate_validator.py", lambda out: summarize_label_certificate(out / "kt_kcss_n3llprime_label_certificate_output.json"), "kt_kcss_n3llprime_label_certificate_output.json"))
VALIDATORS.append(("formalism completion", "kt_kcss_formalism_completion_validator.py", formalism_completion_summary, "kt_kcss_formalism_completion_validator_output.json"))


def main() -> int:
    parser = argparse.ArgumentParser(description="Run all KCSS pure-kT formalism validators.")
    parser.add_argument("--output-dir", default=str(REPO_ROOT / "outputs" / "latest"), help="Directory for generated outputs")
    parser.add_argument("--clean", action="store_true", help="Remove output directory before running")
    parser.add_argument("--python", default=sys.executable, help="Python interpreter to use")
    args = parser.parse_args()

    outdir = Path(args.output_dir).resolve()
    if args.clean and outdir.exists():
        shutil.rmtree(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    env = os.environ.copy()
    env["KCSS_OUTPUT_DIR"] = str(outdir)
    env.setdefault("PYTHONPATH", str(REPO_ROOT))

    results: list[CheckResult] = []
    for name, script, summarizer, _ in VALIDATORS:
        script_path = VALIDATOR_DIR / script
        log_path = outdir / f"{script_path.stem}.log"
        print(f"[run] {name}: {script}")
        proc = subprocess.run(
            [args.python, str(script_path)],
            cwd=str(REPO_ROOT),
            env=env,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
        )
        log_path.write_text(proc.stdout)
        if proc.returncode != 0:
            formal_passed = False
            expected_blocked = False
            detail = f"returncode={proc.returncode}; see {log_path.name}"
        else:
            try:
                formal_passed, expected_blocked, detail = summarizer(outdir)
            except Exception as exc:
                formal_passed = False
                expected_blocked = False
                detail = f"could not summarize output: {exc!r}"
        results.append(CheckResult(name, script, proc.returncode, formal_passed, expected_blocked, detail, str(log_path.relative_to(REPO_ROOT))))

    formal_ok = all(r.returncode == 0 and r.formal_passed for r in results)
    blocked = [r for r in results if r.expected_blocked]
    summary = {
        "formal_suite_passed": formal_ok,
        "num_checks": len(results),
        "num_formal_passed": sum(1 for r in results if r.returncode == 0 and r.formal_passed),
        "expected_blocked_gates": [r.name for r in blocked],
        "implementation_N3LLprime_claim_allowed": False,
        "note": "The suite validates the formal prescription. Physical N3LL-prime implementation remains blocked until external NNLO/N3LO coefficient payloads are staged and the a_s^3 expansion gate passes.",
    }
    (outdir / "run_all_checks_summary.json").write_text(json.dumps({"summary": summary, "results": [asdict(r) for r in results]}, indent=2))
    lines = [
        "KCSS pure-kT formalism validation-suite summary",
        "=================================================",
        f"formal_suite_passed: {formal_ok}",
        f"checks: {summary['num_formal_passed']}/{summary['num_checks']}",
        "",
        "Individual checks:",
    ]
    for r in results:
        status = "PASS" if (r.returncode == 0 and r.formal_passed) else "FAIL"
        if r.expected_blocked:
            status += " (expected implementation gate blocked)"
        lines.append(f"  - {status}: {r.name} :: {r.detail}")
    lines.extend(["", summary["note"]])
    (outdir / "run_all_checks_report.txt").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))
    return 0 if formal_ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
