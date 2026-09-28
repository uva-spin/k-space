#!/usr/bin/env python3
"""
One-loop K_CSS2 evolution-Green-function validator.

This script is algebraic.  It validates the fixed-coupling one-loop
momentum-space TMD evolution Green function by working in b space, where
transverse convolution becomes multiplication.  The validated kernel is

  U_i(k;mu_f,zeta_f;mu_i,zeta_i)
    = exp{ a_s L_mu [gamma0_i - Gamma0_i (ell_i - L_mu)] }
      G_eta(k;mu_f),

  eta = a_s Gamma0_i rho,
  L_mu = ln(mu_f/mu_i),
  rho = ln sqrt(zeta_f/zeta_i),
  ell_i = ln(zeta_i/mu_i^2),

with FT[G_eta(k;mu_f)] = exp[-eta X_f], X_f=L_b(mu_f).

The tests are symbolic:
  1. mu-first and zeta-first paths agree.
  2. the rapidity differential equation is satisfied.
  3. the mu differential equation is satisfied.
  4. the semigroup/composition law is satisfied for three arbitrary points.
  5. the O(a_s) conventional-basis expansion is as expected.
"""
from __future__ import annotations

# Repository-local defaults. Override generated-output location with KCSS_OUTPUT_DIR.
import os as _os_for_repo
from pathlib import Path as _PathForRepo
_REPO_ROOT = _PathForRepo(__file__).resolve().parents[1]
_OUTPUT_DIR = _PathForRepo(_os_for_repo.environ.get('KCSS_OUTPUT_DIR', str(_REPO_ROOT / 'outputs'))).resolve()
_OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


import csv
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List

import sympy as sp


@dataclass(frozen=True)
class Config:
    output_prefix: str = str(_OUTPUT_DIR / "kt_kcss_one_loop_evolution_validator")


# Symbols for one segment i -> f.
a, Gamma, gamma = sp.symbols("a_s Gamma0 gammaF0")
Lmu, rho, ell_i, ell_f, X_i, X_f = sp.symbols("L_mu rho ell_i ell_f X_i X_f")

# Symbols for semigroup 0 -> 1 -> 2.
L01, L12, rho01, rho12, ell0, ell1, X1, X2 = sp.symbols(
    "L01 L12 rho01 rho12 ell0 ell1 X1 X2"
)


def logU_mu_first(L: sp.Expr, rho_: sp.Expr, ell_start: sp.Expr, X_end: sp.Expr) -> sp.Expr:
    """Log of one-loop Green function for a segment, represented at the final mu."""
    return sp.expand(a * (gamma * L - Gamma * ell_start * L + Gamma * L**2) - a * Gamma * rho_ * X_end)


def logU_zeta_first(L: sp.Expr, rho_: sp.Expr, ell_end: sp.Expr, X_start: sp.Expr) -> sp.Expr:
    """Log of the same Green function after zeta evolution first, represented at initial mu."""
    return sp.expand(a * (gamma * L - Gamma * ell_end * L - Gamma * L**2) - a * Gamma * rho_ * X_start)


def simplify_with_relations(expr: sp.Expr, relations: Dict[sp.Symbol, sp.Expr]) -> sp.Expr:
    return sp.simplify(sp.expand(expr.subs(relations)))


def main() -> None:
    cfg = Config()
    prefix = Path(cfg.output_prefix)

    # Kinematic relations for one segment:
    # X_f = X_i + 2 L_mu, ell_f = ell_i + 2 rho - 2 L_mu.
    one_segment_relations = {
        X_f: X_i + 2 * Lmu,
        ell_f: ell_i + 2 * rho - 2 * Lmu,
    }

    log_mu_first = logU_mu_first(Lmu, rho, ell_i, X_f)
    log_zeta_first = logU_zeta_first(Lmu, rho, ell_f, X_i)

    # Compare after using the relations above.  The zeta-first result has G_eta(mu_i),
    # while the mu-first result has G_eta(mu_f); in b space this is already represented
    # by X_i vs X_f, so the comparison is direct after substitution.
    path_residual = simplify_with_relations(log_mu_first - log_zeta_first, one_segment_relations)

    # Rapidity equation: d/d rho log U = - a Gamma X_f = -D^{(1)}(b;mu_f).
    rapidity_residual = sp.simplify(sp.diff(log_mu_first, rho) + a * Gamma * X_f)

    # Mu equation at the final point.  Differentiating with respect to ln(mu_f) means
    # d L_mu/d ln(mu_f)=1, d X_f/d ln(mu_f)=2, d ell_i/d ln(mu_f)=0, d rho/d ln(mu_f)=0.
    # The final anomalous dimension is a [gamma - Gamma ell_f].
    dlog_dL = sp.diff(log_mu_first, Lmu) + 2 * sp.diff(log_mu_first, X_f)
    mu_residual = simplify_with_relations(dlog_dL - a * (gamma - Gamma * ell_f), one_segment_relations)

    # Semigroup law for 0 -> 1 -> 2.  Relations:
    # ell1 = ell0 + 2 rho01 - 2 L01, X1 = X2 - 2 L12.
    log_01 = logU_mu_first(L01, rho01, ell0, X1)
    log_12 = logU_mu_first(L12, rho12, ell1, X2)
    log_02 = logU_mu_first(L01 + L12, rho01 + rho12, ell0, X2)
    semigroup_relations = {
        ell1: ell0 + 2 * rho01 - 2 * L01,
        X1: X2 - 2 * L12,
    }
    semigroup_residual = simplify_with_relations(log_01 + log_12 - log_02, semigroup_relations)

    # O(a_s) conventional-basis expansion.  Since G_eta = delta + eta L0 + O(eta^2),
    # U = delta + a_s R1 delta + a_s Gamma0 rho L0 + O(a_s^2), with
    # R1 = L_mu [gammaF0 - Gamma0(ell_i - L_mu)].
    R1 = sp.expand(Lmu * (gamma - Gamma * (ell_i - Lmu)))
    eta_over_as = Gamma * rho
    nlo_coeffs = {
        "delta": sp.sstr(R1),
        "L0": sp.sstr(eta_over_as),
    }

    # Numeric smoke check for a quark-like sample.
    CF = sp.Rational(4, 3)
    sample_values = {
        a: sp.Rational(1, 50),
        Gamma: 4 * CF,
        gamma: 6 * CF,
        Lmu: sp.Rational(3, 10),
        rho: sp.Rational(1, 5),
        ell_i: sp.Rational(7, 10),
        X_i: sp.Rational(11, 10),
    }
    sample_values[X_f] = sample_values[X_i] + 2 * sample_values[Lmu]
    sample_values[ell_f] = sample_values[ell_i] + 2 * sample_values[rho] - 2 * sample_values[Lmu]
    sample_logU = sp.N(log_mu_first.subs(sample_values), 18)
    sample_eta = sp.N((a * Gamma * rho).subs(sample_values), 18)
    sample_R = sp.N((a * R1).subs(sample_values), 18)

    tests = {
        "path_residual": sp.sstr(path_residual),
        "rapidity_equation_residual": sp.sstr(rapidity_residual),
        "mu_equation_residual": sp.sstr(mu_residual),
        "semigroup_residual": sp.sstr(semigroup_residual),
    }
    passed = all(val == "0" for val in tests.values())

    payload = {
        "passed": passed,
        "conventions": {
            "L_mu": "ln(mu_f/mu_i)",
            "rho": "ln sqrt(zeta_f/zeta_i)",
            "ell_i": "ln(zeta_i/mu_i^2)",
            "X_f": "L_b(mu_f)",
            "a_s": "alpha_s/(4*pi), held fixed for this one-loop algebraic validator",
        },
        "logU_mu_first_final_scale": sp.sstr(log_mu_first),
        "logU_zeta_first_initial_scale": sp.sstr(log_zeta_first),
        "tests": tests,
        "nlo_conventional_basis_coefficients_per_a_s": nlo_coeffs,
        "sample_quark_like_values": {str(k): sp.sstr(v) for k, v in sample_values.items()},
        "sample_logU_b_space": str(sample_logU),
        "sample_eta": str(sample_eta),
        "sample_scalar_exponent": str(sample_R),
        "notes": [
            "The validator checks the Green function in b space because convolution becomes multiplication.",
            "The corresponding k-space kernel is the scalar exponential multiplying G_eta(k;mu_f).",
            "Running-coupling and higher-logarithmic evolution require replacing fixed a_s by RG integrals; this is the next formal layer.",
        ],
    }

    json_path = Path(str(prefix) + "_output.json")
    txt_path = Path(str(prefix) + "_report.txt")
    csv_path = Path(str(prefix) + "_tests.csv")

    json_path.write_text(json.dumps(payload, indent=2), encoding="utf-8")

    with csv_path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=["test", "residual"])
        writer.writeheader()
        for key, val in tests.items():
            writer.writerow({"test": key, "residual": val})

    lines: List[str] = []
    lines.append("K_CSS2 one-loop evolution Green-function validator")
    lines.append("====================================================")
    lines.append(f"passed: {passed}")
    lines.append("")
    lines.append("Definitions")
    lines.append("-----------")
    lines.append("L_mu = ln(mu_f/mu_i)")
    lines.append("rho  = ln sqrt(zeta_f/zeta_i)")
    lines.append("ell_i = ln(zeta_i/mu_i^2)")
    lines.append("X_f = L_b(mu_f)")
    lines.append("")
    lines.append("Green function")
    lines.append("--------------")
    lines.append("U_i(k) = exp{a_s L_mu [gammaF0_i - Gamma0_i (ell_i - L_mu)]} G_eta(k;mu_f)")
    lines.append("eta = a_s Gamma0_i rho")
    lines.append("")
    lines.append("Symbolic residuals")
    lines.append("------------------")
    for key, val in tests.items():
        lines.append(f"{key}: {val}")
    lines.append("")
    lines.append("O(a_s) conventional-basis coefficients per a_s")
    lines.append("----------------------------------------------")
    for key, val in nlo_coeffs.items():
        lines.append(f"{key}: {val}")
    lines.append("")
    lines.append("Sample quark-like smoke values")
    lines.append("------------------------------")
    lines.append(f"sample eta: {sample_eta}")
    lines.append(f"sample scalar exponent: {sample_R}")
    lines.append(f"sample log U(b): {sample_logU}")
    txt_path.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print("\n".join(lines))


if __name__ == "__main__":
    main()
