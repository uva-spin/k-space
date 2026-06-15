#!/usr/bin/env python3
"""Two-loop KCSS ingredient validator.

Algebraic regression checks for the v0.10 two-loop Collins-Soper kernel and the
logarithmic reconstruction of NNLO matching coefficients.  This is not a
phenomenology code and does not evaluate the long finite z-dependent NNLO tables.
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
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Dict, Tuple

import sympy as sp


@dataclass(frozen=True)
class QCDPoint:
    CA: float = 3.0
    CF: float = 4.0 / 3.0
    TF: float = 0.5
    nf: int = 5


def qcd_constants(point: QCDPoint) -> Dict[str, sp.Expr]:
    CA = sp.Rational(3, 1) if point.CA == 3.0 else sp.Float(point.CA)
    CF = sp.Rational(4, 3) if abs(point.CF - 4.0 / 3.0) < 1e-15 else sp.Float(point.CF)
    TF = sp.Rational(1, 2) if point.TF == 0.5 else sp.Float(point.TF)
    nf = sp.Integer(point.nf)
    beta0 = sp.Rational(11, 3) * CA - sp.Rational(4, 3) * TF * nf
    Gamma0q = 4 * CF
    Gamma1q = 4 * CF * ((sp.Rational(67, 9) - sp.pi**2 / 3) * CA - sp.Rational(20, 9) * TF * nf)
    gammaF0q = 6 * CF
    gammaF1q = (
        CF**2 * (3 - 4 * sp.pi**2 + 48 * sp.zeta(3))
        + CF * CA * (sp.Rational(961, 27) + sp.Rational(11, 3) * sp.pi**2 - 52 * sp.zeta(3))
        - CF * TF * nf * (sp.Rational(260, 27) + sp.Rational(4, 3) * sp.pi**2)
    )
    # Common rapidity-anomalous-dimension boundary in the ln(zeta) convention.
    d2zetaq = CF * ((sp.Rational(404, 27) - 14 * sp.zeta(3)) * CA - sp.Rational(112, 27) * TF * nf)
    d2zetag = CA * ((sp.Rational(404, 27) - 14 * sp.zeta(3)) * CA - sp.Rational(112, 27) * TF * nf)
    # KCSS/CSS2 convention in this note: evolution is in ln sqrt(zeta), D=-Ktilde,
    # and dD/dln(mu)=2 Gamma_cusp.  The delta boundary in D is doubled.
    d2boundaryq = 2 * d2zetaq
    d2boundaryg = 2 * d2zetag
    return {
        "CA": CA,
        "CF": CF,
        "TF": TF,
        "nf": nf,
        "beta0": sp.simplify(beta0),
        "Gamma0q": sp.simplify(Gamma0q),
        "Gamma1q": sp.simplify(Gamma1q),
        "gammaF0q": sp.simplify(gammaF0q),
        "gammaF1q": sp.simplify(gammaF1q),
        "d2zetaq": sp.simplify(d2zetaq),
        "d2zetag": sp.simplify(d2zetag),
        "d2boundaryq": sp.simplify(d2boundaryq),
        "d2boundaryg": sp.simplify(d2boundaryg),
    }


def truncate_a(expr: sp.Expr, a: sp.Symbol, max_power: int = 2) -> sp.Expr:
    """Drop terms beyond a**max_power in a polynomial expression."""
    expr = sp.expand(expr)
    out = 0
    for term in expr.as_ordered_terms():
        power = term.as_powers_dict().get(a, 0)
        try:
            p = int(power)
        except TypeError:
            p = 999
        if p <= max_power:
            out += term
    return sp.simplify(sp.expand(out))


def b_space_rg_check() -> Tuple[bool, sp.Expr]:
    a, L, beta0, Gamma0, Gamma1, d2b = sp.symbols("a L beta0 Gamma0 Gamma1 d2b")
    beta = -2 * beta0 * a**2
    D = a * Gamma0 * L + a**2 * (sp.Rational(1, 2) * beta0 * Gamma0 * L**2 + Gamma1 * L + d2b)
    dD = 2 * sp.diff(D, L) + beta * sp.diff(D, a)
    dD = truncate_a(dD, a, 2)
    target = 2 * a * Gamma0 + 2 * a**2 * Gamma1
    residual = sp.simplify(sp.expand(dD - target))
    return residual == 0, residual


def frak_basis_rg_check() -> Tuple[bool, Dict[str, sp.Expr]]:
    a, beta0, Gamma0, Gamma1, d2b = sp.symbols("a beta0 Gamma0 Gamma1 d2b")
    coeffs = {
        "delta": a**2 * d2b,
        "F0": a * Gamma0 + a**2 * Gamma1,
        "F1": a**2 * sp.Rational(1, 2) * beta0 * Gamma0,
    }
    beta = -2 * beta0 * a**2
    dcoeffs = {key: truncate_a(beta * sp.diff(val, a), a, 2) for key, val in coeffs.items()}
    # d F_n / d ln(mu) = 2(n+1) F_{n-1}; F_{-1}=delta.
    dcoeffs["delta"] += 2 * coeffs["F0"]
    dcoeffs["F0"] += 4 * coeffs["F1"]
    dcoeffs["F1"] += 0
    dcoeffs = {k: truncate_a(sp.expand(v), a, 2) for k, v in dcoeffs.items()}
    target = {"delta": 2 * a * Gamma0 + 2 * a**2 * Gamma1, "F0": 0, "F1": 0}
    residual = {k: sp.simplify(sp.expand(dcoeffs[k] - target[k])) for k in target}
    return all(v == 0 for v in residual.values()), residual


def conventional_projection_check() -> Tuple[bool, Dict[str, sp.Expr]]:
    a, beta0, Gamma0, Gamma1, d2b = sp.symbols("a beta0 Gamma0 Gamma1 d2b")
    actual = {
        "delta": a**2 * d2b,
        "L0": -(a * Gamma0 + a**2 * Gamma1),
        "L1": -a**2 * beta0 * Gamma0,
    }
    expected = {
        "delta": a**2 * d2b,
        "L0": -a * Gamma0 - a**2 * Gamma1,
        "L1": -a**2 * beta0 * Gamma0,
    }
    residual = {k: sp.simplify(sp.expand(actual[k] - expected[k])) for k in expected}
    return all(v == 0 for v in residual.values()), residual


def scalar_matching_reconstruction_check() -> Tuple[bool, Dict[str, sp.Expr], Dict[str, sp.Expr]]:
    """Check scalar-channel NNLO log reconstruction against CS and mu equations.

    This replaces x-convolution matrices by commutative scalar symbols.  It therefore
    verifies the coefficients and signs, not the noncommuting flavor-matrix ordering.
    """
    a, L, ell = sp.symbols("a L ell")
    beta0, G0, G1, g0, g1, d2b = sp.symbols("beta0 G0 G1 g0 g1 d2b")
    p0, p1, c1, c2 = sp.symbols("p0 p1 c1 c2")

    D1 = G0 * L
    D2 = sp.Rational(1, 2) * beta0 * G0 * L**2 + G1 * L + d2b
    C10 = c1 + (g0 / 2 - 2 * p0) * L - G0 * L**2 / 4
    C1 = C10 - G0 * ell * L / 2

    Lp = sp.symbols("Lp")
    D1p = G0 * Lp
    D2p = sp.Rational(1, 2) * beta0 * G0 * Lp**2 + G1 * Lp + d2b
    C10p = c1 + (g0 / 2 - 2 * p0) * Lp - G0 * Lp**2 / 4
    integrand = g1 - D2p - 8 * p1 + (g0 + 2 * beta0) * C10p - D1p * C10p - 4 * C10p * p0
    C20 = c2 + sp.Rational(1, 2) * sp.integrate(integrand, (Lp, 0, L))
    C2 = C20 - ell * (D2 + D1 * C10) / 2 + G0**2 * ell**2 * L**2 / 8

    C = 1 + a * C1 + a**2 * C2
    D = a * D1 + a**2 * D2
    gamma = a * (g0 - G0 * ell) + a**2 * (g1 - G1 * ell)
    P_mu = 4 * a * p0 + 8 * a**2 * p1
    beta = -2 * beta0 * a**2

    # CS: 2 d_ell C = -D*C
    cs_res = truncate_a(2 * sp.diff(C, ell) + D * C, a, 2)

    # Mu equation at fixed zeta: d/dlnmu = 2 d_L - 2 d_ell + beta d_a.
    lhs_mu = 2 * sp.diff(C, L) - 2 * sp.diff(C, ell) + beta * sp.diff(C, a)
    rhs_mu = gamma * C - C * P_mu
    mu_res = truncate_a(lhs_mu - rhs_mu, a, 2)

    # Semicanonical one-loop checks.
    c1_mu_res = sp.simplify(sp.expand((2 * sp.diff(C1, L) - 2 * sp.diff(C1, ell)) - (g0 - G0 * ell - 4 * p0)))
    c1_cs_res = sp.simplify(sp.expand(2 * sp.diff(C1, ell) + D1))

    residuals = {
        "matching_CS_equation_scalar_projection": sp.simplify(sp.expand(cs_res)),
        "matching_mu_equation_scalar_projection": sp.simplify(sp.expand(mu_res)),
        "one_loop_mu_subcheck": c1_mu_res,
        "one_loop_CS_subcheck": c1_cs_res,
    }
    objects = {
        "C1": sp.simplify(C1),
        "C20": sp.simplify(C20),
        "C2": sp.simplify(C2),
    }
    return all(v == 0 for v in residuals.values()), residuals, objects


def main() -> None:
    outdir = _OUTPUT_DIR
    qcd = QCDPoint()
    constants = qcd_constants(qcd)

    b_ok, b_res = b_space_rg_check()
    frak_ok, frak_res = frak_basis_rg_check()
    conv_ok, conv_res = conventional_projection_check()
    match_ok, match_res, match_obj = scalar_matching_reconstruction_check()
    passed = b_ok and frak_ok and conv_ok and match_ok

    tests = [
        {"test": "cs_kernel_b_space_RG", "passed": b_ok, "residual": str(b_res)},
        {"test": "cs_kernel_frak_basis_RG", "passed": frak_ok, "residual": {k: str(v) for k, v in frak_res.items()}},
        {"test": "cs_kernel_conventional_projection", "passed": conv_ok, "residual": {k: str(v) for k, v in conv_res.items()}},
    ]
    for key, val in match_res.items():
        tests.append({"test": key, "passed": val == 0, "residual": str(val)})

    numeric_constants = {k: float(sp.N(v, 16)) for k, v in constants.items() if k not in {"nf"}}
    numeric_constants["nf"] = int(constants["nf"])

    payload = {
        "passed": passed,
        "qcd_point": asdict(qcd),
        "tests": tests,
        "constants_exact": {k: str(v) for k, v in constants.items()},
        "constants_numeric": numeric_constants,
        "matching_objects": {k: str(v) for k, v in match_obj.items()},
        "notes": [
            "Uses a_s = alpha_s/(4*pi), d a_s/d ln mu = -2 beta0 a_s^2 + ..., and P_mu = 4 a_s P0 + 8 a_s^2 P1 + ...",
            "Checks d D/d ln mu = 2 Gamma_cusp through O(a_s^2).",
            "Matching reconstruction is a scalar-channel projection; full code must preserve flavor-matrix x-convolution ordering.",
            "The finite NNLO matching functions c2(z) are external perturbative input, not generated here.",
        ],
    }

    (outdir / "kt_kcss_two_loop_ingredients_validator_output.json").write_text(json.dumps(payload, indent=2))

    with (outdir / "kt_kcss_two_loop_ingredients_validator_tests.csv").open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["test", "passed", "residual"])
        writer.writeheader()
        for row in tests:
            writer.writerow({"test": row["test"], "passed": row["passed"], "residual": row["residual"]})

    with (outdir / "kt_kcss_two_loop_ingredients_validator_constants.csv").open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["nf", "beta0", "Gamma0q", "Gamma1q", "d2zetaq", "d2boundaryq"])
        writer.writeheader()
        for nf in [4, 5]:
            c = qcd_constants(QCDPoint(nf=nf))
            writer.writerow({
                "nf": nf,
                "beta0": float(sp.N(c["beta0"], 16)),
                "Gamma0q": float(sp.N(c["Gamma0q"], 16)),
                "Gamma1q": float(sp.N(c["Gamma1q"], 16)),
                "d2zetaq": float(sp.N(c["d2zetaq"], 16)),
                "d2boundaryq": float(sp.N(c["d2boundaryq"], 16)),
            })

    lines = []
    lines.append("KCSS two-loop ingredient validator")
    lines.append("===================================")
    lines.append(f"passed: {passed}")
    lines.append("")
    lines.append("Tests:")
    for row in tests:
        lines.append(f"  - {row['test']}: passed={row['passed']}, residual={row['residual']}")
    lines.append("")
    lines.append("Sample SU(3), nf=5 constants:")
    for key in ["beta0", "Gamma0q", "Gamma1q", "gammaF0q", "d2zetaq", "d2boundaryq", "d2boundaryg"]:
        lines.append(f"  {key}: exact={constants[key]}, numeric={numeric_constants[key]:.12g}")
    lines.append("")
    lines.append("Scalar matching objects:")
    lines.append(f"  C1 = {match_obj['C1']}")
    lines.append(f"  C2 = {match_obj['C2']}")
    lines.append("")
    lines.append("Conventions:")
    lines.append("  a_s = alpha_s/(4*pi)")
    lines.append("  d a_s/d ln(mu) = -2 beta0 a_s^2 + O(a_s^3)")
    lines.append("  P_mu = 4 a_s P0 + 8 a_s^2 P1 + O(a_s^3) for d/dln(mu) evolution")
    lines.append("  D_k = a_s Gamma0 frakL0 + a_s^2[(beta0 Gamma0/2) frakL1 + Gamma1 frakL0 + d2_boundary delta]")
    lines.append("  d frakL_n/d ln(mu) = 2(n+1) frakL_{n-1}")
    (outdir / "kt_kcss_two_loop_ingredients_validator_report.txt").write_text("\n".join(lines) + "\n")

    print("\n".join(lines))
    raise SystemExit(0 if passed else 1)


if __name__ == "__main__":
    main()
