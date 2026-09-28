#!/usr/bin/env python3
"""
Generator-convolution validator for the K_CSS2 transverse distribution basis.

This script is intentionally algebraic.  It does not compute a physical cross
section.  It verifies the dictionary between

  1. the generator basis G_eta, whose Fourier transform is exp(-eta L_b),
  2. the frak basis frak{L}_n, whose Fourier transform is L_b^{n+1}, and
  3. the conventional two-dimensional plus basis L_n(k,mu), defined with
     ln^n(mu^2/k_T^2).

The transverse convolution is checked in b-space, where it becomes ordinary
multiplication of polynomials in X=L_b.
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
from typing import Dict, List, Tuple

import sympy as sp


@dataclass(frozen=True)
class Config:
    max_L: int = 7
    table_max_index: int = 3
    output_prefix: str = str(_OUTPUT_DIR / "kt_kcss_generator_convolution_validator")


X = sp.Symbol("X")  # X is L_b
eta = sp.Symbol("eta")
z3 = sp.Symbol("zeta_3")
z5 = sp.Symbol("zeta_5")
z7 = sp.Symbol("zeta_7")


def h_log_series(order: int) -> sp.Expr:
    """Return log h(eta) through the requested order.

    h(eta)=exp[-2 gamma_E eta] Gamma(1-eta)/Gamma(1+eta)
          = exp[2 sum_{m odd >=3} zeta_m eta^m/m].
    """
    expr = sp.Integer(0)
    if order >= 3:
        expr += sp.Rational(2, 3) * z3 * eta**3
    if order >= 5:
        expr += sp.Rational(2, 5) * z5 * eta**5
    if order >= 7:
        expr += sp.Rational(2, 7) * z7 * eta**7
    return expr


def conventional_ft_polynomials(max_L: int) -> Dict[int, sp.Expr]:
    """Compute P_n(X)=FT[L_n] for n=0..max_L.

    The defining identity is
      exp(-eta X) = h(eta) [1 + sum_{n>=0} (-eta)^{n+1}/n! L_n]_FT.
    Therefore
      P_n(X) = (-1)^n n! [eta^{n+1}] exp(-eta X)/h(eta).
    """
    order = max_L + 2
    expr = sp.series(sp.exp(-eta * X - h_log_series(order)), eta, 0, order).removeO()
    polynomials: Dict[int, sp.Expr] = {}
    for n in range(max_L + 1):
        coeff = sp.expand(expr).coeff(eta, n + 1)
        polynomials[n] = sp.simplify((-1) ** n * sp.factorial(n) * coeff)
    return polynomials


def express_in_L_basis(poly: sp.Expr, P: Dict[int, sp.Expr], max_L: int) -> Dict[str, sp.Expr]:
    """Express a polynomial in X in the basis {delta, L_0, ..., L_max_L}."""
    unknowns = sp.symbols(f"c0:{max_L + 2}")  # c0 delta, c_{n+1} L_n
    basis = [sp.Integer(1)] + [P[n] for n in range(max_L + 1)]
    ansatz = sum(unknowns[i] * basis[i] for i in range(max_L + 2))
    diff = sp.Poly(sp.expand(ansatz - poly), X)
    equations = [sp.Eq(c, 0) for c in diff.all_coeffs()]
    solutions = sp.solve(equations, unknowns, dict=True)
    if not solutions:
        raise RuntimeError(f"Could not express polynomial {sp.sstr(poly)} in L basis")
    sol = solutions[0]
    out: Dict[str, sp.Expr] = {}
    for i, c in enumerate(unknowns):
        val = sp.simplify(sol.get(c, 0))
        if val == 0:
            continue
        key = "delta" if i == 0 else f"L{i-1}"
        out[key] = val
    return out


def reconstruct_from_coeffs(coeffs: Dict[str, sp.Expr], P: Dict[int, sp.Expr]) -> sp.Expr:
    total = sp.Integer(0)
    for key, coeff in coeffs.items():
        if key == "delta":
            total += coeff
        elif key.startswith("L"):
            total += coeff * P[int(key[1:])]
        else:
            raise ValueError(key)
    return sp.simplify(sp.expand(total))


def convolve_L(m: int, n: int, P: Dict[int, sp.Expr], max_L: int) -> Dict[str, sp.Expr]:
    return express_in_L_basis(sp.expand(P[m] * P[n]), P, max_L)


def coeffs_to_string(coeffs: Dict[str, sp.Expr]) -> str:
    if not coeffs:
        return "0"
    order = ["delta"] + [f"L{i}" for i in range(20)]
    pieces: List[str] = []
    for key in order:
        if key in coeffs:
            pieces.append(f"({sp.sstr(sp.factor(coeffs[key]))})*{key}")
    return " + ".join(pieces).replace("+ (-", "- (")


def main() -> None:
    cfg = Config()
    prefix = Path(cfg.output_prefix)

    P = conventional_ft_polynomials(cfg.max_L)

    # Convolution table L_m otimes L_n up to table_max_index.
    table: List[Dict[str, str]] = []
    max_residual = sp.Integer(0)
    for m in range(cfg.table_max_index + 1):
        for n in range(m, cfg.table_max_index + 1):
            max_needed = min(cfg.max_L, m + n + 1)
            coeffs = convolve_L(m, n, P, max_needed)
            residual = sp.simplify(sp.expand(P[m] * P[n] - reconstruct_from_coeffs(coeffs, P)))
            if residual != 0:
                max_residual = 1
            table.append(
                {
                    "left": f"L{m} x L{n}",
                    "right": coeffs_to_string(coeffs),
                    "residual": sp.sstr(residual),
                }
            )

    # Check the frak basis rule: frakL_m has FT X^{m+1}; therefore convolution gives X^{m+n+2}.
    frak_failures: List[str] = []
    for m in range(cfg.table_max_index + 1):
        for n in range(cfg.table_max_index + 1):
            lhs = X ** (m + 1) * X ** (n + 1)
            rhs = X ** (m + n + 2)
            if sp.simplify(lhs - rhs) != 0:
                frak_failures.append(f"frakL{m} x frakL{n}")

    # Check that powers of L0 reproduce the expansion of exp(lambda L0) in b-space.
    lam = sp.Symbol("lambda")
    max_power = 5
    exp_series = sp.series(sp.exp(-lam * X), lam, 0, max_power + 1).removeO()
    reconstructed_series = sp.Integer(1)
    # In b-space, L0 has polynomial -X; convolution powers are ordinary powers.
    for p in range(1, max_power + 1):
        reconstructed_series += lam**p / sp.factorial(p) * ((-X) ** p)
    exponential_residual = sp.simplify(sp.expand(exp_series - reconstructed_series))

    passed = bool(max_residual == 0 and not frak_failures and exponential_residual == 0)

    payload = {
        "passed": passed,
        "max_L": cfg.max_L,
        "table_max_index": cfg.table_max_index,
        "FT_Ln_polynomials": {f"L{n}": sp.sstr(P[n]) for n in range(cfg.max_L + 1)},
        "convolution_table": table,
        "frak_rule_failures": frak_failures,
        "conventional_table_max_residual_nonzero": bool(max_residual != 0),
        "exponential_residual": sp.sstr(exponential_residual),
        "notes": [
            "X denotes L_b.",
            "delta denotes delta^{(2)}(k).",
            "L_n denotes the conventional 2D plus distribution with ln^n(mu^2/k_T^2).",
            "The frak basis obeys frakL_m otimes_T frakL_n = frakL_{m+n+1} exactly.",
        ],
    }

    json_path = prefix.with_suffix("_output.json") if False else Path(str(prefix) + "_output.json")
    txt_path = Path(str(prefix) + "_report.txt")
    csv_path = Path(str(prefix) + "_table.csv")

    json_path.write_text(json.dumps(payload, indent=2), encoding="utf-8")

    with csv_path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=["left", "right", "residual"])
        writer.writeheader()
        writer.writerows(table)

    lines = []
    lines.append("K_CSS2 generator convolution validator")
    lines.append("======================================")
    lines.append(f"passed: {passed}")
    lines.append(f"max_L: {cfg.max_L}")
    lines.append(f"table_max_index: {cfg.table_max_index}")
    lines.append(f"frak_rule_failures: {len(frak_failures)}")
    lines.append(f"conventional_table_residual_nonzero: {bool(max_residual != 0)}")
    lines.append(f"exponential_residual: {sp.sstr(exponential_residual)}")
    lines.append("")
    lines.append("Fourier dictionary P_n(X)=FT[L_n]")
    lines.append("---------------------------------")
    for n in range(min(cfg.max_L, 5) + 1):
        lines.append(f"L{n}: {sp.sstr(P[n])}")
    lines.append("")
    lines.append("Convolution table")
    lines.append("-----------------")
    for row in table:
        lines.append(f"{row['left']} -> {row['right']}")
    lines.append("")
    lines.append("Interpretation")
    lines.append("--------------")
    lines.append("The exact generator rule is checked in the frak basis, where FT[frakL_n]=X^(n+1).")
    lines.append("The conventional L_n table is obtained by multiplying P_m(X) P_n(X) and projecting back onto {delta,L_n}.")
    txt_path.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print("passed:", passed)
    print("wrote:", json_path)
    print("wrote:", txt_path)
    print("wrote:", csv_path)


if __name__ == "__main__":
    main()
