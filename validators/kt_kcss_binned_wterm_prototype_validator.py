#!/usr/bin/env python3
"""
Binned W-term prototype validator for the K_CSS2 pure-k_T formalism.

This is a formal regression validator, not a phenomenology calculation.  It checks that
an evolved perturbative boundary at a profile scale mu_T reproduces the canonical
one-loop CSS2 singular W term at mu=Q after including the one-loop DGLAP conversion of
collinear PDFs.  It also checks the two-loop rapidity-kernel projection used by the
NLL/NNLL-style W-term prototype.
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
import math
from pathlib import Path
from typing import Dict, Iterable, Tuple

import sympy as sp


def trunc_a(expr: sp.Expr, a: sp.Symbol, order: int) -> sp.Expr:
    """Truncate a symbolic series in a through a**(order-1)."""
    return sp.series(sp.expand(expr), a, 0, order).removeO().expand()


def bin_L(n: int, qa: float, qb: float, Q: float) -> float:
    """Bin action of L_n(q,Q) over qa < |q| < qb."""
    if qb <= qa:
        raise ValueError("Need qb > qa")
    if qa == 0.0:
        Lb = math.log(Q * Q / (qb * qb))
        return -(Lb ** (n + 1)) / float(n + 1)
    La = math.log(Q * Q / (qa * qa))
    Lb = math.log(Q * Q / (qb * qb))
    return (La ** (n + 1) - Lb ** (n + 1)) / float(n + 1)


def bin_delta(qa: float, qb: float) -> float:
    """Bin action of delta^(2)(q). The first bin is taken closed at the origin."""
    return 1.0 if qa == 0.0 and qb > 0.0 else 0.0


def poly_to_distribution_coeffs(poly: sp.Expr, X: sp.Symbol) -> Dict[str, sp.Expr]:
    """Map c2 X^2 + c1 X + c0 to A*L1 + B*L0 + C*delta.

    The Fourier dictionary is F[L0] = -X and F[L1] = -X^2/2.
    """
    P = sp.Poly(sp.expand(poly), X)
    c2 = P.coeff_monomial(X ** 2)
    c1 = P.coeff_monomial(X)
    c0 = P.coeff_monomial(1)
    return {
        "L1": sp.simplify(-2 * c2),
        "L0": sp.simplify(-c1),
        "delta": sp.simplify(c0),
    }


def bin_action(coeffs: Dict[str, float], qa: float, qb: float, Q: float) -> float:
    return (
        coeffs.get("L1", 0.0) * bin_L(1, qa, qb, Q)
        + coeffs.get("L0", 0.0) * bin_L(0, qa, qb, Q)
        + coeffs.get("delta", 0.0) * bin_delta(qa, qb)
    )


def main() -> None:
    out_dir = _OUTPUT_DIR
    report_path = out_dir / "kt_kcss_binned_wterm_prototype_validator_report.txt"
    csv_path = out_dir / "kt_kcss_binned_wterm_prototype_validator_tests.csv"
    json_path = out_dir / "kt_kcss_binned_wterm_prototype_validator_output.json"
    bins_path = out_dir / "kt_kcss_binned_wterm_prototype_validator_bins.csv"

    # Symbols.  X = L_b(Q), L = ln(Q/mu_T), and Li = L_b(mu_T) = X - 2L.
    a, X, L = sp.symbols("a X L")
    CF, PA, PB, cA, cB, h = sp.symbols("CF PA PB cA cB h")
    G0, g0 = sp.symbols("G0 g0")
    rho, G1, beta0, d2 = sp.symbols("rho G1 beta0 d2")

    Li = X - 2 * L

    def leg_profile(P, c):
        # Boundary matching at (mu_T,zeta_T=mu_T^2), expressed in b space.
        # C^(1) = c + 3 C_F L_i - C_F L_i^2 - 2 P L_i.
        boundary = c + 3 * CF * Li - CF * Li ** 2 - 2 * P * Li
        # Express the input collinear PDF at mu_T through the PDF at Q.
        # df/dln(mu) = 4 a_s P^(0) \otimes f, hence f(mu_T) = f(Q) - 4 a_s L P \otimes f(Q).
        pdf_to_Q = -4 * P * L
        # One-leg evolution from (mu_T,mu_T^2) to (Q,Q^2) at fixed coupling through O(a).
        evolution = g0 * L + G0 * L ** 2 - G0 * L * X
        return sp.expand(boundary + pdf_to_Q + evolution)

    leg_A_profile = leg_profile(PA, cA)
    leg_B_profile = leg_profile(PB, cB)
    subs_quark = {G0: 4 * CF, g0: 6 * CF}

    leg_A_profile_q = sp.simplify(leg_A_profile.subs(subs_quark))
    leg_A_canonical = cA + 3 * CF * X - CF * X ** 2 - 2 * PA * X
    delta_leg = sp.simplify(sp.expand(leg_A_profile_q - leg_A_canonical))

    pair_profile = sp.simplify((h + leg_A_profile + leg_B_profile).subs(subs_quark))
    pair_canonical = h + cA + cB + 6 * CF * X - 2 * CF * X ** 2 - 2 * (PA + PB) * X
    delta_pair = sp.simplify(sp.expand(pair_profile - pair_canonical))

    # Distribution-basis coefficients of the canonical one-loop singular W term.
    dist_coeffs = poly_to_distribution_coeffs(pair_canonical, X)
    expected_coeffs = {
        "L1": 4 * CF,
        "L0": 2 * (PA + PB) - 6 * CF,
        "delta": h + cA + cB,
    }
    delta_coeffs = {k: sp.simplify(dist_coeffs[k] - expected_coeffs[k]) for k in expected_coeffs}

    # Two-loop rapidity kernel projection for the symmetric pair leg.
    # rho here denotes rho_pair = rho_A + rho_B.
    rapidity_L_basis = {
        "delta": 1 - a ** 2 * rho * d2,
        "L0": a * rho * G0 + a ** 2 * rho * G1,
        "L1": a ** 2 * (rho * beta0 * G0 - rho ** 2 * G0 ** 2),
    }
    rapidity_poly_from_L = trunc_a(
        rapidity_L_basis["delta"]
        + rapidity_L_basis["L0"] * (-X)
        + rapidity_L_basis["L1"] * (-sp.Rational(1, 2) * X ** 2),
        a,
        3,
    )
    D_b = a * G0 * X + a ** 2 * (sp.Rational(1, 2) * beta0 * G0 * X ** 2 + G1 * X + d2)
    rapidity_poly_expected = trunc_a(sp.exp(-rho * D_b), a, 3)
    delta_rapidity_projection = sp.simplify(trunc_a(rapidity_poly_from_L - rapidity_poly_expected, a, 3))

    # Numerical bin checks for several profile scales.  These are smoke tests for the
    # binned distribution action after the algebraic identity above is established.
    values = {
        CF: 4.0 / 3.0,
        PA: 0.37,
        PB: 0.22,
        cA: -0.31,
        cB: 0.14,
        h: 1.7,
    }
    Q = 8.0
    profile_L_values = [0.0, math.log(Q / 2.0), math.log(Q / 3.5)]
    edges = [0.0, 0.25, 0.5, 1.0, 2.0, 4.0, 8.0]

    canonical_coeffs_num = {k: float(sp.N(v.subs(values))) for k, v in expected_coeffs.items()}
    max_bin_delta = 0.0
    bin_rows = []
    for Lval in profile_L_values:
        profile_poly_num = sp.expand(pair_profile.subs(values).subs({L: Lval}))
        coeffs_profile_sym = poly_to_distribution_coeffs(profile_poly_num, X)
        coeffs_profile_num = {k: float(sp.N(v)) for k, v in coeffs_profile_sym.items()}
        for qa, qb in zip(edges[:-1], edges[1:]):
            w_profile = bin_action(coeffs_profile_num, qa, qb, Q)
            w_canonical = bin_action(canonical_coeffs_num, qa, qb, Q)
            diff = w_profile - w_canonical
            max_bin_delta = max(max_bin_delta, abs(diff))
            bin_rows.append(
                {
                    "L_profile": Lval,
                    "qa": qa,
                    "qb": qb,
                    "W_profile_coeff_bin": w_profile,
                    "W_canonical_coeff_bin": w_canonical,
                    "difference": diff,
                }
            )

    tests = [
        ("leg_profile_to_canonical_Oa", delta_leg),
        ("pair_profile_to_canonical_Oa", delta_pair),
        ("dist_coeff_L1", delta_coeffs["L1"]),
        ("dist_coeff_L0", delta_coeffs["L0"]),
        ("dist_coeff_delta", delta_coeffs["delta"]),
        ("two_loop_pair_rapidity_projection", delta_rapidity_projection),
    ]
    passed_symbolic = all(sp.simplify(expr) == 0 for _, expr in tests)
    passed_numeric = max_bin_delta < 5e-13
    passed = bool(passed_symbolic and passed_numeric)

    with csv_path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["test", "residual"])
        writer.writeheader()
        for name, expr in tests:
            writer.writerow({"test": name, "residual": str(sp.simplify(expr))})
        writer.writerow({"test": "max_bin_delta", "residual": f"{max_bin_delta:.17e}"})

    with bins_path.open("w", newline="") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=["L_profile", "qa", "qb", "W_profile_coeff_bin", "W_canonical_coeff_bin", "difference"],
        )
        writer.writeheader()
        writer.writerows(bin_rows)

    payload = {
        "passed": passed,
        "passed_symbolic": bool(passed_symbolic),
        "passed_numeric": bool(passed_numeric),
        "max_bin_delta": max_bin_delta,
        "canonical_distribution_coefficients": {k: str(v) for k, v in expected_coeffs.items()},
        "leg_profile_minus_canonical": str(delta_leg),
        "pair_profile_minus_canonical": str(delta_pair),
        "two_loop_pair_rapidity_projection_residual": str(delta_rapidity_projection),
        "Q_smoke_GeV": Q,
        "profile_L_values": profile_L_values,
        "bin_edges_GeV": edges,
    }
    json_path.write_text(json.dumps(payload, indent=2))

    report = f"""K_CSS2 binned W-term prototype validator
==============================================
passed: {passed}
passed_symbolic: {passed_symbolic}
passed_numeric: {passed_numeric}

Symbolic checks
---------------
leg profile -> canonical O(a): {sp.simplify(delta_leg)}
pair profile -> canonical O(a): {sp.simplify(delta_pair)}
distribution coefficient residuals:
  L1: {delta_coeffs['L1']}
  L0: {delta_coeffs['L0']}
  delta: {delta_coeffs['delta']}
two-loop pair rapidity projection residual: {delta_rapidity_projection}

Smoke-test settings
-------------------
Q = {Q} GeV
profile L values = {profile_L_values}
bin edges = {edges}
max bin-level profile-canonical difference = {max_bin_delta:.17e}

Interpretation
--------------
The O(a_s) fixed-order expansion of the profile-scale K_CSS2 W-term prototype is
independent of the auxiliary profile scale mu_T once one-loop DGLAP conversion of the
collinear PDFs is included.  The resulting bin-integrated distribution coefficients are
4 C_F L_1 + [2(P_A+P_B)-6 C_F] L_0 + (h+c_A+c_B) delta.
"""
    report_path.write_text(report)

    print(report)


if __name__ == "__main__":
    main()
