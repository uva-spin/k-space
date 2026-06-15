#!/usr/bin/env python3
"""
One-loop K_CSS2 binned singular-spectrum validator.

This script implements the bin functionals

    B_n(q_a,q_b;Q) = int_{q_a^2}^{q_b^2} dq_T^2 pi L_n(q,Q)

for the two-dimensional plus distributions used in the project note, and checks the
one-loop diagonal Drell-Yan singular kernel normalization.

It is intentionally minimal: no PDF interface, no Y term, no resummation, and no
nonperturbative model.  The purpose is to catch signs, pi factors, delta support, and
plus-prescription errors before moving to higher-order convolution algebra.
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
import json
import math
from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Iterable, List

CF_DEFAULT = 4.0 / 3.0


def a_s(alpha_s: float) -> float:
    """Project convention: a_s = alpha_s/(4*pi)."""
    return alpha_s / (4.0 * math.pi)


def h_dy_quark_1(cf: float = CF_DEFAULT) -> float:
    """Reduced CSS2 one-loop Drell-Yan hard factor coefficient at mu=Q."""
    return cf * (-16.0 + 7.0 * math.pi**2 / 3.0)


def log_Q2_over_q2(Q: float, q: float) -> float:
    if Q <= 0.0:
        raise ValueError("Q must be positive")
    if q <= 0.0:
        raise ValueError("q must be positive when evaluating an ordinary logarithm")
    return math.log((Q * Q) / (q * q))


def B_n(n: int, qa: float, qb: float, Q: float, *, zero_tol: float = 1e-15) -> float:
    """Bin integral of pi*L_n over q_a^2 <= q_T^2 <= q_b^2.

    For qa > 0 this is the ordinary integral.  For qa = 0 the plus prescription is
    imposed with reference scale Q, so B_n(0,Q;Q)=0.
    """
    if n < 0:
        raise ValueError("n must be nonnegative")
    if Q <= 0.0:
        raise ValueError("Q must be positive")
    if qa < -zero_tol or qb <= 0.0 or qb <= max(qa, 0.0):
        raise ValueError("Require 0 <= qa < qb and qb > 0")

    if abs(qa) <= zero_tol:
        Lb = log_Q2_over_q2(Q, qb)
        return -(Lb ** (n + 1)) / float(n + 1)

    La = log_Q2_over_q2(Q, qa)
    Lb = log_Q2_over_q2(Q, qb)
    return (La ** (n + 1) - Lb ** (n + 1)) / float(n + 1)


def B_delta(qa: float, qb: float, *, zero_tol: float = 1e-15) -> float:
    """Bin integral of pi*delta^(2)(q).  The delta contributes only to the zero bin."""
    if qa < -zero_tol or qb <= 0.0 or qb <= max(qa, 0.0):
        raise ValueError("Require 0 <= qa < qb and qb > 0")
    return 1.0 if abs(qa) <= zero_tol else 0.0


def diagonal_singular_weight(
    qa: float,
    qb: float,
    Q: float,
    alpha_s: float,
    *,
    cf: float = CF_DEFAULT,
    include_delta_constant: bool = True,
) -> float:
    """One-loop diagonal q qbar singular bin weight in the canonical mu=Q check.

    This returns the coefficient multiplying e_q^2 f_q/A f_qbar/B for the one-loop
    singular piece only.  For bins away from the origin it equals

        alpha_s/pi * C_F * [B_1 - 3/2 B_0].

    For the zero bin, include_delta_constant adds the hard and canonical TMD delta
    constants h_q^(1) - (pi^2/3) C_F.
    """
    b0 = B_n(0, qa, qb, Q)
    b1 = B_n(1, qa, qb, Q)
    bd = B_delta(qa, qb)
    delta_coeff = h_dy_quark_1(cf) - (math.pi**2 / 3.0) * cf
    return a_s(alpha_s) * (
        4.0 * cf * b1
        - 6.0 * cf * b0
        + (delta_coeff if include_delta_constant else 0.0) * bd
    )


def ordinary_offzero_weight(qa: float, qb: float, Q: float, alpha_s: float, *, cf: float = CF_DEFAULT) -> float:
    """Ordinary-function integral for qa>0, in dq_T^2 normalization."""
    if qa <= 0.0:
        raise ValueError("ordinary_offzero_weight is only valid for qa > 0")
    return (alpha_s / math.pi) * cf * (B_n(1, qa, qb, Q) - 1.5 * B_n(0, qa, qb, Q))


@dataclass
class BinRow:
    qa: float
    qb: float
    B0: float
    B1: float
    Bdelta: float
    diag_singular_weight: float
    ordinary_offzero_weight: float | None
    offzero_difference: float | None


def make_rows(edges: Iterable[float], Q: float, alpha_s: float, cf: float = CF_DEFAULT) -> List[BinRow]:
    edge_list = list(edges)
    rows: List[BinRow] = []
    for qa, qb in zip(edge_list[:-1], edge_list[1:]):
        diag = diagonal_singular_weight(qa, qb, Q, alpha_s, cf=cf)
        if qa > 0.0:
            ordinary = ordinary_offzero_weight(qa, qb, Q, alpha_s, cf=cf)
            diff = diag - ordinary
        else:
            ordinary = None
            diff = None
        rows.append(
            BinRow(
                qa=qa,
                qb=qb,
                B0=B_n(0, qa, qb, Q),
                B1=B_n(1, qa, qb, Q),
                Bdelta=B_delta(qa, qb),
                diag_singular_weight=diag,
                ordinary_offzero_weight=ordinary,
                offzero_difference=diff,
            )
        )
    return rows


def run_unit_tests(Q: float, alpha_s: float, edges: List[float], tol: float = 5e-13) -> dict:
    tests = {}

    # Plus-prescription null cumulant: int_0^Q L_n = 0.
    plus_null = {f"B{n}(0,Q)": B_n(n, 0.0, Q, Q) for n in range(4)}
    tests["plus_null_max_abs"] = max(abs(v) for v in plus_null.values())
    tests["plus_null_values"] = plus_null

    # Additivity across bins, including the zero bin. Use the first three edges if possible.
    if len(edges) >= 4:
        q0, q1, q2, q3 = edges[:4]
        add_vals = {}
        for n in range(3):
            add_vals[f"B{n}(q0,q2)-B{n}(q0,q1)-B{n}(q1,q2)"] = (
                B_n(n, q0, q2, Q) - B_n(n, q0, q1, Q) - B_n(n, q1, q2, Q)
            )
            add_vals[f"B{n}(q1,q3)-B{n}(q1,q2)-B{n}(q2,q3)"] = (
                B_n(n, q1, q3, Q) - B_n(n, q1, q2, Q) - B_n(n, q2, q3, Q)
            )
        tests["additivity_max_abs"] = max(abs(v) for v in add_vals.values())
        tests["additivity_values"] = add_vals

    # Off-origin diagonal check.
    rows = make_rows(edges, Q, alpha_s)
    diffs = [abs(r.offzero_difference) for r in rows if r.offzero_difference is not None]
    tests["offzero_diag_max_abs_difference"] = max(diffs) if diffs else 0.0

    tests["pass"] = all(
        tests.get(key, 0.0) < tol
        for key in ["plus_null_max_abs", "additivity_max_abs", "offzero_diag_max_abs_difference"]
        if key in tests
    )
    tests["tolerance"] = tol
    return tests


def parse_edges(text: str) -> List[float]:
    edges = [float(x.strip()) for x in text.split(",") if x.strip()]
    if len(edges) < 2:
        raise ValueError("Need at least two bin edges")
    if edges[0] < 0.0:
        raise ValueError("First edge must be nonnegative")
    if any(b <= a for a, b in zip(edges[:-1], edges[1:])):
        raise ValueError("Edges must be strictly increasing")
    return edges


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--Q", type=float, default=8.0, help="Hard scale Q in GeV")
    parser.add_argument("--alpha-s", type=float, default=0.25, help="alpha_s(Q)")
    parser.add_argument("--edges", type=str, default="0,0.25,0.5,1,2,4,8", help="Comma-separated q_T bin edges in GeV")
    parser.add_argument("--out-prefix", type=str, default=str(_OUTPUT_DIR / "kt_kcss_one_loop_validator_report"), help="Output prefix for CSV/JSON files")
    args = parser.parse_args()

    edges = parse_edges(args.edges)
    rows = make_rows(edges, args.Q, args.alpha_s)
    tests = run_unit_tests(args.Q, args.alpha_s, edges)

    out_prefix = Path(args.out_prefix)
    csv_path = out_prefix.with_suffix(".csv")
    json_path = out_prefix.with_suffix(".json")
    txt_path = out_prefix.with_suffix(".txt")

    with csv_path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(asdict(rows[0]).keys()))
        writer.writeheader()
        for row in rows:
            writer.writerow(asdict(row))

    payload = {
        "Q": args.Q,
        "alpha_s": args.alpha_s,
        "a_s": a_s(args.alpha_s),
        "edges": edges,
        "tests": tests,
        "rows": [asdict(r) for r in rows],
    }
    json_path.write_text(json.dumps(payload, indent=2), encoding="utf-8")

    lines = []
    lines.append("K_CSS2 one-loop binned singular-spectrum validator")
    lines.append(f"Q = {args.Q:g} GeV, alpha_s = {args.alpha_s:g}, a_s = {a_s(args.alpha_s):.12g}")
    lines.append(f"Bin edges = {edges}")
    lines.append("")
    lines.append("Unit tests:")
    lines.append(f"  plus null max abs            = {tests['plus_null_max_abs']:.6e}")
    lines.append(f"  additivity max abs           = {tests.get('additivity_max_abs', 0.0):.6e}")
    lines.append(f"  offzero diagonal max abs diff = {tests['offzero_diag_max_abs_difference']:.6e}")
    lines.append(f"  pass                         = {tests['pass']}")
    lines.append("")
    lines.append("Rows:")
    lines.append("  qa       qb       B0              B1              Bdelta  diag_weight       offzero_diff")
    for r in rows:
        diff = "--" if r.offzero_difference is None else f"{r.offzero_difference:.3e}"
        lines.append(
            f"  {r.qa:<8g} {r.qb:<8g} {r.B0:<15.8g} {r.B1:<15.8g} {r.Bdelta:<7.1f} {r.diag_singular_weight:<17.8g} {diff}"
        )
    txt_path.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print("\n".join(lines))
    print(f"\nWrote {csv_path}")
    print(f"Wrote {json_path}")
    print(f"Wrote {txt_path}")


if __name__ == "__main__":
    main()
