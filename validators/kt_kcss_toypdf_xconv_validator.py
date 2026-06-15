#!/usr/bin/env python3
r"""
Flavor-vector toy-PDF x-convolution validator for the one-loop K_CSS2 W term.

This script implements the second executable validation layer in the living note:
  - the longitudinal x-space plus prescription for P_qq^(0),
  - endpoint cancellation in b_{q<-q}=2 P_qq^(0)-3 C_F delta(1-z),
  - full flavor-vector one-leg insertions in the binned one-loop singular W term,
  - duplicate master-luminosity and explicit-leg construction paths,
  - full-bin additivity.

It is not a phenomenology code.  The PDFs are analytic positive toy functions.
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
from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Callable, Dict, Iterable, Mapping, Tuple

from scipy.integrate import quad

CF = 4.0 / 3.0
TF = 0.5
PI = math.pi
Q = 8.0
ALPHA_S = 0.25
A_S = ALPHA_S / (4.0 * PI)
X_A = 0.32
X_B = 0.21
BIN_EDGES = [0.0, 0.25, 0.50, 1.00, 2.00, 4.00, 8.00]
FLAVORS = ("u", "d", "s")
CHARGES = {"u": 2.0 / 3.0, "d": -1.0 / 3.0, "s": -1.0 / 3.0}
HARD_H1 = CF * (-16.0 + 7.0 * PI * PI / 3.0)
X_TESTS = [0.05, 0.12, 0.21, 0.32, 0.55, 0.78]


@dataclass(frozen=True)
class ToyPDF:
    norm: float
    a: float
    b: float
    c: float = 0.0

    def __call__(self, x: float) -> float:
        if not (0.0 < x < 1.0):
            if abs(x - 1.0) < 1e-15:
                return 0.0
            raise ValueError(f"x out of range: {x}")
        return self.norm * x**self.a * (1.0 - x) ** self.b * (1.0 + self.c * x)

    def derivative(self, x: float) -> float:
        val = self(x)
        return val * (self.a / x - self.b / (1.0 - x) + self.c / (1.0 + self.c * x))


PDF_A: Dict[str, ToyPDF] = {
    "u": ToyPDF(1.85, -0.25, 3.20, 0.30),
    "ubar": ToyPDF(0.27, -0.15, 7.00, 0.20),
    "d": ToyPDF(1.08, -0.20, 4.10, 0.10),
    "dbar": ToyPDF(0.33, -0.12, 6.50, 0.10),
    "s": ToyPDF(0.18, -0.10, 7.50, 0.00),
    "sbar": ToyPDF(0.18, -0.10, 7.50, 0.00),
    "g": ToyPDF(3.60, -0.35, 5.00, 0.40),
}
PDF_B: Dict[str, ToyPDF] = {
    "u": ToyPDF(1.30, -0.22, 3.60, 0.15),
    "ubar": ToyPDF(0.42, -0.18, 6.30, 0.10),
    "d": ToyPDF(1.34, -0.19, 3.85, 0.25),
    "dbar": ToyPDF(0.29, -0.10, 6.90, 0.15),
    "s": ToyPDF(0.21, -0.11, 7.30, 0.05),
    "sbar": ToyPDF(0.19, -0.11, 7.40, 0.00),
    "g": ToyPDF(4.05, -0.32, 5.35, 0.20),
}


def anti(q: str) -> str:
    return f"{q}bar"


def B_n(n: int, qa: float, qb: float, q_ref: float = Q) -> float:
    if qa < 0.0 or qb <= qa or qb > q_ref + 1e-14:
        raise ValueError(f"invalid bin: {qa}, {qb}")
    if qa == 0.0:
        if abs(qb - q_ref) < 1e-14:
            return 0.0
        Lb = math.log(q_ref * q_ref / (qb * qb))
        return -(Lb ** (n + 1)) / (n + 1)
    La = math.log(q_ref * q_ref / (qa * qa))
    Lb = math.log(q_ref * q_ref / (qb * qb))
    return (La ** (n + 1) - Lb ** (n + 1)) / (n + 1)


def B_delta(qa: float, qb: float) -> float:
    return 1.0 if qa == 0.0 and qb > 0.0 else 0.0


def R(kernel: Callable[[float], float], pdf: ToyPDF, x: float) -> float:
    def integrand(z: float) -> float:
        return kernel(z) * pdf(x / z) / z
    val, _ = quad(integrand, x, 1.0, epsabs=1e-12, epsrel=1e-11, limit=200)
    return val


def P_basic(pdf: ToyPDF, x: float) -> float:
    """P[f](x)= int_x^1 dz (phi(z)-f(x))/(1-z) + f(x) ln(1-x)."""
    fx = pdf(x)
    fpx = pdf.derivative(x)

    def integrand(z: float) -> float:
        if 1.0 - z < 1e-8:
            return fx + x * fpx
        return (pdf(x / z) / z - fx) / (1.0 - z)

    val, _ = quad(integrand, x, 1.0, epsabs=1e-12, epsrel=1e-11, limit=200)
    return val + fx * math.log(1.0 - x)


def Pqq_decomp(pdf: ToyPDF, x: float) -> float:
    return CF * (2.0 * P_basic(pdf, x) - R(lambda z: 1.0 + z, pdf, x) + 1.5 * pdf(x))


def Pqq_direct(pdf: ToyPDF, x: float) -> float:
    fx = pdf(x)
    fpx = pdf.derivative(x)

    def integrand(z: float) -> float:
        # ((1+z^2)/(1-z))*(phi-f) -> 2*(f+x f') as z -> 1.
        if 1.0 - z < 1e-8:
            return 2.0 * (fx + x * fpx)
        return ((1.0 + z * z) / (1.0 - z)) * (pdf(x / z) / z - fx)

    val, _ = quad(integrand, x, 1.0, epsabs=1e-12, epsrel=1e-11, limit=200)
    return CF * (val + fx * (2.0 * math.log(1.0 - x) + x + 0.5 * x * x))


def bqq_decomp(pdf: ToyPDF, x: float) -> float:
    return CF * (4.0 * P_basic(pdf, x) - 2.0 * R(lambda z: 1.0 + z, pdf, x))


def bqq_uncancelled(pdf: ToyPDF, x: float) -> float:
    return 2.0 * Pqq_decomp(pdf, x) - 3.0 * CF * pdf(x)


def bqg(pdf_g: ToyPDF, x: float) -> float:
    return 2.0 * TF * R(lambda z: z * z + (1.0 - z) ** 2, pdf_g, x)


def cqq(pdf: ToyPDF, x: float) -> float:
    return CF * (-(PI * PI) / 6.0 * pdf(x) + 2.0 * R(lambda z: 1.0 - z, pdf, x))


def cqg(pdf_g: ToyPDF, x: float) -> float:
    return R(lambda z: 2.0 * z * (1.0 - z), pdf_g, x)


def leg_b(pdfset: Mapping[str, ToyPDF], parton: str, x: float) -> float:
    return bqq_decomp(pdfset[parton], x) + bqg(pdfset["g"], x)


def leg_c(pdfset: Mapping[str, ToyPDF], parton: str, x: float) -> float:
    return cqq(pdfset[parton], x) + cqg(pdfset["g"], x)


@dataclass(frozen=True)
class Luminosity:
    phi0: float
    phi_b_A: float
    phi_b_B: float
    phi_c_A: float
    phi_c_B: float


def luminosity(q: str, pdfA: Mapping[str, ToyPDF], pdfB: Mapping[str, ToyPDF], xA: float, xB: float) -> Luminosity:
    aq = anti(q)
    fAq, fAaq = pdfA[q](xA), pdfA[aq](xA)
    fBq, fBaq = pdfB[q](xB), pdfB[aq](xB)
    return Luminosity(
        phi0=fAq * fBaq + fAaq * fBq,
        phi_b_A=leg_b(pdfA, q, xA) * fBaq + leg_b(pdfA, aq, xA) * fBq,
        phi_b_B=fAq * leg_b(pdfB, aq, xB) + fAaq * leg_b(pdfB, q, xB),
        phi_c_A=leg_c(pdfA, q, xA) * fBaq + leg_c(pdfA, aq, xA) * fBq,
        phi_c_B=fAq * leg_c(pdfB, aq, xB) + fAaq * leg_c(pdfB, q, xB),
    )


@dataclass(frozen=True)
class BinResult:
    qa: float
    qb: float
    B0: float
    B1: float
    Bdelta: float
    born: float
    l1: float
    l0: float
    delta_one_loop: float
    total: float


def W_master(qa: float, qb: float, pdfA: Mapping[str, ToyPDF]=PDF_A, pdfB: Mapping[str, ToyPDF]=PDF_B, xA: float=X_A, xB: float=X_B) -> BinResult:
    B0, B1, Bd = B_n(0, qa, qb), B_n(1, qa, qb), B_delta(qa, qb)
    born = l1 = l0 = delta1 = 0.0
    for q in FLAVORS:
        e2 = CHARGES[q] ** 2
        lum = luminosity(q, pdfA, pdfB, xA, xB)
        born += e2 * lum.phi0 * Bd
        l1 += e2 * A_S * 4.0 * CF * lum.phi0 * B1
        l0 += e2 * A_S * (lum.phi_b_A + lum.phi_b_B) * B0
        delta1 += e2 * A_S * (HARD_H1 * lum.phi0 + lum.phi_c_A + lum.phi_c_B) * Bd
    return BinResult(qa, qb, B0, B1, Bd, born, l1, l0, delta1, born + l1 + l0 + delta1)


def channel_terms(pdfA: Mapping[str, ToyPDF], pdfB: Mapping[str, ToyPDF], xA: float, xB: float, pA: str, pB: str) -> Tuple[float, float, float]:
    # Returns (phi0, b-insertion sum, c-insertion sum) for one explicit channel pA * pB.
    phi0 = pdfA[pA](xA) * pdfB[pB](xB)
    bsum = leg_b(pdfA, pA, xA) * pdfB[pB](xB) + pdfA[pA](xA) * leg_b(pdfB, pB, xB)
    csum = leg_c(pdfA, pA, xA) * pdfB[pB](xB) + pdfA[pA](xA) * leg_c(pdfB, pB, xB)
    return phi0, bsum, csum


def W_explicit(qa: float, qb: float, pdfA: Mapping[str, ToyPDF]=PDF_A, pdfB: Mapping[str, ToyPDF]=PDF_B, xA: float=X_A, xB: float=X_B) -> BinResult:
    B0, B1, Bd = B_n(0, qa, qb), B_n(1, qa, qb), B_delta(qa, qb)
    born = l1 = l0 = delta1 = 0.0
    for q in FLAVORS:
        e2 = CHARGES[q] ** 2
        # q_A * qbar_B and qbar_A * q_B
        p0a, ba, ca = channel_terms(pdfA, pdfB, xA, xB, q, anti(q))
        p0b, bb, cb = channel_terms(pdfA, pdfB, xA, xB, anti(q), q)
        phi0 = p0a + p0b
        bsum = ba + bb
        csum = ca + cb
        born += e2 * phi0 * Bd
        l1 += e2 * A_S * 4.0 * CF * phi0 * B1
        l0 += e2 * A_S * bsum * B0
        delta1 += e2 * A_S * (HARD_H1 * phi0 + csum) * Bd
    return BinResult(qa, qb, B0, B1, Bd, born, l1, l0, delta1, born + l1 + l0 + delta1)


def run_tests() -> Dict[str, float | bool]:
    pdfsets = [PDF_A, PDF_B]
    partons = ["u", "ubar", "d", "dbar", "s", "sbar"]

    deltas_P = []
    deltas_b = []
    for pdfset in pdfsets:
        for parton in partons:
            pdf = pdfset[parton]
            for x in X_TESTS:
                deltas_P.append(abs(Pqq_direct(pdf, x) - Pqq_decomp(pdf, x)))
                deltas_b.append(abs(bqq_uncancelled(pdf, x) - bqq_decomp(pdf, x)))
    delta_P = max(deltas_P)
    delta_b = max(deltas_b)

    delta_W = 0.0
    for i in range(len(BIN_EDGES)-1):
        qa, qb = BIN_EDGES[i], BIN_EDGES[i+1]
        delta_W = max(delta_W, abs(W_master(qa, qb).total - W_explicit(qa, qb).total))

    delta_add = 0.0
    # adjacent-bin additivity checks
    for i in range(len(BIN_EDGES)-2):
        qa, qb, qc = BIN_EDGES[i], BIN_EDGES[i+1], BIN_EDGES[i+2]
        lhs = W_master(qa, qc).total
        rhs = W_master(qa, qb).total + W_master(qb, qc).total
        delta_add = max(delta_add, abs(lhs-rhs))

    # Full cumulant additivity over all bins
    grid_total = sum(W_master(BIN_EDGES[i], BIN_EDGES[i+1]).total for i in range(len(BIN_EDGES)-1))
    cumulant = W_master(0.0, Q).total
    delta_full = abs(grid_total - cumulant)

    delta_swap = 0.0
    for i in range(len(BIN_EDGES)-1):
        qa, qb = BIN_EDGES[i], BIN_EDGES[i+1]
        orig = W_master(qa, qb, PDF_A, PDF_B, X_A, X_B).total
        swap = W_master(qa, qb, PDF_B, PDF_A, X_B, X_A).total
        delta_swap = max(delta_swap, abs(orig-swap))

    plus_cumulant_null = abs(B_n(0, 0.0, Q)) + abs(B_n(1, 0.0, Q))

    passed = all(val < 1e-10 for val in [delta_P, delta_b, delta_W, delta_add, delta_full, delta_swap, plus_cumulant_null])
    return {
        "passed": passed,
        "Delta_P_direct_vs_decomp": delta_P,
        "Delta_b_endpoint_cancel": delta_b,
        "Delta_W_master_vs_explicit_leg": delta_W,
        "Delta_add_adjacent_bins": delta_add,
        "Delta_full_cumulant_bins": delta_full,
        "Delta_A_B_swap": delta_swap,
        "plus_cumulant_null": plus_cumulant_null,
    }


def write_outputs(outdir: Path = _OUTPUT_DIR) -> None:
    rows = [W_master(BIN_EDGES[i], BIN_EDGES[i+1]) for i in range(len(BIN_EDGES)-1)]
    tests = run_tests()
    params = {
        "PDF_A": {k: asdict(v) for k, v in PDF_A.items()},
        "PDF_B": {k: asdict(v) for k, v in PDF_B.items()},
    }
    lums = {q: asdict(luminosity(q, PDF_A, PDF_B, X_A, X_B)) for q in FLAVORS}

    payload = {
        "metadata": {
            "description": "flavor-vector toy-PDF x-convolution validator for one-loop K_CSS2",
            "Q_GeV": Q,
            "alpha_s": ALPHA_S,
            "a_s": A_S,
            "x_A": X_A,
            "x_B": X_B,
            "bin_edges_GeV": BIN_EDGES,
            "x_tests": X_TESTS,
            "hard_h1": HARD_H1,
            "flavors": FLAVORS,
        },
        "toy_pdf_parameters": params,
        "flavor_luminosities": lums,
        "tests": tests,
        "bins": [asdict(row) for row in rows],
        "cumulant_0_to_Q": asdict(W_master(0.0, Q)),
    }

    json_path = outdir / "kt_kcss_toypdf_xconv_validator_output.json"
    csv_path = outdir / "kt_kcss_toypdf_xconv_validator_report.csv"
    txt_path = outdir / "kt_kcss_toypdf_xconv_validator_report.txt"
    json_path.write_text(json.dumps(payload, indent=2, sort_keys=True))

    with csv_path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["qa", "qb", "B0", "B1", "Bdelta", "born", "l1", "l0", "delta_one_loop", "total"])
        writer.writeheader()
        for row in rows:
            writer.writerow(asdict(row))

    lines = []
    lines.append("Flavor-vector toy-PDF x-convolution validator")
    lines.append("================================================")
    lines.append(f"Q = {Q:.8g} GeV, alpha_s = {ALPHA_S:.8g}, a_s = {A_S:.12g}")
    lines.append(f"x_A = {X_A:.8g}, x_B = {X_B:.8g}")
    lines.append(f"bin edges [GeV] = {BIN_EDGES}")
    lines.append("")
    lines.append("Regression tests:")
    for key, value in tests.items():
        if isinstance(value, bool):
            lines.append(f"  {key}: {value}")
        else:
            lines.append(f"  {key}: {value:.16e}")
    lines.append("")
    lines.append("Bin totals:")
    for row in rows:
        lines.append(
            f"  [{row.qa:5.2f}, {row.qb:5.2f}]  "
            f"Born={row.born:+.12e}  L1={row.l1:+.12e}  L0={row.l0:+.12e}  "
            f"Delta1={row.delta_one_loop:+.12e}  Total={row.total:+.12e}"
        )
    c = W_master(0.0, Q)
    lines.append("")
    lines.append(
        "Cumulant [0,Q]: "
        f"Born={c.born:+.12e}  L1={c.l1:+.12e}  L0={c.l0:+.12e}  "
        f"Delta1={c.delta_one_loop:+.12e}  Total={c.total:+.12e}"
    )
    lines.append("")
    lines.append("Note: these are algebraic regression tests with analytic toy PDFs, not physical predictions.")
    txt_path.write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    write_outputs()
