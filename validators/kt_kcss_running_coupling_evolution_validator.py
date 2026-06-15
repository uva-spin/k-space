#!/usr/bin/env python3
"""
Symbolic validator for the running-coupling K_CSS2 evolution kernel.

The checks are performed through O(a_i^2), the order retained in the v0.11
running-coupling endpoint kernel.  This is a formal regression validator, not a
phenomenology calculation.
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
import sympy as sp


def main() -> None:
    a, L, rho, ell, X = sp.symbols('a L rho ell X')
    L1, L2, r1, r2 = sp.symbols('L1 L2 r1 r2')
    G0, G1, g0, g1, b0, d2 = sp.symbols('G0 G1 g0 g1 b0 d2')
    L0dist, L1dist, delta = sp.symbols('L0dist L1dist delta')
    fr0, fr1 = sp.symbols('fr0 fr1')

    def trunc(expr):
        return sp.series(sp.expand(expr), a, 0, 3).removeO().expand()

    def D_b(aa, XX):
        return trunc(aa*G0*XX + aa**2*(sp.Rational(1, 2)*b0*G0*XX**2 + G1*XX + d2))

    def A_Gamma(aa, LL):
        return trunc(aa*G0*LL + aa**2*(G1*LL - b0*G0*LL**2))

    def A_gamma(aa, LL):
        return trunc(aa*g0*LL + aa**2*(g1*LL - b0*g0*LL**2))

    def S_Gamma(aa, LL):
        return trunc(sp.Rational(1, 2)*aa*G0*LL**2 + aa**2*(sp.Rational(1, 2)*G1*LL**2 - sp.Rational(2, 3)*b0*G0*LL**3))

    def a_final(aa, LL):
        return trunc(aa - 2*b0*aa**2*LL)

    def R_exp(aa, LL, eell):
        return trunc(A_gamma(aa, LL) - eell*A_Gamma(aa, LL) + 2*S_Gamma(aa, LL))

    def endpoint_exp(aa, LL, rr, eell, XX):
        af = a_final(aa, LL)
        return trunc(R_exp(aa, LL, eell) - rr*D_b(af, XX + 2*LL))

    af = a_final(a, L)
    D_i = D_b(a, X)
    D_f = D_b(af, X + 2*L)
    AG = A_Gamma(a, L)
    R_mu = R_exp(a, L, ell)
    E = trunc(R_mu - rho*D_f)

    # Integrated CS-kernel RG identity and path independence.
    delta_D = sp.simplify(trunc(D_f - D_i - 2*AG))
    R_zeta_f = R_exp(a, L, ell + 2*rho)
    delta_path = sp.simplify(trunc((R_mu - rho*D_f) - (R_zeta_f - rho*D_i)))

    # Fixed-coupling one-loop reduction.
    E_fixed_expected = a*L*(g0 - G0*(ell - L)) - a*rho*G0*(X + 2*L)
    E_fixed_actual = trunc(E.subs({b0: 0, G1: 0, g1: 0, d2: 0}))
    delta_fixed = sp.simplify(trunc(E_fixed_actual - E_fixed_expected))

    # Rapidity endpoint equation at the exponent level: dE/drho = -D(mu_f).
    delta_rho = sp.simplify(trunc(sp.diff(E, rho) + D_f))

    # Mu endpoint equation at fixed zeta_f.  Since ell_f=ell+2rho-2L,
    # dE/dL should equal gamma_i(a_f)-Gamma_i(a_f)*ell_f through O(a^2).
    gamma_f = trunc((af*g0 + af**2*g1) - (af*G0 + af**2*G1)*(ell + 2*rho - 2*L))
    delta_mu = sp.simplify(trunc(sp.diff(E, L) - gamma_f))

    # Semigroup through O(a^2).
    E10 = endpoint_exp(a, L1, r1, ell, X)
    a1 = a_final(a, L1)
    ell1 = ell + 2*r1 - 2*L1
    X1 = X + 2*L1
    E21 = endpoint_exp(a1, L2, r2, ell1, X1)
    E20 = endpoint_exp(a, L1 + L2, r1 + r2, ell, X)
    delta_semi = sp.simplify(trunc(E20 - E10 - E21))

    # Projection from frak basis to conventional basis.
    U_frak = (delta
              - a*rho*G0*fr0
              + a**2*((sp.Rational(1, 2)*rho**2*G0**2 - sp.Rational(1, 2)*rho*b0*G0)*fr1
                      - rho*G1*fr0 - rho*d2*delta))
    U_conv_from_frak = sp.expand(U_frak.subs({fr0: -L0dist, fr1: -2*L1dist}))
    U_conv_expected = (delta + a*rho*G0*L0dist
                       + a**2*(rho*G1*L0dist + (rho*b0*G0 - rho**2*G0**2)*L1dist - rho*d2*delta))
    delta_proj = sp.simplify(sp.expand(U_conv_from_frak - U_conv_expected))

    tests = {
        'fixed': delta_fixed,
        'D': delta_D,
        'path': delta_path,
        'rho': delta_rho,
        'mu': delta_mu,
        'semi': delta_semi,
        'rapidity_proj': delta_proj,
    }
    rows = []
    for name, expr in tests.items():
        simplified = sp.simplify(expr)
        rows.append({'test': name, 'delta_expr': str(simplified), 'passed': simplified == 0})

    # Quark-like smoke point. gamma1 is an explicit structural input, not a claim.
    CA = 3.0
    CF = 4.0 / 3.0
    TF = 0.5
    nf = 5.0
    beta0 = (11*CA - 4*TF*nf) / 3
    Gamma0 = 4*CF
    Gamma1 = 4*CF*((67/9 - math.pi**2/3)*CA - (20/9)*TF*nf)
    gamma0 = 6*CF
    gamma1 = 17.25
    d2num = CF*(CA*(808/27 - 28*float(sp.zeta(3).evalf())) - (224/27)*TF*nf)
    vals = {a: 0.018, L: 0.58, rho: 0.24, ell: 0.31, X: -0.52,
            G0: Gamma0, G1: Gamma1, g0: gamma0, g1: gamma1, b0: beta0, d2: d2num}
    R_q = float(R_mu.subs(vals))
    af_num = float(af.subs(vals))
    eta_q = float((rho*(af*G0 + af**2*G1)).subs(vals))
    kappa_q = float((rho*af**2*b0*G0/2).subs(vals))

    passed = all(row['passed'] for row in rows)
    output = {
        'passed': passed,
        'tests': rows,
        'smoke': {
            'R_q': R_q,
            'eta_q': eta_q,
            'kappa_q': kappa_q,
            'a_f': af_num,
            'parameters': {'a_i': vals[a], 'L': vals[L], 'rho': vals[rho], 'ell_i': vals[ell], 'L_b_i': vals[X]},
        },
        'constants': {
            'nf': nf,
            'CF': CF,
            'CA': CA,
            'TF': TF,
            'beta0': beta0,
            'Gamma0q': Gamma0,
            'Gamma1q': Gamma1,
            'gamma0q': gamma0,
            'gamma1q_smoke_input': gamma1,
            'd2q': d2num,
        },
        'note': 'All algebraic identities are checked through O(a_i^2). The smoke gamma1 input is structural only.',
    }

    base = _OUTPUT_DIR
    json_path = base / 'kt_kcss_running_coupling_evolution_validator_output.json'
    csv_path = base / 'kt_kcss_running_coupling_evolution_validator_tests.csv'
    txt_path = base / 'kt_kcss_running_coupling_evolution_validator_report.txt'
    json_path.write_text(json.dumps(output, indent=2, sort_keys=True))
    with csv_path.open('w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=['test', 'delta_expr', 'passed'])
        writer.writeheader()
        writer.writerows(rows)
    lines = [
        'K_CSS2 running-coupling evolution symbolic validator',
        '=' * 60,
        f'passed: {passed}',
        '',
        'Tests:',
    ]
    for row in rows:
        lines.append(f"  {row['test']}: delta={row['delta_expr']}, passed={row['passed']}")
    lines += [
        '',
        'Quark-like smoke point:',
        f'  R_q = {R_q:.12e}',
        f'  eta_q = {eta_q:.12e}',
        f'  kappa_q = {kappa_q:.12e}',
        '',
        output['note'],
    ]
    txt_path.write_text('\n'.join(lines) + '\n')
    print(txt_path.read_text())


if __name__ == '__main__':
    main()
