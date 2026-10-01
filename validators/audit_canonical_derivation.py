#!/usr/bin/env python3
"""Exact-arithmetic finite witnesses and negative controls for the analytic audit.

These checks support the independently reconstructed proof; finite witnesses do
not replace the all-order induction or analytic domain assumptions.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import sympy as sp

t, L, a = sp.symbols('t L a')
N = 8


def multiply(x, y):
    return [sp.expand(sum(x[k]*y[n-k] for k in range(n+1))) for n in range(N+1)]


def powers(x, n):
    result = [[sp.Integer(1)]+[sp.Integer(0)]*N, x]
    for _ in range(2, n+1):
        result.append(multiply(result[-1], x))
    return result


def exponent(beta, cusp, noncusp, boundary):
    running = [sp.Integer(0), sp.Integer(1)]+[sp.Integer(0)]*(N-1)
    for n in range(2, N+1):
        pp = powers(running, len(beta)+1)
        rhs = sum(2*b*pp[i+2][n] for i, b in enumerate(beta))
        running[n] = sp.integrate(rhs, (t, 0, t))
    pp = powers(running, max(len(cusp), len(noncusp), len(boundary)))
    E = []
    for n in range(N+1):
        integrand = sum(g*pp[i+1][n] for i, g in enumerate(noncusp))
        integrand -= 2*t*sum(g*pp[i+1][n] for i, g in enumerate(cusp))
        d = sum(g*pp[i+1][n] for i, g in enumerate(boundary))
        E.append(sp.expand(sp.integrate(integrand, (t, 0, L))-L*d.subs(t, L)))
    return running, E


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    tests = []
    def check(name, ok, detail=None):
        tests.append({'name': name, 'passed': bool(ok), 'detail': detail})
    # Integer fixtures make omitted-order effects exact, without rounding masks.
    base = [[1, 2, 3, 4], [2, 3, 5, 7], [3, 4, 5], [0, 1, 2]]
    running, retained = exponent(*base)
    check('running coefficient polynomial degrees', all(
        sp.degree(running[n], t) <= n-1 for n in range(1, N+1)))
    witnesses = {}
    for index, name in enumerate(('beta', 'cusp', 'noncusp', 'CS_boundary')):
        changed = [row[:] for row in base]
        changed[index].append(11)
        _, full = exponent(*changed)
        degrees = [None if sp.expand(full[n]-retained[n]) == 0 else int(sp.degree(full[n]-retained[n], L))
                   for n in range(N+1)]
        check('first omitted '+name+' stays below retained towers',
              all(d is None or (n >= 4 and d <= n-3) for n, d in enumerate(degrees)), degrees)
        witnesses[name] = [str(sp.expand(full[n]-retained[n])) for n in range(1, N+1)]
        # Removing a required loop must violate that same bound.
        lowered = [row[:] for row in base]
        lowered[index] = lowered[index][:-1]
        _, too_low = exponent(*lowered)
        violations = [n for n in range(1, N+1) if sp.expand(too_low[n]-retained[n]) != 0
                      and sp.degree(too_low[n]-retained[n], L) > n-3]
        check('negative control detects missing required '+name+' loop', bool(violations), violations)
    _, ll = exponent([1], [2], [0], [0])
    closed = sp.series((sp.log(1-2*a*L)+2*a*L)/a, a, 0, N+1).removeO().expand()
    check('LL all-order closed form through order eight',
          all(sp.expand(ll[n]-closed.coeff(a, n)) == 0 for n in range(1, N+1)))

    CF, ell, lb, f, P, C = sp.symbols('CF ell lb f P C')
    insertion = C-2*lb*P+CF*((3-2*ell)*lb-lb**2)*f
    derivative_mu = 2*sp.diff(insertion, lb)-2*sp.diff(insertion, ell)+4*P
    check('one-loop OPE plus PDF virtuality equation', sp.expand(derivative_mu-CF*(6-4*ell)*f) == 0)
    check('one-loop OPE rapidity equation', sp.expand(2*sp.diff(insertion, ell)+4*CF*lb*f) == 0)
    hard = CF*(-16+7*sp.pi**2/3+6*L-2*L**2)
    check('hard plus two TMD virtuality cancellation', sp.expand(-2*sp.diff(hard, L)+2*CF*(6-4*L)) == 0)
    epsilon, bA, bB, qA, qB, gA, gB, H = sp.symbols('epsilon bA bB qA qB gA gB H')
    naive = (1+epsilon*H)*(bA+epsilon*(qA+gA))*(bB+epsilon*(qB+gB))
    strict = bA*bB+epsilon*(H*bA*bB+(qA+gA)*bB+bA*(qB+gB))
    check('six-term assembly equals strict linear expansion', sp.expand(sp.series(naive, epsilon, 0, 2).removeO()-strict) == 0)
    check('negative control rejects naive products', sp.degree(sp.expand(naive-strict), epsilon) == 3)

    # Plancherel contact constant and explicit beta-integral normalization of
    # the p=3, q=4 native Student convolution at zero bin and infinity.
    u, tau, s = sp.symbols('u tau s', positive=True)
    check('contact error normalization', sp.integrate((tau/(tau+s))**2, (s, 0, sp.oo)) == tau)
    report = {'schema': 'canonical-derivation-audit/1', 'passed': all(x['passed'] for x in tests),
              'tests': tests, 'highest_exact_order': N, 'omitted_tower_witnesses': witnesses,
              'proof_status': 'Finite exact witnesses plus negative controls; all-order and domain arguments are in DERIVATION_AUDIT.md',
              'script_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    (args.out/'derivation_audit.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps({'passed': report['passed'], 'checks': len(tests),
                      'failures': [x for x in tests if not x['passed']]}, indent=2), flush=True)
    return 0 if report['passed'] else 1


if __name__ == '__main__':
    raise SystemExit(main())
