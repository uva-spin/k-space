#!/usr/bin/env python3
"""Independent primary-source cusp cross-check; never evaluates downloaded code."""
import argparse,json,hashlib,re,sys
from pathlib import Path
import sympy as s
from sympy.parsing.mathematica import parse_mathematica
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from css2_coefficients import qcd_coefficients
p=argparse.ArgumentParser();p.add_argument('--sources',type=Path,required=True);p.add_argument('--out',type=Path,required=True);args=p.parse_args()
path=args.sources/'hkm_cusp/anc/HKM_GammaC.txt';raw=path.read_text()
# Only this elementary expression grammar is accepted, rather than evaluating
# arbitrary Mathematica or Python. Reject unexpected tokens before parsing.
words=set(re.findall(r'[A-Za-z][A-Za-z0-9_]*',raw))
assert words <= {'aS','CA','CR','CF','NF','TF','C4A','C4F','Z'}
assert not re.search(r'[^\w\s+*/^().,\[\]-]',raw)
e=s.expand(parse_mathematica(raw))
Z=s.Function('Z');e=e.replace(lambda z:z.func==Z,lambda z:s.zeta(z.args[0]))
N=s.Integer(3)
# HKM printed color invariants use N_A; ancillary requires N_R=N_F.
sub={s.Symbol('CA'):N,s.Symbol('CF'):(N*N-1)/(2*N),s.Symbol('CR'):(N*N-1)/(2*N),
     s.Symbol('TF'):s.Rational(1,2),s.Symbol('C4A'):N*(N*N+6)/48*(N*N-1)/N,
     s.Symbol('C4F'):(N**4-6*N*N+18)/(96*N*N)*(N*N-1)/N}
rows=[]
for nf in (3,4,5,6):
 original=e.subs(sub).subs(s.Symbol('NF'),nf)
 # HKM uses alpha_s/pi = 4*a, MPS implementation uses alpha_s/(4*pi).
 values=[float((original.coeff(s.Symbol('aS'),i)*4**i).evalf(40)) for i in range(1,5)]
 got=list(qcd_coefficients(nf).cusp)
 residual=max(abs(a-b) for a,b in zip(values,got))
 rows.append({'nf':nf,'hkm_40_digit_evaluation':values,'implementation_mps':got,'max_abs_difference':residual,'passed':residual<2e-9})
result={'primary_ancillary_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),'comparison':rows,'passed':all(r['passed'] for r in rows)}
args.out.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
raise SystemExit(0 if result['passed'] else 1)
