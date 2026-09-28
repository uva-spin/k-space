#!/usr/bin/env bash
# Optional source acquisition for future physical coefficient validation.
# Store ancillary PDFs and source archives outside this code repository.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
if [[ "$#" -ne 1 ]]; then
  echo "Usage: bash scripts/stage_public_coefficients.sh /absolute/external/staging-directory" >&2
  exit 2
fi
DEST="$(python3 - "$1" "$REPO_ROOT" <<'PYCODE'
from pathlib import Path
import sys
p=Path(sys.argv[1]).expanduser().resolve()
r=Path(sys.argv[2]).resolve()
if p.is_relative_to(r) or r.is_relative_to(p):
    raise SystemExit('Choose a staging directory outside and not an ancestor of the repository')
print(p)
PYCODE
)"
mkdir -p "$DEST/2006.05329v2" "$DEST/2012.03256v2"

fetch() {
  local url="$1"
  local out="$2"
  if [[ -s "$out" ]]; then
    echo "exists: $out"
    return 0
  fi
  echo "fetch: $url -> $out"
  curl -fL --retry 4 --retry-delay 2 "$url" -o "$out"
}

# Direct ancillary URLs documented by the arXiv records. If arXiv changes these paths,
# use the fallback e-print source packages below and extract the anc/ directory manually.
EMV_BASE="https://arxiv.org/src/2006.05329v2/anc"
LYZZ_BASE="https://arxiv.org/src/2012.03256v2/anc"

# EMV2020: Transverse momentum dependent PDFs at N^3LO.
for f in TMDPDF.m PTBF.m PTBF_ZExpansion.m PTBF_ZbExpansion.m PT_SoftFunction.m PT_SoftFunction_Renormalized.m info.txt; do
  fetch "$EMV_BASE/$f" "$DEST/2006.05329v2/$f" || true
done

# LYZZ2021: Unpolarized Quark and Gluon TMD PDFs and FFs at N^3LO.
for f in BeamFunction.m BeamfunctionN.m FFMatchingKernels.m FFMatchingKernelsN.m TMDFF.m TMDFFN.m TMDPDF.m TMDPDFN.m ancillary_readme.pdf resummedFFsinglets.m softfunction.m; do
  fetch "$LYZZ_BASE/$f" "$DEST/2012.03256v2/$f" || true
done

# Fallback: source packages. These should contain the TeX source and ancillary directory.
fetch "https://arxiv.org/e-print/2006.05329v2" "$DEST/2006.05329v2/e-print-2006.05329v2.tar" || true
fetch "https://arxiv.org/e-print/2012.03256v2" "$DEST/2012.03256v2/e-print-2012.03256v2.tar" || true

# Exclude the ledger itself on repeated runs; use portable SHA256 and relative names.
python3 - "$DEST" <<'PYCODE'
import hashlib
from pathlib import Path
import sys
root=Path(sys.argv[1])
ledger=root/'SHA256SUMS.txt'
rows=[]
for path in sorted(root.rglob('*')):
    if path==ledger or path.is_symlink() or not path.is_file():
        continue
    rows.append(hashlib.sha256(path.read_bytes()).hexdigest()+'  '+path.relative_to(root).as_posix())
ledger.write_text('\n'.join(rows)+'\n')
PYCODE

echo "Wrote hash ledger: $DEST/SHA256SUMS.txt"
