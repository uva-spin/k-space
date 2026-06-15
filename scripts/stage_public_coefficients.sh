#!/usr/bin/env bash
# Stage public high-order TMD matching coefficient ancillary files for the KCSS N3LL' certificate.
# Run this in an environment with network access. The current ChatGPT code sandbox has no live network.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
DEST="${1:-$REPO_ROOT/coefficient_sources}"
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

# Record hashes for all files that were actually acquired.
(
  cd "$DEST"
  find . -type f -maxdepth 3 -print0 | sort -z | xargs -0 sha256sum > SHA256SUMS.txt
)

echo "Wrote hash ledger: $DEST/SHA256SUMS.txt"
