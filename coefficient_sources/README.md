# Physical coefficient source acquisition

The scientific checks use the declared input in `data/coefficients_demo.json`
and explicitly labeled formal fixtures. Optional ancillary acquisition is separate:

```bash
bash scripts/stage_public_coefficients.sh /absolute/external/staging-directory
```

The destination must be outside this repository because upstream source archives
and ancillary PDFs belong to the source/writing workspace. Acquisition and hashes
establish file provenance; physical validation requires parsing, convention checks,
and independent perturbative comparisons. See [coefficients](../docs/COEFFICIENTS.md).

The script excludes its checksum ledger from the hashed payload on repeated runs.
Acquisition is optional and is never called by verification. Check each upstream
source's license and attribution terms before redistribution or coefficient reuse;
the repository license does not grant rights to downloaded third-party payloads.
