# Coefficient source staging area

This directory is intentionally empty in the repository archive.

Use

```bash
bash scripts/stage_public_coefficients.sh
```

in a networked environment to stage the public ancillary coefficient files needed for a future implementation-level \(N^3\mathrm{LL}'\) claim.  The formal validation suite will pass without these files, but the physical coefficient gate will report `physics_ready=False` until mandatory source-stamped payloads are present and parsed.
