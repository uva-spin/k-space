# Validation checks

This checklist summarizes the formal checks implemented in the `validators/` directory.  The checks are algebraic, symbolic, or infrastructure-oriented.  They do not constitute a physical TMD extraction or a physical cross-section prediction.

## Expected global status

Running

```bash
python run_all_checks.py --clean
```

should produce

```text
formal_suite_passed: True
checks: 12/12
```

The physical coefficient and implementation-level label gates should remain blocked unless external NNLO/N$^3$LO coefficient payloads have been staged and validated.

## Distribution and one-loop checks

1. **Two-dimensional plus-bin functionals**
   - Checks $\mathcal B_n(0,Q;Q)=0$.
   - Checks bin additivity.
   - Checks the off-origin relation to the ordinary one-loop singular kernel.

2. **Toy-PDF $x$-convolution layer**
   - Checks the longitudinal plus prescription for $P_{qq}^{(0)}$.
   - Checks flavor-vector luminosity assembly.
   - Checks bin additivity after longitudinal and transverse distributions are combined.

## Convolution algebra and evolution

3. **Generator algebra**
   - Checks $\mathcal G_\eta\otimes_T\mathcal G_\xi=\mathcal G_{\eta+\xi}$.
   - Generates the $\mathcal L_m\otimes_T\mathcal L_n$ table.

4. **One-loop evolution Green function**
   - Checks path independence in $(\mu,\zeta)$.
   - Checks rapidity and $\mu$ evolution equations.
   - Checks semigroup composition.

5. **Two-loop ingredients**
   - Checks the two-loop Collins-Soper-kernel projection.
   - Checks the scalar NNLO logarithmic reconstruction skeleton.

6. **Running-coupling evolution**
   - Checks endpoint RG-integral identities through $\mathcal O(a_s^2)$.
   - Checks projection back to the rapidity exponent.

## W-term and accuracy labeling

7. **Binned W-term prototype**
   - Checks one-loop profile-scale cancellation.
   - Checks bin-level equality of profile and canonical forms.

8. **Accuracy ingredient manifest**
   - Encodes loop-order requirements for $N^k\mathrm{LL}'$ bookkeeping.
   - For $N^3\mathrm{LL}'$, requires four-loop cusp/beta, three-loop noncusp/CS/hard/matching inputs, and three-loop DGLAP splitting kernels.

9. **N$^3$LO matching-import interface**
   - Validates the coefficient schema and projection map with synthetic mock data.
   - Explicitly does not authorize mock coefficients for physical claims.

10. **Physical coefficient source gate**
    - Audits whether expected external coefficient payloads have been staged.
    - Expected status in the formalism repository: infrastructure passes, physics gate is blocked.

11. **$N^3\mathrm{LL}'$ label certificate**
    - Confirms formal prescription readiness.
    - Blocks implementation-level claim until physical payloads and $\mathcal O(a_s^3)$ expansion validation pass.

12. **Formalism completion**
    - Checks that the completed formalism note contains generic process, nonperturbative module, regulated normalization, $W+Y$, claim-discipline, validation, and closure sections.
