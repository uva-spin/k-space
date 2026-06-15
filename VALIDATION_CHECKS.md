# Validation checklist

This file records the checks that support the formalism paper.  They are deliberately algebraic and infrastructure-oriented; they do not constitute a physical data extraction.

## Distribution and one-loop checks

1. **Two-dimensional plus-bin functionals**
   - Checks \(\mathcal B_n(0,Q;Q)=0\).
   - Checks bin additivity.
   - Checks the off-origin relation to the ordinary one-loop singular kernel.

2. **Toy-PDF \(x\)-convolution layer**
   - Checks the longitudinal plus prescription for \(P_{qq}^{(0)}\).
   - Checks flavor-vector luminosity assembly.
   - Checks bin additivity after longitudinal and transverse distributions are combined.

## Convolution algebra and evolution

3. **Generator algebra**
   - Checks \(\mathcal G_\eta\otimes_T\mathcal G_\xi=\mathcal G_{\eta+\xi}\).
   - Generates the \(\mathcal L_m\otimes_T\mathcal L_n\) table.

4. **One-loop evolution Green function**
   - Checks path independence in \((\mu,\zeta)\).
   - Checks rapidity and \(\mu\)-evolution equations.
   - Checks semigroup composition.

5. **Two-loop ingredients**
   - Checks the two-loop Collins--Soper kernel projection.
   - Checks the scalar NNLO logarithmic reconstruction skeleton.

6. **Running-coupling evolution**
   - Checks endpoint RG-integral identities through \(\mathcal O(a_s^2)\).
   - Checks projection back to the rapidity exponent.

## W-term and accuracy labeling

7. **Binned W-term prototype**
   - Checks one-loop profile-scale cancellation.
   - Checks bin-level equality of profile and canonical forms.

8. **Accuracy ingredient manifest**
   - Encodes the loop-order requirements for \(N^k\mathrm{LL}'\).
   - For \(N^3\mathrm{LL}'\), requires four-loop cusp/beta, three-loop noncusp/CS/matching/hard inputs, and three-loop DGLAP splitting kernels.

9. **N\(^3\)LO matching-import interface**
   - Validates the schema and projection map with mock data.
   - Explicitly does not authorize use of mock coefficients for physical claims.

10. **Physical coefficient source gate**
    - Audits whether expected external coefficient payloads have been staged.
    - Expected status in this repository: infrastructure passes, physics gate is blocked.

11. **\(N^3\mathrm{LL}'\) label certificate**
    - Confirms formal prescription readiness.
    - Blocks implementation-level claim until physical payloads and \(a_s^3\) expansion validation pass.

12. **Formalism completion**
    - Checks that the completed formalism note contains generic process, nonperturbative module, regulated normalization, \(W+Y\), claim-discipline, validation, and closure sections.
