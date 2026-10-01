# Method-paper claim and scope

The paper establishes a native transverse-momentum representation and numerical
realization of the Collins/CSS2 leading-power singular TMD operator on a declared
domain. The canonical branch carries N³LL evolution using four-loop beta and cusp
and three-loop noncusp and canonical Collins-Soper inputs. Its action on an
admitted momentum-space input is evaluated by momentum-space convolutions and
elliptic resolvent solves; the production action does not Fourier- or
Hankel-transform the fitted input.

The current physical boundary assembly is strict NLO. The accurate label is
**N³LL evolution with strict NLO hard and OPE boundaries**. The package and paper
do not claim an N³LL′ spectrum. Collins/CSS2 factorization and its RG equations
are cited inputs rather than a new factorization proof.

## Branches that must remain separate

The uncapped spectral-canonical branch is the subject of the retained-tower,
bounded-symbol, and resolvent theorems. Its analytic endpoint is zero at spectral
`s=0`; its infrared completion is endpoint independent.

The capped and noncanonical branches are process-level diagnostics. They have
finite, generally nonzero endpoint values and use separate profile assumptions.
They do not inherit the canonical zero endpoint or the canonical N³LL tower
theorem. They may appear only as bounded physical demonstrations with their own
tests and qualifications.

The [canonical physical audit](CANONICAL_PHYSICAL_AUDIT.md) now provides a
separate path that combines the actual uncapped canonical operator with physical
PDFs and all six strict-NLO components. Its independent continuum reference,
two-leg algebraic-tail inputs, gradients, contacts, signed bins and PDF-endpoint
bound directly test this common prescription. That audit has its own frozen
numerical acceptance criteria; it is not inferred from the capped demonstrator.

## Publication claims

The paper may state that:

- the imported CSS2 singular operator has the stated native-kT representation on
  its admitted distribution/operator domain;
- the spectral-canonical solution retains the conventional N³LL evolution towers
  with the declared 4/4/3/3-loop inputs and weak-coupling window;
- origin contacts, ordinary inputs, and finite-bin actions are represented by the
  contact-aware functional calculus;
- the fitted input can remain in momentum space during the production operator
  action;
- actual PDFs, charge/flavor sums, and strict-NLO hard/OPE insertions can be used
  in a bounded Drell-Yan demonstrator.

Claims about N³LL′ boundaries, complete W+Y matching, massive threshold
completion, a global extraction, or a representation advantage in fitted results
belong to later work.

The machine-readable statement is
[`manifests/method_paper_accuracy.json`](../manifests/method_paper_accuracy.json).
The executable physical canonical check is
[`validators/validate_physical_canonical.py`](../validators/validate_physical_canonical.py).
