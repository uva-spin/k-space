# Coefficients and input scope

The package contains two declared numerical inputs. Both use massless fixed-flavor
QCD with a=alpha_s/(4π), the CSS2 D=-K convention, and derivatives with respect
to ln(mu).

`data/coefficients_demo.json` is the lower-order input to the original
unexpanded demonstration.

| Vector | Supplied coefficients |
|---|---|
| `beta` | beta0, beta1 |
| `cusp` | Gamma0, Gamma1 for the quark representation |
| `noncusp` | gammaF0 |
| `d` | d1=0, d2 |

The reference initialization is mu0=2 GeV and a(mu0)=0.02. The code's default RG
integration interval is t=ln(nu/mu0) in [0,1200], with DOP853, rtol=2×10⁻¹²,
and atol=2×10⁻¹⁴. Values outside the solved interval require a separately extended
solve rather than extrapolation.

`data/coefficients_css2_nf5_n3ll.json` supplies the physical nf=5 evolution
vectors at 4/4/3/3 loops: beta and quark cusp through four loops, and TMD
noncusp and canonical CS boundary through three loops. Its source versions,
hashes, equation locators, and convention changes are pinned in
`manifests/css2_sources.json`. The `physical_canonical` suite compares the JSON
against the independent formula implementation in `src/css2_coefficients.py`
and exercises it in the canonical operator.

The `Coefficients` container accepts up to four beta and cusp entries and three
noncusp/CS-boundary entries. N3LL metadata requires exactly 4/4/3/3 entries.
**Array length and metadata alone do not validate a physics payload.** Source
qualification, convention checks, independent coefficient comparisons, and RG
tests provide the evidence for the evolution input.

The current physical input supports N³LL evolution with strict NLO hard and OPE
boundaries. An instantiated N³LL′ prediction additionally needs three-loop hard
and TMD matching boundaries,
matching-log collinear evolution, and process/observable validation. The
`mock_import` suite validates a record schema with explicitly synthetic values;
those values are never fed into the numerical QCD demo.

The [reference guide](REFERENCES.md) lists published sources. This update does not
fetch external coefficients during tests or execute downloaded symbolic programs.
A later importer should verify source bytes, parse a declared symbolic subset,
record transformations, retain x-space distributions, and compare independent
physical expansions before using any accuracy label for an actual spectrum.
