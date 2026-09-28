# Coefficients and input scope

`data/coefficients_demo.json` is the numerical input to the unexpanded demonstration.
It records massless fixed-flavor QCD with n_f=5 and a=alpha_s/(4π), using the CSS2
D=-K and d/d ln(mu) convention.

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

The `Coefficients` container accepts up to four beta and cusp entries and three
noncusp/CS-boundary entries. **Array length and the
`physical_n3llprime_payload` metadata field do not validate a physics payload.**
The shipped field is false; no operator routine converts it into a certificate.

An instantiated N³LL′ prediction additionally needs convention-checked physical
coefficients at the required orders, three-loop hard and TMD matching boundaries,
matching-log collinear evolution, and process/observable validation. The
`mock_import` suite validates a record schema with explicitly synthetic values;
those values are never fed into the numerical QCD demo.

The [reference guide](REFERENCES.md) lists published sources. This update does not
fetch external coefficients during tests or execute downloaded symbolic programs.
A later importer should verify source bytes, parse a declared symbolic subset,
record transformations, retain x-space distributions, and compare independent
physical expansions before using any accuracy label for an actual spectrum.
