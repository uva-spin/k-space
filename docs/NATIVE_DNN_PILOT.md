# Native-kT mixture-density DNN pilot

This implementation learns an x-dependent nonperturbative boundary directly in
kT while retaining the declared E866 massless nf5 W-only prescription. It is a
new model analysis, not a reproduction of the reference paper's b-space prior.
The first experimental sample remains the 41 registered E866 proton-proton bins.

## Learned input and fixed physics

At Q0=5 GeV the shared boundary is

```math
S_\theta(x,k)=\sum_{i=1}^{6} a_i(x;\theta)
 \frac{e^{-k^2/v_i}}{\pi v_i},\qquad
 a_i=\operatorname{softmax}_i\!\left(\operatorname{MLP}
   \left[\frac{\ln(x/0.1)}{\ln4}\right]\right).
```

The network has three tanh hidden layers with 12 units each, and six output
logits: 414 trainable parameters. The basis variances are 0.08, 0.20, 0.50,
1.25, 3.125 and 7.8125 GeV^2. Each input is positive and normalized with the
2*pi*k*dk measure. This admits non-Gaussian shapes and x-kT correlations, but is
narrower than an arbitrary positive kT function. The fixed basis and its tails
are explicit model assumptions. They do not imply positivity of evolved TMDs;
no evolved function or cross section is clipped.

The same boundary is shared across both proton legs and all active nf5 quark
channels, including massless charm and bottom. No separate heavy-flavor input
is inferred.
The collinear PDF, perturbative coupling, strict Born/hard/qq/qg insertions and
profile evolution remain fixed. The NP-CS coefficient is fixed at
0.15790301183804753 GeV^2, the conditional Gaussian diagnostic's value. This is
a registered choice, not an independent determination of a universal CS kernel.
The profile and measurement limitations in the E866 baseline remain in force.

The pair input is evaluated by the exact native-kT convolution identity

```math
G(\cdot;v_i)*G(\cdot;v_j)=G(\cdot;v_i+v_j),\qquad
v_{ij}(Q)=v_i+v_j+8g\ln(Q/Q_0).
```

The fixed quadratic CS contribution is therefore a Gaussian convolution in
momentum space. No learned input or production training gradient is Fourier
transformed. Analytic transforms appear only in the independent validation path.

## Native operator and bin conversion

The operator retains xF quadrature nodes until after evaluating the learned
weights at xA and xB. The previous x-independent Gaussian shortcut is not used.
Each strict insertion has its own signed rational coefficients and explicit
contact term. Shared shifted finite-volume elliptic solves act on the native
pair basis. Linearity permits these responses to be precomputed; neural training
then contracts the physical response tensors with the bilinear leg weights.
Their reverse-mode gradients are analytic and tested against finite differences.

The registered conditional experimental estimator is

```math
I_r=\frac{E_{r,\mathrm{tab}}}
 {\bar q_r\,\Delta q_r\,\Delta x_{F,r}}
 \int_r dQ^2\,dx_F\,dq\;
 \frac{q}{\sqrt{Q^2+q^2+s x_F^2/4}} I_Q(Q,x_F,q).
```

The tabulated mean energies are preserved, despite their unresolved physical
interpretation. This convention is an explicit analysis hypothesis, not a claim
that the experiment's original estimator code has been authenticated. Mass is
integrated, not averaged. Cell weights integrate q/E exactly, including the
origin. This adapter is distinct from a uniform average of the physical invariant.

## Optimization and regularization

The experimental objective is half the diagonal data chi-square, a 6.5% Gaussian
normalization prior, weight regularization and a functional smoothness penalty:

```math
L=\frac12\sum_r\left[\frac{N I_r-d_r}{\sigma_r}\right]^2
 +\frac12\left[\frac{N-1}{0.065}\right]^2
 +\frac{\lambda}{2}\|\theta\|^2
 +\frac{\eta}{2}\operatorname{mean}_j\sum_i
       \left[D^2_{\ln x}a_i(x_j)\right]^2.
```

N=exp(eta_N) is restricted to [0.5,1.5]. Statistical and point-systematic errors
are combined in quadrature and treated as diagonal. The functional penalty uses
257 logarithmic points on x in [0.04,0.70], centered second differences and
coefficient eta=0.1. Weight regularization includes biases, excludes N, and uses
lambda=1e-5 in noiseless closure and 1e-4 in the experimental pilot. These are
model-prior choices, not calibrated confidence prescriptions.

Initial unregularized-in-function fits exploited x variation unresolved by the
coarser Q/xF grid. Those results are rejected and retained externally. The
functional prior and denser 8-by-16 longitudinal quadrature address that failure;
12-by-24 independent reference quadrature checks the final learned inputs.
The numerical tolerances are not relaxed to admit a smaller chi-square.

Optimization uses Adam followed by L-BFGS-B. Seeds and convergence diagnostics
are retained. The best seed is chosen using training objective only. No number
of neural parameters is subtracted from the data count to report a misleading
chi-square per degree of freedom. Initialization spread is not an experimental
uncertainty band or a calibrated confidence interval.

## Closure and validation

Independent continuum predictions from an analytic, x-correlated non-Gaussian
mixture generate the noiseless closure targets. The network receives their
observables, not the generating formula. All ten bins in the 7.2-8.7 GeV mass
window are withheld; 31 bins train the network. The training and held-out
predictions, input profiles and weak sensitivity directions are reported.
A good cross-section reconstruction does not establish unique full-TMD recovery.

A local sensitivity matrix uses eight log-x hat functions for five relative
mixture logits. Its singular values describe identifiable and weak directions in
those declared function coordinates; neural weight redundancy is not mistaken
for physical information. The global x/CS ambiguity is controlled in this pilot
by fixing the CS coefficient, not by claiming that a DNN removes it.

The checks include every basis pair and separate insertion, direct kT pair
convolution, input normalization, contact endpoints, full neural gradients and
post-training integration refinement. The measured numerical errors are finite
checks on this model class, not global continuum certificates. Flexible basis,
CS, PDF, measurement and regularizer variations remain distinct systematic
studies; this first pilot does not provide a full global uncertainty extraction.

## Running and preserving results

The local runtime and PDF resources are inherited from the completed resolution
and prescription runs. `build_native_dnn_kernel.py` creates the native response
and independent reference tensors. `run_native_dnn.py` runs closure before data;
`polish_dnn.py` continues saved, numerically validated experimental fits with the
same objective and tolerances. The continuation retains the previous optimizer
record. Neither an iteration limit nor a small change in chi-square is silently
reclassified as convergence. The external run's reproduction script creates a
fresh sibling directory and leaves accepted evidence intact.

`reconstruct_native_tmd.py` exports each strict single-leg matching component for
all active quark flavors at mu=Q and zeta=Q^2, excluding the hard factor. Those
values have an independent reference check. `audit_dnn_pilot.py` reports input
recovery, initialization spread and regularizer-grid refinement. A change to the
mean-energy conversion in that audit is reconstruction-only with normalization
profiled; it is not a full neural refit or a calibrated systematic uncertainty.
