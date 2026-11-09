# Sampling development progress

Directive: `PhaseSpaceTools_codex_directive.md`, supplied 10 November 2026.
Development branch: `codex/exact-quantum-sampling`. Starting revision:
`b180b50f0eddfed02d2f19728c6cbb5d15865b36` (v1.0.7).
No repository-specific AGENTS.md exists. Global Julia protections apply.

## Baseline

Julia 1.13.3, installed juliaup channel `1.13`; unmodified `Pkg.test()` passed,
including Aqua. The old tests did not detect the squeezing phase, Crescent
Q broadening, missing Bogoliubov P vacuum contribution, or erroneous signed
Fock number variance. Existing Julia channels/default were preserved.

## Stage 0 (complete)

Corrected the W squeezing rotation, Crescent Q convolution, Bogoliubov P
covariance/classicality check, exact Fock vacuum, signed Fock positive-P
factorial and number-variance tests, correlated doubled-coordinate wording,
Q ordering indices, and two positive-W convolution normalisation factors.
Statistical regressions check real and imaginary parts at six estimated
standard errors, with a fixed baseline seed.

Validation: Julia 1.13.3 `Pkg.test()` passed, including 104 scientific
regressions and Aqua.

## Scope

Stages 0–6 and the independent Wplus track are complete. Release, merge and
publication decisions remain separate; no release or merge was performed.

## Stage 1 (complete)

Added RNG-first methods with default wrappers; matrix Q sampling and in-place
entry points; legacy return adapters; finite-parameter/count checks; bounded
rejection diagnostics; pure-overlap capability dispatch and stable logarithmic
coherent, Fock and squeezed overlaps. Thermal positive-P now uses its narrow
ordinary P. Minimum Julia compatibility is explicitly 1.10 (the oldest installed
channel to be tested), rather than an unverified Julia 1.0 claim.

Validation: full Julia 1.13.3 package tests passed, including 102 infrastructure checks.

## Stage 2 (complete)

Added displacement (including the pure overlap phase), rotations, classical
mixtures, product states, finite number mixtures, photon-added thermal Q and
photon-subtracted thermal P/W/Q/positive-P. Probability laws are cached.
Tests use explicit MersenneTwister RNGs and independently known Gamma moments.

Validation: full Julia 1.13.3 tests passed, including 104 exact-family checks.

## Stage 3 (complete)

Implemented finite pure superpositions with cached Gram normalisation, |c|
proposal weights, scaled logarithmic amplitudes and bounded rejection; coherent,
Fock/product/displaced/rotated and pure squeezed components; cats with stable
parity/small-odd proposals; compass, sparse Fock, NOON, photon-added coherent
and squeezed-cat constructors. Small-odd envelope and accepted density are
independently derived in docs/src/exact-sampling.md. Finite-Fock validation uses
ladder matrices/exponentials independently of sampler formulas.

Validation: full Julia 1.13.3 suite passed with 189 interference/overlap checks;
90 additional independent pointwise odd-cat density/bound checks passed in the
warm session. Exact odd parity compares Float64 π explicitly (the Irrational
π equality edge case was caught and fixed).

## Stage 4 (complete)

Added explicit unitary passive mixing and the O(M) collective-orbital exact Q
sampler, plus specialised occupied-mode positive-P with zero orthogonal vacuum.
Generic Q conversion remains selectable. Tests cover complex orbitals, vacuum,
single/two/multiple modes, cross and fourth moments and an independent binomial
Hilbert expansion. No per-sample number rescaling is present.

Validation: full Julia 1.13.3 suite passed, including 349 multimode checks.

## Stage 5 (complete)

Added exact Q loss and Gaussian phase diffusion, scalar/per-mode parameters,
common versus independent phase noise, and generic positive-P inheritance.
Validation compares lossy Fock states to binomial mixtures and lossy cats to
independent Kraus density matrices; first/anomalous/higher phase factors and
identity/complete-loss limits are covered.

Validation: full Julia 1.13.3 suite passed, including 99 channel checks.

## Stage 6 (complete)

Added general Gaussian covariance validation/cached factors, semidefinite P,
squeezed-thermal/two-mode/multimode Bogoliubov constructors and the missing
legacy two-mode W path. Added exact centred photon-added Gaussian Q through a
nonnegative monomial mixture; its derivation and finite-density-matrix cutoff
checks are included. The user specifically asked about quadrature units:
explicit tests/documentation reconcile [x,p]=i with X=√2x,Y=√2p and [X,Y]=2i.

Validation: full Julia 1.13.3 tests passed, including 169 Gaussian/quadrature
and 77 photon-added Gaussian checks. Independent density-matrix cutoffs were
increased to 64/96 levels to achieve <1e-7 fourth-moment convergence; sampler
code itself has no cutoff.

## Separate Wplus track (complete)

Added generic positive-P convolution and exact Fock binomial–Gamma sampling.
Removed the finite window/asymptotic n=320 switch from the default positive-W
path; retained every exported Laguerre helper. An independent density identity
and symmetric first/fourth moments cover n=0,1,40,319,320,321,1000,10000.
Projection is explicitly Gaussian-smoothed W, never an ordinary-W initialiser.
Weighted ordinary W remains deferred as directed.

Validation: full Julia 1.13.3 tests passed, including 140 exact positive-W checks.

## Final validation (complete)

Enabled all Aqua checks, including previously disabled undefined-export and
compatibility checks. Parameters 0.12 reexports an undefined `@pack`; explicit
reexport of its seven implemented APIs (and module binding) preserves working
legacy names and removes that nonexistent symbol. Added Statistics compatibility
and CI coverage for Julia 1.10 and 1.13. Benchmarks and executable examples have completed; no default Julia channel or
installed runtime was changed.


Full final package tests pass on Julia 1.10.12 and 1.13.3, with every Aqua check
enabled. Statistical tests use explicit seeded RNGs; default-RNG API wrappers
retain argument/return tests. Signed Fock number variance uses its estimator
uncertainty. Documenter doctests and local HTML build pass without warnings.
All 23 examples run. Benchmark report and all 66 measurements are committed.
See VALIDATION.md for counts and benchmark/REPORT.md for variance/time comparisons.

## Resume / remaining limitations

No required implementation stage remains open. Review or release/merge is the
next human decision; no network publication was performed. The development
branch has separate commits for stages 0–6 and positive-W. Known limits:

- Ordinary Fock W remains explicitly approximate for every n>0.
- Generic superposition Gram preparation supports the documented coherent,
  Fock, product and squeezed component families. Passive/collective inputs
  expose pointwise overlaps but require further pairwise-Gram implementations
  to become components; unsupported cases fail informatively.
- Near-destructive general-phase sums can have low acceptance. Finite-precision
  overlap/Gram limits and configurable proposal exhaustion are explicit.
- The public passive constructor accepts explicit unitary matrices; structured
  operators can extend the existing mul! sampling contract.
- Photon addition to Gaussian states is centred, single-mode only.
- Weighted ordinary W, direct general Gaussian positive-P factorisation, GPU
  support and automatic differentiation are deferred/optional as directed.
- Initial estimator variance does not establish dynamical positive-P stability.
- Local tests passed; GitHub-hosted CI was configured but not triggered here.
