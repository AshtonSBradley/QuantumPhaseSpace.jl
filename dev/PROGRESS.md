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

## Remaining

Stages 1–6, separate Wplus track, examples, benchmarks and final validation.
No release, merge or publication is authorised.

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
