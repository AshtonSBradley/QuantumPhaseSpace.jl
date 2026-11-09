# Quantum-state sampling change report

Implemented on `codex/exact-quantum-sampling`, starting at
`b180b50f0eddfed02d2f19728c6cbb5d15865b36` (v1.0.7). Release version, default
Julia channel and existing Julia installations are unchanged. No merge, push,
release or publication was performed.

## Completed work

- **0 — correctness:** consistent complex squeezing phases; exact Crescent Q;
  classicality-checked Bogoliubov P covariance; exact vacuum W; signed complex
  factorial/number-variance regressions; corrected ordering and normalisations.
- **1 — infrastructure:** explicit RNGs and compatible default wrappers;
  modes × samples Q matrices; in-place entry point; parameter/count checks;
  bounded rejection diagnostics; pure-overlap capability; stable logarithmic
  overlaps; low-noise thermal positive-P.
- **2 — exact families:** displacement with its overlap phase, rotations,
  mixtures/products, finite number mixtures, photon-added thermal Q, and
  photon-subtracted thermal P/W/Q/positive-P.
- **3 — coherent sums:** analytic Gram preparation, scaled rejection with |c|
  proposal weights; coherent/Fock/product/displaced/squeezed components; stable
  even/odd/general cats, compass states, sparse Fock codes, multimode NOON,
  exact photon-added coherent expansions and squeezed cats.
- **4 — multimode:** explicit passive unitaries and an O(M) fixed-number
  condensate Q sampler; specialised occupied-mode positive-P. Samples are
  never renormalised to a fixed norm or particle number.
- **5 — channels:** exact pure-loss Q convolution and Gaussian phase diffusion,
  with explicit common versus independent mode phases and per-mode parameters.
- **6 — Gaussian:** uncertainty validation and cached covariance factors;
  semidefinite classical P limits; squeezed thermal, two-mode and Bogoliubov
  constructors; missing two-mode W; exact centred photon-added Gaussian Q.
- **Wplus — separate commits:** generic positive-P convolution and exact Fock
  binomial–Gamma sampling, without a finite window or n=320 asymptotic switch.
  Existing Laguerre helpers remain exported. Doubled positive-W is kept distinct
  from ordinary W in code, support tables and derivations.

## API and conventions

Existing public sampling calls and single-mode vector returns remain. The
legacy two-mode squeezed methods retain four-vector returns. New multimode
states return a pair of matrices; `sampleQ` always returns one modes × samples
matrix. Pure and mixed sampling capabilities remain distinct. Some pointwise
pure overlaps (passive/collective states) do not yet have pairwise Gram machinery;
unsupported superposition components fail with a clear error.

The Gaussian API uses `a=(x+ip)/√2`, `[x,p]=i`, vacuum V=I/2. In the alternative
convention, X=√2x,Y=√2p, `[X,Y]=2i`, vacuum VXY=I: supply `VXY/2`. Explicit
commutator, reconstruction and vacuum covariance regressions pass. Single-mode
squeezing is `D(β)S(ϵ)` with a negative centred anomalous moment; the historical
two-mode phase is separately documented.

Julia compatibility is now explicitly 1.10+, tested locally on 1.10.12 and
1.13.3. Full Aqua checks are enabled; the pre-existing undefined Parameters
`@pack` re-export was removed while preserving implemented exports.

## Validation and benchmarks

Each final package run passes **1,572 assertions** on both tested Julia versions,
including full Aqua. Tests check complex first/anomalous/factorial moments,
cross-mode correlations, acceptance laws, exact limits, seeded reproducibility,
parameter validity and zero batches. Independent finite-Fock exponentials,
Kraus density matrices and cutoff-convergence checks validate selected states.
The added-Gaussian comparison uses 64/96 levels; production sampling introduces
no such cutoff. Documenter doctests and local HTML build pass without warnings.
All 23 runnable Q/positive-P examples produce finite complete batches.

`benchmark/REPORT.md` and `benchmark/results.tsv` contain 66 warmed measurements
at N=100,10,000,100,000, including preparation, allocations, acceptance and
estimator variance/time. All occupation means are within 2.45 measured standard
errors of their analytic targets. The 16-mode condensate specialisation reduces
mode-1 occupation variance by about 45× relative to generic Q conversion. Direct
Fock positive-W allocates much less memory; convolution is faster in these local
runs. These comparisons concern initial sampling, not dynamical stability.

## Limits and deferred work

Ordinary Fock W remains an explicit n>0 large-occupation approximation. Finite
Fock sums are exact only for the expansion supplied; external truncation error
must be stated separately. Near-cancelling generic sums can be inefficient or
numerically unresolved. Passive structured operators have an internal mul!
pathway but no public constructor yet. Gaussian photon addition is centred and
single-mode only. Weighted ordinary W, direct general Gaussian positive-P
factorisation, GPU/AD support and positive-W dynamics remain deferred or optional
as directed. No positive-W projection is used as an ordinary-W initialiser.

The archived user directive is `SAMPLING_DIRECTIVE.md`; `PROGRESS.md` records the
stage history, and `VALIDATION.md` records final check totals. No required stage
remains open. GitHub-hosted CI is configured for 1.10/1.13 but was not triggered
from this local branch.
