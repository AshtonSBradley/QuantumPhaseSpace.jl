# Codex directive: PhaseSpaceTools.jl sampling development

## Objective

Implement the work below in priority order in https://github.com/AshtonSBradley/PhaseSpaceTools.jl. The main aim is reliable, reusable quantum-state initialisation, with particular emphasis on exact Husimi-Q sampling and the corresponding positive-P initial conditions.

Work through the implementation, tests, documentation and benchmarks. Complete each stage before moving to the next dependent stage. Use reviewable commits on a dedicated development branch, following the repository's own instructions. Make routine implementation decisions autonomously. Record genuine blockers and continue independent work where possible.

The source review behind this directive inspected master at commit b180b50f0eddfed02d2f19728c6cbb5d15865b36, with Project.toml version 1.0.7. Reinspect the current checkout: some findings may already be fixed. The earlier review checked selected formulas independently in Python but did not run the Julia package tests. Establish an actual Julia baseline.

This directive authorises implementation work. Release and merge decisions are separate.

## Scientific constraints

1. Keep ordinary Wigner W, Husimi Q, Glauber P, positive-P and doubled positive-W distinct in the API, documentation and tests.
2. Positive-W development is wanted, but on its own track. Its samples are not interchangeable with ordinary W samples. Do not use a positive-W sampler as a fallback for an ordinary-W method or as a default initialiser for ordinary truncated-Wigner dynamics.
3. For an ordinary phase-plane sample, the second returned coordinate is its complex conjugate. In doubled representations, the coordinates need not be conjugates and are generally correlated. Do not describe them as statistically independent.
4. An exact Q initial-state sampler does not establish low sampling variance in positive-P dynamics. Measure initial estimator variance separately from sampler speed.
5. Classical mixtures and coherent superpositions are different operations. Preserve interference phases when sampling a superposition.
6. An exact sampler for a supplied finite Fock expansion is exact for that expansion. If an infinite state was truncated, report the cutoff and truncation error separately.
7. Keep the current approximate large-n Fock W sampler explicitly approximate for n > 0. Ordinary Wigner negativity cannot be removed by relabelling a doubled-space sampler.
8. Do not normalise individual stochastic field samples to a prescribed particle number. A fixed-number quantum state generally has fluctuating phase-space sample norms.
9. Preserve existing public entry points and return conventions unless a documented migration is necessary. Add generic machinery incrementally.

## Conventions

Use circular complex Gaussian noise

    ξ = (η₁ + iη₂)/√2,   ηⱼ ~ Normal(0,1),

so E[ξ] = E[ξ²] = 0 and E[|ξ|²] = 1. In Julia, randn(rng, ComplexF64) has this convention.

Throughout this directive, Gamma(k, θ) means shape k and scale θ. A uniformly random phase lies in [0, 2π).

For Gaussian-state covariance matrices use

    a = (x + ip)/√2,   [x,p] = i,   V_vac = I/2.

Document coordinate ordering, covariance definitions and squeezing conventions explicitly.

The established mode-by-mode Q-to-positive-P conversion is

    μ ~ Q,   γ ~ CN(0, I), independent,
    α = μ + γ,
    αplus = conj(μ) - conj(γ).

Use specialised positive-P samplers when they have lower variance, and retain this construction as an exact generic fallback.

---

## Priority 0 — establish the baseline and fix correctness

Read the repository instructions and inspect src, tests, documentation, dependencies and CI. Run the package tests before modifying the code and record the Julia version and starting commit.

Reproduce and fix the following issues if they remain present.

### 0.1 Squeezing convention

The reviewed Squeezed Q and W methods used opposite rotation phases. Define the state explicitly; the proposed convention is

    |ψ⟩ = D(β) S(ε)|0⟩,
    S(ε) = exp[(conj(ε)a² - ε a†²)/2],
    ε = r exp(iθ).

Then

    ⟨a⟩ = β,
    ⟨δa†δa⟩ = sinh²(r),
    ⟨δa²⟩ = -exp(iθ)sinh(r)cosh(r).

Make Q, W and positive-P agree on these moments after the appropriate ordering corrections. Test nonzero displacement and several genuinely complex squeezing phases. Occupation-only tests cannot detect the phase error.

### 0.2 Crescent Q sampling

At zero phase diffusion, Crescent should reduce to its declared squeezed state. The reviewed Q formula doubled the squeezed W noise, which does not give the Q distribution of the same state.

For this family, whose ordinary W function is nonnegative, an exact construction is:

    zW ~ W of the declared phase-diffused squeezed state,
    zQ = zW + ξ/√2.

Alternatively use a correctly constructed squeezed Q sample and apply the random phase rotation. Isotropic Gaussian convolution commutes with that rotation.

Define the phase-diffusion model precisely. Preserve equivalent legacy parameter behaviour where practical. Check the zero-diffusion reduction and nonzero-diffusion moments.

### 0.3 Bogoliubov Glauber P

The reviewed Glauber-P method directly transformed thermal P samples and missed the vacuum contribution. For the convention

    a = u b + conj(v)b†,   |u|² - |v|² = 1,

the target moments are

    ⟨a†a⟩ = (|u|² + |v|²)nbar + |v|²,
    ⟨a²⟩ = (2nbar + 1)u conj(v).

An ordinary nonnegative Gaussian P sampler exists only when V - I/2 is positive semidefinite. Sample that covariance when permitted; reject unsupported nonclassical cases with an informative error. Positive-P remains available. Do not patch this by adding arbitrary positive noise.

### 0.4 Vacuum, variance tests and documentation

- Make wigner(Fock(0), N) use the exact vacuum Gaussian.
- Preserve the large-n approximation for n > 0 and document its limitations.
- Fix the Fock positive-P number-variance test:

      Var(n) = E[(αplus)² α²] + E[αplus α] - E[αplus α]².

  Test the signed/complex quantities before extracting their physically real values. Do not use absolute values to conceal a systematic error.
- Add direct factorial-moment tests. For Fock n, ⟨a†²a²⟩ = n(n-1).
- Correct the positive-P independence wording and any ordering-index errors.
- Check normalisation constants in the phase-space derivations.

Completion: baseline issues explained; deterministic regression tests and seeded statistical tests pass; documentation states consistent conventions.

---

## Priority 1 — add the minimum reusable sampling infrastructure

Implement only the structure needed by subsequent stages; avoid a wholesale rewrite.

- Explicit RNG arguments, with backward-compatible default-RNG wrappers.
- A consistent internal multimode layout, preferably modes × samples for each coordinate array. Preserve the current public single-mode and existing two-mode interfaces through adapters.
- Reusable sampler objects where preparation is nontrivial: categorical tables, Gram matrices and covariance factors.
- In-place/batched paths where they materially reduce allocations.
- Validated dimensions, finite parameters, nonnegative occupations, integer photon counts and sample counts.
- Diagnostics for rejection samplers: proposals, acceptance rate and a configurable attempt limit. On exhaustion, raise an informative error; never return an incomplete batch as though it were complete.
- Capability dispatch that distinguishes pure-state overlap evaluation from mixed-state sampling.
- Stable complex coherent-state overlap evaluation, preferably log magnitude plus phase or an equivalent scaled representation. Avoid underflow in high Fock occupations and widely separated components.

Do not require GPU support or automatic differentiation in this stage. Keep Float64 defaults; introduce parametric numeric types where straightforward without delaying the sampling work.

Add the low-noise thermal positive-P specialisation using its ordinary P samples:

    α = β + √nbar ξ,
    αplus = conj(α).

Completion: old public calls still work; seeded calls are reproducible; multimode data flow and normalisation conventions are documented.

---

## Priority 2 — implement the inexpensive exact Q families

### 2.1 Displacement, rotation and mixtures

Add composable operations rather than duplicating state formulas.

- Displaced(state, β): z → z + β for Q.
- Rotated(state, θ), with R(θ) = exp(iθ n): z → exp(iθ)z.
- Mixture(states, probabilities): draw a component and sample it.
- ProductState: independent component Q samples, with an explicit mode layout.

For pure states, displacement must also transform the overlap phase:

    ⟨z|D(β)|ψ⟩
      = exp[(β conj(z) - conj(β) z)/2] ⟨z-β|ψ⟩.

A missing phase does not affect a displaced state's Q density by itself, but corrupts subsequent coherent superpositions.

### 2.2 Displaced Fock states and number mixtures

For D(β)|n⟩:

    X ~ Gamma(n+1, 1),
    z = β + √X exp(iθ).

For ρ = Σₙ pₙ |n⟩⟨n|, draw n from p, then use the Fock Q sampler. Support finite probability vectors first. An optional distribution-valued number law can be added without imposing an arbitrary cutoff.

### 2.3 Photon-added thermal states

Define

    ρ_add ∝ (a†)^m ρ_th a^m,

where m ≥ 0 is an integer and nbar ≥ 0 is the mean occupation of the original thermal state.

The exact Q sampler is

    X ~ Gamma(m+1, nbar+1),
    z = √X exp(iθ).

Validation:

    ⟨n⟩ = (m+1)nbar + m,
    E_Q[|z|²] = (m+1)(nbar+1),
    E_Q[|z|⁴] = (m+1)(m+2)(nbar+1)².

At nbar = 0 this reduces to the Fock-m Q sampler.

### 2.4 Photon-subtracted thermal states

Define

    ρ_sub ∝ a^m ρ_th (a†)^m.

For m > 0 require nbar > 0; conditioning on photon subtraction from the vacuum is undefined. Handle m = 0 as the original thermal state.

Sample its nonnegative Glauber P distribution:

    I ~ Gamma(m+1, nbar),
    aP = √I exp(iθ).

Then

    Q:          z = aP + ξ,
    ordinary W: z = aP + ξ/√2,
    positive-P: (aP, conj(aP)).

This family genuinely admits these ordinary P and W samplers.

Validation:

    ⟨n⟩ = (m+1)nbar,
    ⟨a†²a²⟩ = (m+1)(m+2)nbar².

Completion: analytic moment tests, parameter limits, examples and representation-support entries are present for every new family.

---

## Priority 3 — build the finite-superposition Q engine

This is the central generic extension.

Let the components |ψⱼ⟩ be normalised pure states, with samplable Q functions and evaluable complex overlaps

    Aⱼ(z) = ⟨z|ψⱼ⟩.

Define

    |ψ⟩ = Σⱼ cⱼ|ψⱼ⟩ / √Z,
    Z = Σⱼₖ conj(cⱼ)cₖ ⟨ψⱼ|ψₖ⟩,
    C = Σⱼ |cⱼ|.

Sample as follows:

1. Choose j with probability |cⱼ|/C.
2. Draw z from the Q function of |ψⱼ⟩.
3. Accept with probability

       |Σⱼ cⱼ Aⱼ(z)|²
       ------------------------------- .
       C Σⱼ |cⱼ| |Aⱼ(z)|²

Cauchy–Schwarz bounds this probability by one. The overall acceptance is Z/C².

Implement this directly for one or multiple modes. Use scaled amplitudes so common tiny factors cancel numerically. Validate the Gram normalisation and reject the zero state. Treat nearly cancelling components carefully, expose low acceptance and use stable state-specific formulas where needed.

Do not substitute component probabilities |cⱼ|² for the specified proposal weights without deriving a corresponding valid envelope. Do not substitute a classical mixture for the coherent superposition.

Initial supported components:

- Coherent states.
- Fock states and products of Fock states.
- Displaced versions of those pure states.
- Pure squeezed Gaussian states after their complex overlap convention is validated.

Test invariance under a common nonzero scaling of all coefficients, global phase, and permutation of components. Test interference-sensitive moments against independent finite-Fock calculations at manageable cutoffs.

### 3.1 Cat states

Add convenience constructors for even, odd and general relative-phase cats,

    |Cφ⟩ ∝ |α⟩ + exp(iφ)|-α⟩.

For the direct two-Gaussian proposal:

    z = ±α + ξ, with equiprobable signs,
    u = Re(conj(α)z),  v = Im(conj(α)z),
    p_accept = [1 + cos(2v+φ)/cosh(2u)]/2.

The overall acceptance is

    [1 + exp(-2|α|²)cosφ]/2.

Use stable sech evaluation. For even and odd cats, use exact parity formulas rather than relying on sin(π) being numerically zero.

The even cat has acceptance at least 1/2. The basic odd-cat proposal becomes inefficient as |α| → 0. Do not replace a finite-amplitude odd cat by Fock(1) silently.

An exact alternative proposal for a small odd cat is available:

    A = |α|² > 0,
    X ~ Gamma(2, 2),   z = √X exp(iθ),
    t = conj(α)z,
    p_accept = exp(-|z|²/2 - 2A) |sinh(t)/t|².

Use a stable analytic continuation sinh(t)/t → 1 at t = 0. Its mean acceptance is

    [1 - exp(-2A)] / [8 A exp(A)],

which tends to 1/4 as A → 0. Independently verify this bound and distribution before using it, then select between the two exact proposals using their expected acceptance rates. Use expm1 for small-A normalisation. The state at exactly α = 0 with odd parity is undefined unless the caller explicitly requests its Fock(1) limit.

For a general phase φ near destructive cancellation, retain diagnostics and document the efficiency limitation unless a separately verified proposal is added.

Cat validation, with d = exp(-2|α|²) and D = 1 + d cosφ:

    ⟨a⟩ = -i α d sinφ / D,
    ⟨a²⟩ = α²,
    ⟨n⟩ = |α|²(1 - d cosφ)/D.

For Q, E[|z|²] = ⟨n⟩ + 1.

### 3.2 Other immediate applications

Implement and document:

- Compass states and arbitrary finite coherent superpositions.
- Sparse single-mode Fock superpositions, including small bosonic-code examples.
- Sparse multimode Fock superpositions.
- NOON states, (|N,0⟩ + exp(iφ)|0,N⟩)/√2 for N ≥ 1.
- Photon-added coherent states.
- Squeezed cats once the squeezed overlaps are verified.

The generic engine gives exactly 1/2 acceptance for NOON states, independent of N. Test phase-sensitive coherence at small N; very high-order moment tests at large N may have large variance despite efficient generation of Q samples.

For photon-added coherent states use the exact finite expansion

    (a†)^m |β⟩
      = D(β) Σₖ₌₀ᵐ binomial(m,k) conj(β)^(m-k) √(k!) |k⟩.

This requires m+1 Fock components, not a truncation of the coherent state. For m = 1, the proposed superposition rejection scheme has acceptance

    (1 + |β|²)/(1 + |β|)² ≥ 1/2.

Completion: generic rejection proof documented; expected acceptance rates verified; no Fock truncation hidden in the analytic cat or photon-added coherent samplers; Q-to-positive-P conversion works for all new states.

---

## Priority 4 — add passive multimode transformations and fixed-number condensates

For a passive mode transformation a_out = U a_in with unitary U, Q samples transform as z_out = U z_in. Support an explicit matrix first, with a pathway for user-supplied structured operators.

This allows entangled states made by mixing product Fock, coherent and other supported inputs on beam splitters or interferometers. Include vacuum modes explicitly in the Q input.

Add a fixed-number condensate constructor for

    |N:u⟩ = (Σⱼ uⱼ aⱼ†)^N |0⟩ / √(N!),
    u†u = 1.

An exact O(M) Q sampler is

    X ~ Gamma(N+1, 1),
    θ ~ Uniform(0,2π),
    ξ ~ CN(0,I), independent,
    z = u √X exp(iθ) + ξ - u(u†ξ).

This samples the occupied collective mode as Fock N and the orthogonal modes as vacuum.

Validate

    E[z] = 0,
    E[z z†] = I + N u u†,

and selected higher moments against a small-mode Hilbert-space calculation. Include complex orbitals, N = 0, a single mode, and two-mode examples.

A specialised positive-P implementation can sample only the occupied collective Fock mode and use deterministic zero positive-P vacuum amplitudes in orthogonal modes. Derive and benchmark this against the generic Q conversion; the two valid representations can have different estimator variances.

Completion: mode convention, vacuum treatment and examples are documented; no per-trajectory particle-number renormalisation is introduced.

---

## Priority 5 — add exact loss and phase-diffusion wrappers

For any supported Q sampler, implement pure loss with 0 ≤ η ≤ 1:

    z_out = √η z_in + √(1-η) ξ.

For multiple modes, use independent environment noise and explicit mode-dependent transmissivities when supported.

This is an exact Q representation of the loss channel, including decoherence of superpositions. It is not equivalent to merely shrinking the coherent amplitudes in a pure cat constructor.

Add phase diffusion using

    θ ~ Normal(0, σ²),
    z_out = exp(iθ)z_in.

Distinguish independent mode phases from a common global phase. Define that choice in the API. A single-mode Q moment E[(z*)^p z^q] acquires the factor exp[-(p-q)²σ²/2].

Test:

- Identity and complete-loss limits.
- Coherent states remain coherent under pure loss.
- Output occupations scale by η under pure loss.
- A lossy Fock state agrees with the corresponding binomial number mixture.
- Lossy cats agree with an independent density-matrix calculation.
- Phase-diffusion factors agree for first, anomalous second and selected higher moments.

Completion: examples cover imperfect cat preparation, number states with loss, and phase-diffused squeezed states.

---

## Priority 6 — consolidate Gaussian sampling and add the next non-Gaussian extension

Introduce a general multimode Gaussian state representation with a displacement vector and real symmetric quadrature covariance V. Validate the uncertainty relation V + iΩ/2 ≥ 0, using numerical tolerances appropriate to the input precision.

Sampling covariance:

| Representation | Covariance |
|---|---|
| Ordinary W | V |
| Q | V + I/2 |
| Ordinary Glauber P | V - I/2, only when positive semidefinite |

Cache factors and support semidefinite classical P limits. Provide constructors for displaced squeezed thermal states, two-mode squeezed thermal states and multimode Bogoliubov inputs. Reuse existing state names where possible and fill the current two-mode-squeezed ordinary-W gap if still present.

Retain Q conversion for generic positive-P. A direct Gaussian positive-P factorisation is optional and should be introduced only with an independent derivation and estimator-variance comparison.

A useful subsequent non-Gaussian extension is photon addition to an undisplaced Gaussian state:

    Q_added(z) ∝ |z|^(2m) Q_Gaussian(z).

After diagonalising its two-dimensional Q covariance, expand (x₁²+x₂²)^m into nonnegative monomials. This gives a finite mixture of squared-coordinate gamma distributions.

For one added photon, if the principal-coordinate variances are v₁ and v₂:

1. Choose coordinate j with probability vⱼ/(v₁+v₂).
2. For that coordinate, draw xⱼ² ~ Gamma(3/2, 2vⱼ), with equiprobable sign.
3. Draw the other coordinate from its original zero-mean Gaussian.
4. Rotate back to the original axes.

Independently verify this construction before implementation. It covers photon-added squeezed vacuum and squeezed thermal states. Do not apply the centred formula unchanged to displaced Gaussian states.

Completion: covariance/uncertainty tests, cross-mode and anomalous moments, Gaussian limits and low-dimensional Hilbert-space checks pass.

---

## Separate track Wplus — exact doubled positive-W sampling

Keep this work in separately reviewable commits, clearly labelled as positive-W. It must not block completion of the main Q/positive-P stages.

### Wplus 1: generic convolution

Given valid positive-P samples, draw independent η = ξ/√2 per mode:

    αW = αP + η,
    αWplus = αPplus + conj(η).

This samples a doubled positive-W representation. It does not create ordinary W samples and does not establish a practical positive-W dynamical scheme.

### Wplus 2: direct exact Fock sampler

Replace or supplement the existing finite-window Laguerre rejection sampler with

    K ~ Binomial(n, 2/3),
    X | K ~ Gamma(K+1, 3/2),
    μ = √X exp(iθ),
    γ ~ CN(0,1), independent,
    α = μ + γ,
    αplus = conj(μ) - conj(γ).

The binomial-gamma mixture follows by expanding the implemented Laguerre radial density. Verify the identity independently, then benchmark it against the generic positive-P-to-positive-W convolution.

Cover n = 0, 1, moderate n, values around the former n = 320 switch, and large n. Check symmetrically ordered moments:

    E[α αplus] = n + 1/2,
    E[α² (αplus)²] = n² + n + 1/2.

Remove the rejection window and asymptotic switch from the default exact Fock positive-W path once validated. Preserve or deprecate exported Laguerre helpers deliberately; do not remove public symbols incidentally.

### Wplus 3: document the distinction from ordinary W

For the canonical Q-derived construction,

    μ = [α + conj(αplus)]/2

has density

    Fμ = Q * G_(1/2) = W * G_1,
    G_v(z) = exp(-|z|²/v)/(πv).

Thus the natural projection produces a Gaussian-smoothed W, not W. Recovering W by deconvolution is a separate, noise-amplifying inverse problem.

No implementation should silently project doubled samples and pass them to an ordinary-W solver.

---

## Deferred research track — weighted ordinary W

Keep compact complex-P to weighted ordinary-W sampling outside the initial production scope. Study Drummond and Opanchuk (2020), Sec. V C, before proposing an implementation.

This route produces ordinary complex phase-plane coordinates with signed/complex weights. It is distinct from doubled positive-W and needs:

- A return type explicitly containing weights.
- A correct weighted-estimator interface.
- Convergence and weight-variance diagnostics.
- Benchmarks against direct signed-W sampling and existing approximations.
- A compact starting representation; arbitrary Q-derived positive-P samples can produce unusable weight tails.
- Separate reporting of initial-state sampling error and subsequent dynamical approximation error.

Do not claim a general unweighted positive-W-to-W reduction.

---

## Validation and benchmarking requirements

Use independent mathematical targets rather than tests that merely duplicate implementation.

1. Retain deterministic tests for dimensions, parameter validity, constructor reductions and normalisation conventions.
2. Use explicit seeded RNGs for statistical tests.
3. Test complex first moments, anomalous moments, occupations, factorial moments and cross-mode correlations.
4. Choose tolerances from estimator uncertainty. Keep large, noisy high-order studies in benchmarks rather than making CI flaky.
5. Compare a modest selection of new non-Gaussian states against independent finite-Fock density-matrix calculations, with controlled cutoff convergence.
6. Test zero samples, vacuum and zero-operation limits, complex amplitudes, nonzero phases, large occupations and near cancellation.
7. Report runtime, allocations, rejection acceptance, preparation cost, and error per unit time for selected physical estimators.
8. Benchmark representative small, moderate and large batches. Use the same physical observable and representation when making efficiency comparisons.
9. Distinguish theoretical exactness, floating-point limitations, explicit approximation and Hilbert-space truncation.
10. Run package tests and the relevant documentation/quality checks. Reconcile the declared Julia compatibility with the versions actually tested.

No production sampler should rely on an undocumented finite spatial window, silently discard interference, or conceal a failed normalisation.

## Required deliverables

- Working implementation on a development branch, with separate commits for the main stages and positive-W track.
- Tests and documentation for every public feature.
- A state × representation support table marking exact, approximate and unsupported combinations.
- Examples for a cat, photon-added thermal state, photon-added coherent state, NOON state, fixed-N condensate and lossy cat.
- A benchmark report including estimator variance, not only generation speed.
- A concise change report: completed stages, validation results, API decisions, remaining limitations and deferred research.
- A progress record in the repository so a later Codex session can resume without repeating finished work.

## Sources and derivation status

Use the following as primary references. Recheck normalisation and phase conventions against the implementation rather than copying formulas blindly.

- Repository snapshot reviewed:
  https://github.com/AshtonSBradley/PhaseSpaceTools.jl/tree/b180b50f0eddfed02d2f19728c6cbb5d15865b36
- Package phase-space derivations:
  https://ashtonsbradley.github.io/PhaseSpaceTools.jl/dev/phasespace/
- Olsen and Bradley, Numerical representation of quantum states in the positive-P and Wigner representations, Optics Communications 282, 3924–3929 (2009):
  https://doi.org/10.1016/j.optcom.2009.06.033
- Olsen, Lewis-Swan and Bradley, thermal-state erratum (2016):
  https://doi.org/10.1016/j.optcom.2016.02.068
- Weedbrook et al., Gaussian quantum information, Rev. Mod. Phys. 84, 621 (2012):
  https://doi.org/10.1103/RevModPhys.84.621
- Barnett et al., Statistics of photon-subtracted and photon-added states, Phys. Rev. A 98, 013809 (2018):
  https://arxiv.org/abs/1806.09475
- Drummond and Opanchuk, Initial states for quantum field simulations in phase space, Phys. Rev. Research 2, 033304 (2020), especially Sec. V C:
  https://doi.org/10.1103/PhysRevResearch.2.033304
- Meena, Kumar and Banerjee, Husimi phase distribution in non-Gaussian operations; Q loss-channel convolution, Eq. (35) in the preprint:
  https://arxiv.org/abs/2503.22298
- Drummond, Higher-order stochastic differential equations and the positive Wigner function, Phys. Rev. A 96, 062104 (2017), for the separate dynamical research background:
  https://doi.org/10.1103/PhysRevA.96.062104

The finite-superposition rejection envelope, the specialised cat acceptance formulas, the small-odd-cat alternative proposal, the binomial-gamma positive-W sampler, and the explicit collective-orbital sampler are constructions derived during the planning discussion. Treat them as formulas to verify independently, not as claims of established novelty. Preserve the derivations and validation evidence with the implementation.

