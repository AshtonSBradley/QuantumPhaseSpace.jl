# Exact Q sampling and quantum superpositions

Start with the [cat-state Q tutorial](@ref cat-q-tutorial) for runnable examples
and sampled even/odd/mixture comparisons. This page develops the rejection
bounds, overlap conventions and extensions behind those examples.

## Pure overlaps and coherent sums

`Superposition(states,c)` forms a coherent sum, not a mixture. All components
must supply normalised pure overlaps `Aⱼ(z)=⟨z|ψⱼ⟩`. Coherent, Fock and pure
squeezed states, their products, displacements and rotations are supported.
For `D(β)S(ϵ)|0⟩` we use
`S(ϵ)=exp[(conj(ϵ)a²-ϵa†²)/2]`, `ϵ=r exp(iθ)`, so the centred anomalous
moment is `-exp(iθ)sinh(r)cosh(r)` in every representation.

Writing `C=Σ|cⱼ|` and `Z=c†Gc`, the proposal is
`g(z)=Σ|cⱼ| |Aⱼ(z)|²/(C πᴹ)`. Cauchy–Schwarz gives

```math
\left|\sum_j c_j A_j\right|^2
\le C\sum_j |c_j||A_j|^2.
```

Acceptance is the ratio of the left and right sides. Integrating accepted
proposals gives exactly the target Q density times `Z/C²`; thus
`state.expected_acceptance` predicts efficiency. Coefficients are rescaled to
`C=1` on construction. Common amplitude log scales cancel in rejection tests.
The cached Gram matrix comes from analytic Bargmann Gaussian integration and
normalised polynomial recurrences, evaluated at higher precision during
preparation. This uses no Fock cutoff. Near-zero states with extremely small
acceptance remain inefficient and can exhaust the configurable attempt limit.
Pure-overlap capability alone does not guarantee that every pairwise Gram
integral is implemented. Passive and collective-orbital states have pointwise
overlaps but are not currently supported as generic superposition components.
For operations wrapping an entire sum, expand the operation componentwise before
nesting it in another sum. Unsupported Gram preparation raises an informative
error. Finite precision can make extremely ill-conditioned Gram problems unresolved;
such cases raise an error rather than inventing a normalisation.

Displacement retains the phase
`(β conj(z)-conj(β)z)/2` in `logoverlap`. Omitting it changes a subsequent
superposition, even though it leaves an individual Q density unchanged.

## [Cats, including the small odd-cat limit](@id small-odd-cat)

`Cat(α,φ)`, `EvenCat(α)` and `OddCat(α)` use the two coherent components
analytically, without truncation. Their direct two-Gaussian proposal has
acceptance `[1+exp(-2|α|²)cosφ]/2`. Exact parity arithmetic avoids a spurious
`sin(π)` phase. Even cats have acceptance at least one half. General phases
close to destructive cancellation may have very low acceptance.

For an odd cat, let `A=|α|²>0`, `t=conj(α)z`. The alternate proposal is
`q(z)=|z|² exp(-|z|²/2)/(4π)`, obtained from `|z|²~Gamma(2,2)` and uniform phase.
Its rejection probability is

```math
p(z)=e^{-|z|^2/2-2A}\left|\frac{\sinh t}{t}\right|^2.
```

To verify the envelope independently, use
`sinh(t)/t=∫₀¹cosh(st)ds`, so `|sinh(t)/t|≤exp(|t|)`. Hence
`log p≤-|z|²/2-2A+2√A|z|=-(|z|-2√A)²/2≤0`.
Multiplying `q p` gives
`exp(-|z|²-2A)|sinh t|²/(4πA)`. Comparing with the normalised odd-cat Q
function yields mean acceptance `(1-exp(-2A))exp(-A)/(8A)`, tending to 1/4.
The implementation selects whichever exact proposal has higher expected
acceptance and evaluates `sinh(t)/t` by analytic continuation near zero.
`OddCat(0)` raises an error; `Fock(1)` is an explicit limiting choice.

## Finite expansions and multimode applications

`FockSuperposition(n,c)` accepts sparse single-mode integer occupations or a
modes × components integer matrix. It is exact for that expansion. For an
approximation to an infinite state, report the supplied cutoff and omitted
probability independently. `NOON(N;phase)` uses two orthogonal product Fock
components and has exactly 1/2 acceptance at every N. Efficient sample generation
does not imply low variance for high-order coherence estimators.

`PhotonAddedCoherent(m,β)` uses the exact identity

```math
(a^\dagger)^m|\beta\rangle = D(\beta)\sum_{k=0}^m
 {m\choose k}(\beta^*)^{m-k}\sqrt{k!}|k\rangle.
```

There are exactly `m+1` components and no coherent-state cutoff. For one added
photon the acceptance is `(1+|β|²)/(1+|β|)²≥1/2`.
`SqueezedCat` uses displaced pure squeezed components with validated complex
overlaps. All these states inherit the exact generic Q-to-positive-P conversion;
its estimator variance must still be measured for the observable of interest.

```@docs
Superposition
CoherentSuperposition
Compass
FockSuperposition
NOON
PhotonAddedCoherent
SqueezedCat
Cat
EvenCat
OddCat
```

## Independent validation and references

The test suite constructs annihilation matrices, displacement/squeezing matrix
exponentials, parity Fock expansions and operator moments independently of the
samplers. The odd-cat accepted density is also compared pointwise with a
separately normalised odd Fock series. Production sampling does not use these
finite matrices. Tests use modest parameters with negligible tails and explicitly
check normalisation; convergence checks compare multiple cutoffs.

References for representation conventions and state statistics:

- [Olsen and Bradley (2009)](https://arxiv.org/abs/0906.3905), with the
  [thermal-state erratum (2016)](https://doi.org/10.1016/j.optcom.2016.02.068).
- [Weedbrook et al. (2012)](https://arxiv.org/abs/1110.3234); our quadratures
  use `[x,p]=i` and vacuum covariance `I/2`, rather than that review's `2i`.
- [Barnett et al. (2018)](https://arxiv.org/abs/1806.09475), for factorial
  statistics of photon addition/subtraction and loss.
- [Meena, Kumar and Banerjee (2025), Eq. (35)](https://arxiv.org/html/2503.22298v1),
  for the Q loss-channel convolution.
- [Drummond and Opanchuk (2020)](https://doi.org/10.1103/PhysRevResearch.2.033304),
  for broader sampling and weighted methods. Weighted ordinary W is deferred.

The rejection envelopes, small-odd proposal, binomial–Gamma positive-W law and
collective-orbital construction are derived and verified here; no novelty claim
is implied by their inclusion.

## Passive transformations and fixed-number condensates

`Passive(input,U)` uses `a_out=U*a_in` with a verified square unitary matrix.
Supply vacuum modes explicitly, e.g. `ProductState(Fock(3),Fock(0))`. Each
column of the Q matrix transforms by U. Internally `mul!` separates the mode
operator from sample generation, allowing future structured-operator adapters.

For `FixedNumberCondensate(N,u)`, `u†u=1`, draw a collective Fock Q amplitude
`c=√X exp(iθ)`, `X~Gamma(N+1,1)`, and independent `ξ~CN(0,I)`. Then
`z=u*c+(I-u*u†)ξ`. Orthogonality and `E|c|²=N+1` give
`E[z z†]=I+N*u*u†`. No trajectory is rescaled to a fixed norm.

The specialised positive-P draws one collective Fock positive-P pair and sets
`α=u*a`, `αplus=conj(u)*ap`; orthogonal modes have deterministic zero positive-P
vacuum. Normal ordering makes every orthogonal-vacuum contribution vanish,
proving the construction. `method=:generic` uses independent noise in every
mode through Q conversion. Both are exact, with potentially different estimator
variances; benchmark them using the same observable.

```@docs
Passive
FixedNumberCondensate
```

## Preparation channels

`Loss(state,η)` implements `z_out=√η*z_in+√(1-η)*ξ` independently per mode,
with scalar or per-mode transmissivity and `E|ξ|²=1`. This is the exact Q
channel, including mixed-state decoherence. Shrinking a pure cat's amplitude
alone would omit the environment's information about its components.

`PhaseDiffused(state,σ;common=true)` applies one Gaussian phase of variance σ²
to all modes. `common=false` draws independent phases with scalar or per-mode
standard deviations. In one mode the Q moment `(z*)^p z^q` acquires
`exp[-(p-q)²σ²/2]`. `Crescent(β,ϵ,q)` retains its legacy signed q parameter,
with physical standard deviation `abs(q)`, and equals phase diffusion of the
declared squeezed state. These channels sample mixed states and do not expose
pure-state overlaps for coherent superpositions.

```@docs
Loss
PhaseDiffused
```

## General Gaussian states and quadrature units

`Gaussian(β,V)` orders quadratures `(x₁,p₁,x₂,p₂,...)`, with
`x=(a+a†)/√2`, `p=(a-a†)/(i√2)` and `[x,p]=i`. V is the symmetrised
**centred** covariance. Vacuum V is I/2. The equally valid convention
`X=a+a†=√2x`, `Y=(a-a†)/i=√2p` has `[X,Y]=2i`, vacuum covariance I,
and `a=(X+iY)/2`. Convert an `(X,Y)` covariance by `V=VXY/2` before passing
it to this API; the complex displacement β is unchanged. Regression tests check
both commutators, reconstructions of a, and vacuum sample variances.

The W covariance is V, Q covariance V+I/2, and nonnegative ordinary P covariance
V-I/2. The uncertainty test is `V+iΩ/2≥0`, where Ω has 2×2 blocks
`[0 1;-1 0]`. Eigenvalue tolerances scale with the input precision and covariance
norm. Covariance factors are prepared once, including zero eigenvalues in
classical P limits. Coherent/thermal/squeezed legacy fast paths retain their
public calls; tests compare their conventions to the general covariance model.

`SqueezedThermal` follows the single-mode negative anomalous-moment convention.
`TwoModeSqueezedThermal` retains the legacy `SqueezedTwoMode(r,ϕ)` convention
`⟨a₁a₂⟩=exp(2iϕ)sinh(r)cosh(r)` for vacuum input. `GaussianBogoliubov` uses
`a=U*b+conj(V)*b†` and checks both canonical identities. Every Gaussian state
has Q-derived positive-P; ordinary P rejects nonclassical covariances.

For `PhotonAddedGaussian` the state must be centred and single mode. Diagonalise
its quadrature Q covariance with principal variances v₁,v₂. Expanding
`(x₁²+x₂²)^m` gives nonnegative weights
`binomial(m,k) E[x₁^(2k)] E[x₂^(2(m-k))]`. Each component draws independent
signed squared coordinates from `Gamma(k+1/2,2v₁)` and
`Gamma(m-k+1/2,2v₂)`, then rotates back and divides by √2. This follows directly
by multiplying each normal density by its even monomial. For m=1, choosing the
inflated coordinate has probability vⱼ/(v₁+v₂), with shape 3/2 on that coordinate
and 1/2 on the other. The common factor 2^-m from `|z|^(2m)` cancels in
normalisation. Tests independently apply creation operators to finite density
matrices and verify cutoff convergence. The centred formula is not applied to
displaced inputs.

```@docs
Gaussian
SqueezedThermal
TwoModeSqueezedThermal
GaussianBogoliubov
PhotonAddedGaussian
```

## Separate positive-W track

Generic `positiveW` convolves a valid positive-P pair with independent η=ξ/√2:
`αW=αP+η`, `αWplus=αPplus+conj(η)`. This is an exact doubled representation.
It supplies initial sampling, not a practical positive-W dynamical scheme.

For Fock states set `K~Binomial(n,2/3)`, `X|K~Gamma(K+1,3/2)` and
`μ=√X exp(iθ)`, then add/subtract independent unit circular γ in the canonical
pair. To verify the mixture, expand
`L_n(-4X/3)=Σ binomial(n,k)(4X/3)^k/k!`. Transforming the old radial law by
`X=r²` yields `2 exp(-2X/3)L_n(-4X/3)/3^(n+1)`. Termwise comparison with
Gamma(k+1,3/2) gives weights `binomial(n,k)(2/3)^k(1/3)^(n-k)`, summing to one.
There is no finite rejection window or change of approximation at n=320.
The exported Laguerre helpers and generic finite-window `reject` remain available
for compatibility, and are not used by the default exact Fock positive-W path.

For the canonical Q-derived construction the projection
`μ=(α+conj(αplus))/2` has density `Q*G_(1/2)=W*G_1`, where
`G_v(z)=exp(-|z|²/v)/(πv)`. It is a smoothed W, not ordinary W.
Recovering W requires a separate noise-amplifying deconvolution. No projection
is used to initialise ordinary-W dynamics. Weighted ordinary W, requiring
explicit signed/complex weights and variance/convergence diagnostics, remains
a deferred research track; arbitrary Q-derived positive-P may have unusable
weight tails in such a construction.
