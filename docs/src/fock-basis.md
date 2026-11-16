# [Fock-basis examples](@id fock-basis-examples)

Use [`FockSuperposition`](@ref) for amplitudes in a finite number-state
expansion, and [`NumberMixture`](@ref) for probabilities of individual number
states. Both can be sampled in Husimi Q, positive-P and doubled positive-W.
The samplers are exact for the state supplied. For a projection of an
infinite-dimensional state, check the omitted probability and the cutoff
error in the observables of interest separately.

## A phase-sensitive number-state coherence

For $(|0\rangle+e^{i\phi}|3\rangle)/\sqrt{2}$, the radial Q distribution is
the same as for an equal mixture of $|0\rangle$ and $|3\rangle$. The
angular distribution retains the relative phase, and
$\langle a^3\rangle=e^{i\phi}\sqrt{6}/2$. The incoherent mixture has
$\langle a^3\rangle=0$. Powers of $a$ alone have the same symbol in Q
and positive-P, so both ensembles can estimate this moment directly.
The standard errors below refer to the real and imaginary components.

```@example fock_superposition
using QuantumPhaseSpace, Random, Statistics

rng = MersenneTwister(421)
phase = π/3
state = FockSuperposition([0, 3], [1, cis(phase)])
N = 50_000

Z, stats = sampleQ(rng, state, N; diagnostics=true)
α, _ = positiveP(rng, state, N)

function summarize(values)
    return (; estimate=mean(values),
            se_real=std(real.(values))/sqrt(length(values)),
            se_imag=std(imag.(values))/sqrt(length(values)))
end

(; Q=summarize(vec(Z).^3), positiveP=summarize(α.^3),
   exact=sqrt(6)*cis(phase)/2, acceptance=acceptance_rate(stats))
```

The expected acceptance is $1/2$ for these two equal-magnitude,
orthogonal components. See the [rejection rule and normalization](@ref finite-fock-superpositions)
for general coefficients.

## Kerr evolution of a coherent input

For $H/\hbar=\chi\hat n(\hat n-1)/2$ and $\tau=\chi t$, the coherent
input $|\beta\rangle$ has coefficients

```math
c_n(\tau)=e^{-|\beta|^2/2}\frac{\beta^n}{\sqrt{n!}}
           e^{-i\tau n(n-1)/2}.
```

The phases evolve while the number probabilities stay Poisson with mean
$\nu=|\beta|^2$. With $\beta=2e^{0.35i}$, retaining $0\le n\le27$
discards $5.02\times10^{-15}$ of the probability at every time. Build a
fresh ensemble from the coefficients at the time of interest. At
$\tau=2\pi/3$, $\langle a^3\rangle=\beta^3$, while
$|\langle a\rangle/\beta|=e^{-3\nu/2}$.

```@example kerr_snapshot
using QuantumPhaseSpace, Random, Statistics

rng = MersenneTwister(422)
β = 2cis(0.35)
τ = 2π/3
n = 0:27
coherent = exp(-abs2(β)/2) .* cumprod([1; β ./ sqrt.(1:last(n))])
coefficients = coherent .* cis.(-τ .* n .* (n .- 1) ./ 2)
state = FockSuperposition(n, coefficients)

Z, stats = sampleQ(rng, state, 50_000; diagnostics=true)
z = vec(Z)
α, _ = positiveP(rng, state, length(z))

(; Q_third=mean((z ./ β).^3), positiveP_third=mean((α ./ β).^3),
   exact_third=1, acceptance=acceptance_rate(stats))
```

The acceptance stays constant in time because it depends on the
coefficient magnitudes. The finite sum gives the state at a specified
time. It can also serve as a reference for testing stochastic integration
of the Kerr equations.

From the package checkout, run `julia --project=. examples/kerr_evolution.jl`
to compare eight normally ordered moments at four times. The example
checks the projected state against the full coherent-input solution,

```math
\langle(a^\dagger)^p a^q\rangle_\tau
=(\beta^*)^p\beta^q
 \exp\!\left(\frac{i\tau}{2}[p(p-1)-q(q-1)]\right)
 \exp\!\left(\nu[e^{i\tau(p-q)}-1]\right).
```

It compares the Q and positive-P samples with ladder-matrix moments in
the projected state and reports the sampling uncertainty and cutoff
error separately.

## Canonical thermal state in the number basis

For a number-diagonal Hamiltonian with energies $E_n$, the Gibbs
probabilities are $p_n\propto e^{-\beta_T(E_n-E_{\min})}$, where
$\beta_T$ is inverse temperature in the units of $H$. The energy
shift cancels on normalization. For a Kerr oscillator with
$E_n=\omega n+\chi n(n-1)/2$, the weights are generally distinct from
the geometric probabilities of [`Thermal`](@ref), which describes a harmonic
mode.

```@example canonical_thermal
using QuantumPhaseSpace, Random, Statistics

rng = MersenneTwister(423)
n = 0:35
βT = ω = χ = 1.0
energies = @. ω*n + χ*n*(n-1)/2
weights = exp.(-βT .* (energies .- minimum(energies)))
p = weights ./ sum(weights)
state = NumberMixture(p)

N = 50_000
z = vec(sampleQ(rng, state, N))
α, αplus = positiveP(rng, state, N)
R = abs2.(z)

(; Q=mean(R.^2 .- 4R .+ 2), positiveP=mean((αplus .* α).^2),
   exact=sum(p .* n .* (n .- 1)))
```

The observable here is $\langle\hat n(\hat n-1)\rangle$. The Q estimator
includes the antinormal-ordering correction. For these parameters, the
omitted Gibbs probability above $n=35$ is less than $5\times10^{-290}$.
The bound follows by bounding the omitted weights by a geometric
series with ratio $e^{-\beta_T(\omega+36\chi)}$.

For a finite Hamiltonian with off-diagonal matrix elements, diagonalize
it, form each occupied
eigenvector as a [`FockSuperposition`](@ref), and combine those states with
[`Mixture`](@ref) using the Gibbs probabilities. The coherences inside each
eigenstate then enter its Q sampler.

From the package checkout, run `julia --project=. examples/canonical_thermal.jl`
to compare the canonical Kerr, harmonic and displaced harmonic
oscillators with their reference moments.
