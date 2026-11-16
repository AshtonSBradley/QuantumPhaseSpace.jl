# State and representation support

This table covers bosonic states. See [fermionic truncated Wigner](@ref fermionic-twa)
for the Slater and free-thermal initializers.

**E**: exact probability law, subject to floating-point accuracy. **A**: explicit
large-occupation approximation. **C**: exact only for classical covariance
`V-I/2≥0`, otherwise an informative error. **—**: no public method supplied.
Positive-W is a doubled representation, not an ordinary phase-plane Wigner law.

| State / operation | Ordinary P | Q | Ordinary W | Positive-P | Positive-W |
|:--|:--:|:--:|:--:|:--:|:--:|
| Coherent | E | E | E | E | E |
| Fock(0) | — | E | E | E | E |
| Fock(n>0) | — | E | A | E | E |
| Thermal | E | E | E | E | E |
| Squeezed | — | E | E | E | E |
| SqueezedTwoMode | — | E | E | E | E |
| Crescent | — | E | E | E | E |
| Bogoliubov | C | E | E | E | E |
| Gaussian / SqueezedThermal / TwoModeSqueezedThermal / GaussianBogoliubov | C | E | E | E | E |
| Displaced / Rotated | — | E | — | E | E |
| Mixture / NumberMixture | — | E | — | E | E |
| ProductState | — | E | — | E | E |
| PhotonAddedThermal | — | E | — | E | E |
| PhotonSubtractedThermal | E | E | E | E | E |
| Superposition / CoherentSuperposition / FockSuperposition | — | E | — | E | E |
| Cat / EvenCat / OddCat / Compass / SqueezedCat | — | E | — | E | E |
| PhotonAddedCoherent | — | E | — | E | E |
| NOON | — | E | — | E | E |
| Passive | — | E | — | E | E |
| FixedNumberCondensate | — | E | — | E | E |
| Loss / PhaseDiffused | — | E | — | E | E |
| PhotonAddedGaussian (centred, single mode) | — | E | — | E | E |

Wrappers require a supported input Q sampler. Convenience constructors may return
one of the underlying state types. A dash records implemented API support,
even when an exceptional parameter limit would admit an ordinary P or W law.
Use `Gaussian(existing_state)` for available classical Gaussian P limits.
`SqueezedThermal` follows `D(β)S(ϵ)` ordering. The two-mode phase convention is
specified separately in its docstring.

For n>0 the Fock W sampler is the pre-existing radial Gaussian approximation,
most useful for n≫1; it cannot reproduce a negative Wigner function as an exact
unweighted law. Its n=0 path is exact. No positive-W samples are projected and
passed off as ordinary W, or used as ordinary truncated-Wigner initial data.

Ordinary P/Q/W return conjugate coordinate pairs. Positive-P and positive-W
coordinates generally are correlated and need not be conjugates. Specialised
coherent/thermal/photon-subtracted positive-P paths and the collective condensate
path reduce initial estimator noise when applicable. Exact Q generation alone
makes no promise about positive-P estimator variance or later dynamics.

Finite Fock expansions are exact for the coefficients supplied. Truncation of an
infinite state is the caller's explicit approximation: record its cutoff and
omitted tail probability separately. The analytic cat, photon-added coherent,
and collective-orbital samplers introduce no Hilbert-space cutoff.

## Examples

Start with the [cat-state Q tutorial](@ref cat-q-tutorial) for a visual comparison
of even and odd cats with the classical mixture of ±α. The
[Fock-basis examples](@ref fock-basis-examples) cover phase-sensitive
moments, Kerr revivals and canonical thermal sampling. The
[positive-W example](@ref positive-w-sampling) uses exact number-state
draws and Weyl symbols for intensity moments.

The executable `examples/exact_sampling.jl` demonstrates a cat, photon-added
thermal and coherent states, NOON, a complex fixed-N condensate, a lossy cat,
phase-diffused squeezing and photon-added squeezing. It also runs the helper
families: displaced Fock states, rotations, finite number and classical mixtures,
products, photon subtraction, compass and sparse-code sums, squeezed cats,
passive mixing, and the general Gaussian constructors. For example:

```julia
using QuantumPhaseSpace, Random
rng = MersenneTwister(12)
cat = EvenCat(1.2+0.4im)
Z, stats = sampleQ(rng, Loss(cat,0.65),10_000;diagnostics=true)
α, αplus = positiveP(rng,NOON(2;phase=π/3),10_000)
u = [sqrt(0.3),sqrt(0.7)*cis(0.4)]
field, fieldplus = positiveP(rng,FixedNumberCondensate(20,u),10_000)
```
