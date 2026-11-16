# [Fermionic truncated Wigner](@id fermionic-twa)

[`SlaterDeterminant`](@ref) and [`FermionThermal`](@ref) initialize the
fermionic truncated Wigner approximation (fTWA). They match the exact quantum
means and symmetrized covariances of the chosen bilinears. The Gaussian ansatz
fixes higher moments of the sampled ensemble. In the canonical thermal case,
the package draws an occupation pattern first and uses the Slater ensemble for
that pattern. [Davidson, Sels and Polkovnikov](https://arxiv.org/abs/1604.08664).

Each array has shape `modes × modes × samples`. The Hermitian array `ρ`
represents `cᵢ†cⱼ - δᵢⱼ/2`. Read occupation in mode `i` as `real(ρ[i,i,k]) + 1/2`
on draw `k`. With `pairing=true`, the result also includes `τ` for `cᵢcⱼ`.
`τ` is antisymmetric and its ensemble mean vanishes for the number-conserving
states here. Its fluctuations include both occupied pairs and empty pairs.

```@example fermionic_twa
using QuantumPhaseSpace, Random, Statistics

rng = MersenneTwister(136)
slater = SlaterDeterminant([1, 0, 1, 0])
ρ = wigner(rng, slater, 10_000)

energies = [-1.0, -0.25, 0.25, 1.0]
thermal = FermionThermal(energies, 2.0; particles=2)
samples = wigner(rng, thermal, 10_000; pairing=true)

(; slater_occupation=mean(real.(ρ[1,1,:])) + 1/2,
   slater_hopping_variance=mean(abs2, ρ[1,2,:]),
   thermal_occupation=mean(real.(samples.ρ[1,1,:])) + 1/2)
```

For a Slater determinant, pass occupied orbitals as columns to sample in a
different one-body basis. Pass a Hermitian one-body Hamiltonian to
`FermionThermal` to sample in the Hamiltonian's input basis. For a mode
transformation `c → U*c`, these conventions give
`ρ → conj(U)*ρ*transpose(U)` and `τ → U*τ*transpose(U)`.

```@docs
SlaterDeterminant
FermionThermal
```
