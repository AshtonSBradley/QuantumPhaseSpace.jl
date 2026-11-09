# Distributions
The phase-space distributions that may be sampled for these quantum states are

## Husimi-Q
```@docs
husimiQ
```
## Glauber-P
```@docs
glauberP
```
## Positive-P
```@docs
positiveP
```
## Wigner
```@docs
wigner
```
## Positive-W
```@docs
positiveW
```

## Reproducible and multimode sampling

Every sampler accepts `sampler(rng, state, N)` with an explicit `AbstractRNG`;
the original `sampler(state, N)` uses `Random.default_rng()`. The historical
single-mode vectors and `SqueezedTwoMode` four-vector tuple are preserved.
New multimode states return a pair of modes × samples matrices.

`sampleQ(rng,state,N)` always returns one modes × samples matrix, even for a
single mode. `sampleQ!` fills a caller-owned matrix. Both accept zero columns.
The defaults use Float64 and circular `randn(rng,ComplexF64)` noise with
`E[abs2(ξ)]=1`. A trajectory's norm is never fixed to a particle count.

```@docs
nmodes
sampleQ
sampleQ!
SamplingDiagnostics
SamplingError
acceptance_rate
overlapcapability
PureOverlap
SamplingOnly
logoverlap
```
