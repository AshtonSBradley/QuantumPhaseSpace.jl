# Sampling benchmark report

Measured 10 November 2026 on Apple M4 Pro (aarch64 macOS), Julia 1.13.3, one Julia thread, Float64. Source: sampling implementation after stages 0–6 and Wplus. No GPU was used. Measured package versions: Distributions 0.25.134, GSL 1.0.1, Parameters 0.12.3, SpecialFunctions 2.9.1.

## Method

`benchmark/run.jl` warms every constructor and sampling callback, then records the median elapsed time and allocated bytes of five runs for batches of 100, 10,000 and 100,000. RNG construction and output allocation are included in sampling time. Preparation is measured separately after warm-up. Estimator calculations are outside the timer. Values reported as zero preparation time are below timer resolution, not a claim of zero cost. These are local measurements, not performance guarantees.

The physical estimator is occupation in mode 1. Q uses `|z₁|²−1`, positive-P uses `αplus₁ α₁`, and positive-W uses `αplus₁ α₁−1/2`. Real and imaginary means and variances remain separate. The TSV includes preparation time/bytes, batch time/bytes, acceptance, target, means, variances, standard error and `variance_real * seconds / N`. This last quantity predicts the real-estimator mean-square error after one second of sampling, assuming independent initial samples; lower is better. It describes initial sampling only, not later stochastic dynamics.

## Large batches

N=100,000. P denotes positive-P; W+ denotes doubled positive-W. All exact rejection algorithms have attempt limits. Non-rejection rows have acceptance 1.

| State | Representation / method | Batch ms | Bytes/sample | Acceptance | Preparation µs | Var(real n̂) | Var(imag n̂) | Variance × time/sample |
|:--|:--|--:|--:|--:|--:|--:|--:|--:|
| even_cat | Q / default | 6.786 | 32.3 | 0.5204 | 0.000 | 4.314 | 0 | 2.93e-07 |
| even_cat | P / generic | 7.494 | 64.4 | 0.5204 | 0.000 | 5.316 | 4.934 | 3.98e-07 |
| small_odd_cat | Q / default | 38.688 | 32.3 | 0.2489 | 0.000 | 1.998 | 0 | 7.73e-07 |
| small_odd_cat | P / generic | 40.537 | 64.4 | 0.2489 | 0.000 | 3.009 | 3.974 | 1.22e-06 |
| added_thermal | Q / default | 2.384 | 32.3 | 1.0000 | 0.041 | 8.695 | 0 | 2.07e-07 |
| added_thermal | P / generic | 2.768 | 64.4 | 1.0000 | 0.000 | 9.716 | 10.17 | 2.69e-07 |
| added_coherent | Q / default | 18.643 | 32.3 | 0.5051 | 18.250 | 4.12 | 0 | 7.68e-07 |
| added_coherent | P / generic | 19.526 | 64.4 | 0.5051 | 13.542 | 5.127 | 6.329 | 1e-06 |
| NOON_2 | Q / default | 32.359 | 64.4 | 0.4992 | 18.625 | 3.027 | 0 | 9.8e-07 |
| NOON_2 | P / generic | 33.678 | 128.7 | 0.4992 | 20.959 | 4.032 | 4.02 | 1.36e-06 |
| condensate_30_M16 | P / specialised | 8.248 | 640.8 | 1.0000 | 0.042 | 0.1248 | 0.2421 | 1.03e-08 |
| condensate_30_M16 | P / generic | 43.028 | 1024.5 | 1.0000 | 0.042 | 5.63 | 5.74 | 2.42e-06 |
| thermal | P / specialised | 0.893 | 64.4 | 1.0000 | 0.000 | 9.89 | 0 | 8.84e-08 |
| thermal | P / generic | 1.398 | 128.7 | 1.0000 | 0.000 | 18.18 | 8.252 | 2.54e-07 |
| Fock_0 | Wplus / direct | 3.288 | 32.3 | 1.0000 | 0.000 | 3.249 | 2.988 | 1.07e-07 |
| Fock_0 | Wplus / convolution | 3.180 | 192.9 | 1.0000 | 0.000 | 3.235 | 2.975 | 1.03e-07 |
| Fock_1 | Wplus / direct | 5.131 | 32.3 | 1.0000 | 0.041 | 5.221 | 4.98 | 2.68e-07 |
| Fock_1 | Wplus / convolution | 3.925 | 192.9 | 1.0000 | 0.000 | 5.246 | 4.992 | 2.06e-07 |
| Fock_320 | Wplus / direct | 8.774 | 32.3 | 1.0000 | 0.000 | 641 | 637.8 | 5.62e-05 |
| Fock_320 | Wplus / convolution | 3.781 | 192.9 | 1.0000 | 0.000 | 640.8 | 639 | 2.42e-05 |
| Fock_10000 | Wplus / direct | 7.410 | 32.3 | 1.0000 | 0.000 | 1.997e+04 | 1.994e+04 | 0.00148 |
| Fock_10000 | Wplus / convolution | 3.787 | 192.9 | 1.0000 | 0.041 | 2.005e+04 | 1.982e+04 | 0.000759 |

## Same-representation comparisons

- The collective-mode specialised positive-P sampler reduces mode-1 occupation variance by 45.1× relative to generic Q conversion, and the measured variance–time cost by 235.4× for this N=30, M=16 example. Its zero orthogonal positive-P vacuum is physically valid; no trajectory norm is fixed.
- Thermal specialised positive-P reduces occupation variance by 1.84× and the measured variance–time cost by 2.88× in this displaced nbar=3 case.
- Direct Fock positive-W uses about six times fewer allocated bytes than convolution in these batches. Convolution is faster in the measured cases, especially at high n; their estimator variances agree within sampling uncertainty. The exact direct law is the default to remove the old spatial-window and asymptotic approximations, not because it wins every runtime comparison.
- In this archived benchmark, the broader Gamma proposal for small odd cats gives acceptance approximately 1/4 even at |α|=0.0001. Even-cat, photon-added coherent and NOON acceptance agree with their analytic expectations. Photon-added coherent preparation includes a small exact finite Gram problem; no coherent-state truncation is hidden.

## Batch-size scaling

| State | Representation / method | 100 samples (ms) | 10,000 (ms) | 100,000 (ms) |
|:--|:--|--:|--:|--:|
| even_cat | Q / default | 0.0135 | 0.660 | 6.786 |
| even_cat | P / generic | 0.0146 | 0.747 | 7.494 |
| small_odd_cat | Q / default | 0.0499 | 4.022 | 38.688 |
| small_odd_cat | P / generic | 0.0481 | 3.841 | 40.537 |
| added_thermal | Q / default | 0.0094 | 0.210 | 2.384 |
| added_thermal | P / generic | 0.0100 | 0.287 | 2.768 |
| added_coherent | Q / default | 0.0274 | 1.883 | 18.643 |
| added_coherent | P / generic | 0.0275 | 1.977 | 19.526 |
| NOON_2 | Q / default | 0.0429 | 3.177 | 32.359 |
| NOON_2 | P / generic | 0.0399 | 3.348 | 33.678 |
| condensate_30_M16 | P / specialised | 0.0114 | 0.475 | 8.248 |
| condensate_30_M16 | P / generic | 0.0268 | 2.536 | 43.028 |
| thermal | P / specialised | 0.0078 | 0.062 | 0.893 |
| thermal | P / generic | 0.0085 | 0.177 | 1.398 |
| Fock_0 | Wplus / direct | 0.0103 | 0.334 | 3.288 |
| Fock_0 | Wplus / convolution | 0.0109 | 0.361 | 3.180 |
| Fock_1 | Wplus / direct | 0.0125 | 0.525 | 5.131 |
| Fock_1 | Wplus / convolution | 0.0110 | 0.384 | 3.925 |
| Fock_320 | Wplus / direct | 0.0216 | 0.895 | 8.774 |
| Fock_320 | Wplus / convolution | 0.0104 | 0.350 | 3.781 |
| Fock_10000 | Wplus / direct | 0.0138 | 0.721 | 7.410 |
| Fock_10000 | Wplus / convolution | 0.0102 | 0.385 | 3.787 |

## Accuracy and limits

All 66 measured means were within 2.44 estimated standard errors of the independent analytic occupation targets. This is a seeded validation observation, not a universal accuracy bound. The unit tests separately check complex anomalous, factorial and cross-mode moments and selected finite-Fock density matrices.

Comparisons across Q/P/W+ rows show different estimators of the same physical occupation; efficiency claims above compare the same representation. A fast exact Q sampler does not guarantee small positive-P initial variance, and initial variance does not predict time-dependent positive-P stability. High-order observables can have much larger variance than the occupation benchmark.

All laws in this benchmark are mathematically exact for the supplied states. Float64 precision, near cancellation and attempt-limit exhaustion remain numerical limits. The separate ordinary Fock W sampler is approximate for n>0 and is intentionally absent from exact positive-W speed comparisons. Weighted ordinary W, GPU sampling, automatic differentiation and direct general Gaussian positive-P factorisation remain outside this implementation.

## Reproduce

From a juliaup-managed persistent Julia session in the project, after `Pkg.instantiate()`:

```julia
include("benchmark/run.jl")
run_sampling_benchmarks()
```

The complete inspectable measurements are in [`results.tsv`](results.tsv). The runnable state examples are in [`../examples/exact_sampling.jl`](../examples/exact_sampling.jl).

## Direct Gaussian positive-P for two-mode squeezing

Run `julia --project=. benchmark/gaussian_positive_p.jl` to reproduce the
comparison of `positiveP(rng,SqueezedTwoMode(r,0),N;method)` for
`:canonical` and `:drummond`. Measured on Apple M4 Pro, Julia 1.13.2,
one Julia thread, Float64. Times are medians of nine warmed calls at
N=200,000, including output allocation. RNG creation and a garbage
collection occur before each timed call. Estimator standard errors are measured on a separate batch of the
same size. D is the mean of the two modewise values of
`E|α-conj(α⁺)|²`. The benchmark prints both modewise values and both
occupation estimates.

Each SE cell gives real / imaginary standard errors. The exact number-difference
variance is zero in the squeezed vacuum. Its trajectory estimator with the
normal-ordering correction is `(n₁-n₂)^2+n₁+n₂`. Its imaginary component,
and that of `ab`, vanish to roundoff for the direct sampler at phase zero.

| r | Method | D | SE of ⟨n₁⟩ | SE of ⟨ab⟩ | SE of Var(n₁−n₂) | ms | Bytes/sample |
| ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: |
| 0 | Canonical | 4.00 | .00316 / .00316 | .00317 / .00316 | .0224 / .0223 | 5.36 | 225 |
| 0 | Drummond | 0 | 0 / 0 | 0 / 0 | 0 / 0 | 8.17 | 128 |
| 0.5 | Canonical | 4.00 | .00362 / .00357 | .00372 / .00347 | .0274 / .0273 | 5.16 | 225 |
| 0.5 | Drummond | .634 | .00102 / .000822 | .00131 / <1e−14 | .00292 / <1e−14 | 8.03 | 128 |
| 1 | Canonical | 4.00 | .00579 / .00488 | .00609 / .00452 | .0487 / .0473 | 5.32 | 225 |
| 1 | Drummond | .867 | .00360 / .00185 | .00406 / <1e−14 | .0128 / <1e−14 | 8.03 | 128 |
| 2 | Canonical | 4.00 | .0319 / .0119 | .0324 / .0105 | .282 / .258 | 5.21 | 225 |
| 2 | Drummond | .984 | .0299 / .00571 | .0305 / <1e−14 | .117 / <1e−14 | 8.06 | 128 |
| 3 | Canonical | 4.00 | .228 / .0319 | .228 / .0277 | 2.02 / 1.81 | 5.24 | 225 |
| 3 | Drummond | 1.00 | .225 / .0158 | .225 / <1e−14 | .907 / <1e−14 | 8.14 | 128 |

At r=0.5 the direct sampler cut the real estimator variances for occupation,
anomalous correlation and number-difference variance by factors of about
12.5, 8.0 and 88, respectively. The factors at r=2 were 1.14, 1.13 and 5.8.
The direct sampler took about 1.5 times as long to draw this batch and
allocated 43% fewer bytes. These timings cover initial-state
sampling only.

High-order convergence was measured as relative complex RMSE over 32
independent batches for each method. The exact values for the normally
ordered occupation estimator `(α⁺α)^4` and the pair estimators
`(α₁α₂)^k` are `4! sinh(r)^8` and `k! (sinh(r)cosh(r))^k`,
respectively. Each cell gives canonical / Drummond RMSE.

| r | N | ⟨(a†)^4 a^4⟩ | ⟨(ab)^4⟩ | ⟨(ab)^6⟩ |
| ---: | ---: | ---: | ---: | ---: |
| 0.5 | 1,000 | 171 / 1.20 | 10.4 / .190 | 162 / .419 |
| 0.5 | 10,000 | 59.9 / .532 | 4.05 / .0800 | 161 / .249 |
| 0.5 | 100,000 | 15.6 / .176 | 1.06 / .0244 | 18.7 / .0936 |
| 2 | 1,000 | .340 / .206 | .268 / .185 | .815 / .399 |
| 2 | 10,000 | .103 / .0822 | .114 / .0716 | .766 / .206 |
| 2 | 100,000 | .0346 / .0261 | .0296 / .0264 | .120 / .0830 |

High-order errors vary across finite batches because rare, large polynomial
values contribute strongly. Strong squeezing also requires subtracting
large contributions when estimating a small number-difference variance.
The sampler operates on Float64 covariances and the comparison covers
initial states only.
