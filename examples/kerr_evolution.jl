using QuantumPhaseSpace: FockSuperposition, sampleQ, positiveP, acceptance_rate
using Distributions: Poisson, pdf, quantile, ccdf
using LinearAlgebra: diagm, dot, norm
using Random: MersenneTwister
using Statistics: mean, std
using Printf: @printf, @sprintf

"""
    kerr_coherent_state(β, τ; tail_tol=1e-14)

Prepare the projected coherent state evolved under H/ℏ = χ n(n-1)/2,
where τ = χt. The discarded Poisson probability is independent of τ.
The returned sampler uses the normalized projection onto 0:nmax.
"""
function kerr_coherent_state(β::Number, τ::Real; tail_tol=1e-14)
    β, τ = ComplexF64(β), Float64(τ)
    ν = abs2(β)
    isfinite(ν) && isfinite(τ) || throw(ArgumentError("finite occupation and time required"))
    0 < 1-tail_tol < 1 || throw(ArgumentError("1 - tail_tol must lie strictly between zero and one"))

    τ = rem2pi(τ, RoundNearest)
    photons = Poisson(ν)
    nmax = Int(quantile(photons, 1-tail_tol))
    # Quantile selection rounds the requested cumulative probability.
    while ccdf(photons, nmax) > tail_tol
        nmax += 1
    end
    occupations = 0:nmax
    θ = angle(β)
    ket = [sqrt(pdf(photons, n))*cis(n*θ - τ*n*(n-1)/2) for n in occupations]
    ket ./= norm(ket)
    state = FockSuperposition(occupations, ket)
    return (; state, ket, nmax, discarded_probability=ccdf(photons, nmax))
end

# With n=m+q and n'=m+p in the matrix element, the remaining m-sum is
# the generating function of a Poisson random variable with mean |β|².
function kerr_moment(β, τ, p, q)
    τ = rem2pi(Float64(τ), RoundNearest)
    phase = cis(τ*(p*(p-1) - q*(q-1))/2)
    return conj(β)^p * β^q * phase * exp(abs2(β)*expm1(im*τ*(p-q)))
end

# Commute the normally ordered monomial to antinormal order for Q samples.
function kerr_q_symbol(z, p, q)
    r = min(p, q)
    radius2 = abs2(z)
    polynomial = coefficient = 1.0
    for k in 1:r
        coefficient *= -(p-k+1)*(q-k+1)/k
        polynomial = muladd(polynomial, radius2, coefficient)
    end
    return conj(z)^(p-r) * z^(q-r) * polynomial
end

function compare_kerr_moment(samples, target)
    estimate = mean(samples)
    real_se = std(real.(samples))/sqrt(length(samples))
    imag_se = std(imag.(samples))/sqrt(length(samples))
    error = estimate-target
    atol = 64eps(Float64)*max(1, abs(target))
    largest_z = max(abs(real(error))/max(real_se, atol),
                    abs(imag(error))/max(imag_se, atol))
    passed = abs(real(error)) <= 6real_se+atol && abs(imag(error)) <= 6imag_se+atol
    return (; estimate, real_se, imag_se, largest_z, passed)
end

kerr_complex(z) = @sprintf("%.6f%+.6fi", real(z), imag(z))

"""
    check_kerr_evolution(; β=2cis(0.35), times=(0.37, 2π/3, π, 2π), N=100_000,
                           rng=MersenneTwister(20261112), tail_tol=1e-14)

Compare fresh Husimi-Q and positive-P samples with normally ordered moments
of the full Kerr-evolved coherent state. An independent ladder-matrix
calculation measures the moment error from the Fock cutoff. Return the
samples, statistical errors, and acceptance diagnostics for each time.
"""
function check_kerr_evolution(; β=2cis(0.35), times=(0.37, 2π/3, π, 2π),
                               N=100_000, rng=MersenneTwister(20261112),
                               tail_tol=1e-14)
    N isa Integer && N >= 2 || throw(ArgumentError("N must be an integer ≥ 2"))
    β = ComplexF64(β)
    orders = ((0,1), (0,2), (0,3), (1,1), (1,2), (1,3), (2,2), (3,3))
    results = NamedTuple[]
    println("Kerr H/ℏ = χ n(n−1)/2, β = ", kerr_complex(β), ", ", N, " draws per representation")
    for τ in times
        (; state, ket, nmax, discarded_probability) = kerr_coherent_state(β, τ; tail_tol)
        a = diagm(1 => sqrt.(1.0:nmax))
        moments = NamedTuple[]
        Z, diagnostics = sampleQ(rng, state, N; diagnostics=true)
        z = vec(Z)
        α, αplus = positiveP(rng, state, N)
        rate = acceptance_rate(diagnostics)
        @printf("\nτ = %.8f, nmax = %d, discarded mass = %.3g, Q acceptance = %.3f%% (expected %.3f%%)\n",
            τ, nmax, discarded_probability, 100rate, 100state.expected_acceptance)
        println("(p,q)       analytic                     Q estimate                    +P estimate")
        for (p, q) in orders
            exact = kerr_moment(β, τ, p, q)
            projected = dot(ket, a'^p * a^q * ket)
            cutoff_error = abs(projected-exact)
            cutoff_error <= 1e-9*max(1, abs(β)^(p+q)) ||
                error("Fock cutoff moment check failed at τ=$τ for order ($p,$q)")
            Q = compare_kerr_moment(kerr_q_symbol.(z, p, q), projected)
            P = compare_kerr_moment(αplus.^p .* α.^q, projected)
            Q.passed && P.passed ||
                error("Sample moment exceeded six standard errors at τ=$τ for order ($p,$q)")
            @printf("(%d,%d)   %24s   %24s   %24s\n",
                p, q, kerr_complex(exact), kerr_complex(Q.estimate), kerr_complex(P.estimate))
            push!(moments, (; p, q, exact, projected, cutoff_error, Q, P))
        end
        largest_z = maximum(m -> max(m.Q.largest_z, m.P.largest_z), moments)
        cutoff_error = maximum(m -> m.cutoff_error, moments)
        @printf("Largest component error: %.3f standard errors; largest absolute cutoff error: %.3g\n",
            largest_z, cutoff_error)
        push!(results, (; τ, state, z, α, αplus, diagnostics, moments, nmax,
                        discarded_probability, largest_z, cutoff_error))
    end
    return results
end

if abspath(PROGRAM_FILE) == @__FILE__
    check_kerr_evolution()
end
