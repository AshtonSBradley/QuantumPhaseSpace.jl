using QuantumPhaseSpace: FockSuperposition, Mixture, NumberMixture, sampleQ, positiveP
using LinearAlgebra: Diagonal, Hermitian, I, diag, diagm, eigen, isdiag, ishermitian, tr
using Random: MersenneTwister
using Statistics: mean, std
using Printf: @printf, @sprintf

"""
    canonical_fock_state(H, βT)

Prepare the Gibbs state of a finite Hermitian Hamiltonian in the single-mode
Fock basis |0⟩, …, |K-1⟩. βT is the inverse temperature in the units of H.
Return both the density matrix and the state to pass to sampleQ or positiveP.
"""
function canonical_fock_state(H::AbstractMatrix, βT::Real)
    size(H, 1) > 0 && ishermitian(H) && all(isfinite, H) ||
        throw(ArgumentError("H must be a finite, nonempty Hermitian matrix"))
    βT = Float64(βT)
    isfinite(βT) && βT > 0 || throw(ArgumentError("βT must be finite and positive"))

    if isdiag(H)
        energies = Float64.(real.(diag(H)))
        weights = boltzmann_probabilities(energies, βT)
        return (; state=NumberMixture(weights), ρ=Diagonal(weights))
    end

    energies, eigenvectors = eigen(Hermitian(Matrix{ComplexF64}(H)))
    weights = boltzmann_probabilities(energies, βT)
    occupied = findall(>(0), weights)
    occupations = 0:(size(H, 1)-1)
    states = [FockSuperposition(occupations, view(eigenvectors, :, j)) for j in occupied]
    ρ = eigenvectors * Diagonal(weights) * eigenvectors'
    return (; state=Mixture(states, weights[occupied]), ρ)
end

function boltzmann_probabilities(energies, βT)
    all(isfinite, energies) || throw(ArgumentError("energies must be finite in Float64"))
    weights = exp.(-βT .* (energies .- minimum(energies)))
    return weights ./ sum(weights)
end

function compare_thermal_moment(values, target)
    estimate = mean(values)
    errors = (std(real.(values)), std(imag.(values))) ./ sqrt(length(values))
    δ = estimate-target
    atol = 64eps(Float64)*max(1, abs(target))
    largest_z = max(abs(real(δ))/max(errors[1], atol),
                    abs(imag(δ))/max(errors[2], atol))
    passed = abs(real(δ)) <= 6errors[1]+atol && abs(imag(δ)) <= 6errors[2]+atol
    return (; estimate, real_se=errors[1], imag_se=errors[2], largest_z, passed)
end

thermal_complex(z) = @sprintf("%.6f%+.6fi", real(z), imag(z))

"""
    check_canonical_thermal(; N=100_000, K=36, rng=MersenneTwister(20261112))

Check Q and positive-P samples of the canonical Kerr oscillator
H/ℏ = ω n + χ n(n-1)/2 and of harmonic thermal states. The benchmark
sets ℏ = 1 and βT = ω = χ = 1, so the Kerr vacuum is the unique ground
state. Compare against a converged Kerr number sum and closed-form
harmonic moments. The occupation cutoff is K-1.
"""
function check_canonical_thermal(; N=100_000, K=36, rng=MersenneTwister(20261112))
    N isa Integer && N >= 2 || throw(ArgumentError("N must be an integer ≥ 2"))
    K isa Integer && K >= 2 || throw(ArgumentError("K must be an integer ≥ 2"))
    βT = 1.0
    ω, χ = 1.0, 1.0
    n = Float64.(0:K-1)
    kerr_energies = @. ω*n + χ*n*(n-1)/2
    a = diagm(1 => sqrt.(n[2:end]))
    n_op = a'*a
    operators = (a, a^2, n_op, a'*a^2, a'^2*a^2)
    labels = ("a", "a²", "n", "a†a²", "n(n−1)")

    # The Kerr spectrum grows quadratically. A longer independent sum
    # gives the infinite-basis reference to working precision here.
    n_ref = Float64.(0:4K)
    w_ref = @. exp(-βT*ω*n_ref) * exp(-βT*χ*n_ref*(n_ref-1)/2)
    Z_ref = sum(w_ref)
    kerr_n = sum(n_ref .* w_ref)/Z_ref
    kerr_f2 = sum(n_ref .* (n_ref .- 1) .* w_ref)/Z_ref
    # Above the cutoff, successive Boltzmann weight ratios decrease.
    first_omitted_energy = ω*K + χ*K*(K-1)/2
    first_omitted_gap = ω + χ*K
    kerr_tail_bound = exp(-βT*first_omitted_energy)/(-expm1(-βT*first_omitted_gap)*Z_ref)

    nbar = 1/expm1(βT)
    displacement = 0.65+0.35im
    b = a-displacement*I
    d2 = abs2(displacement)
    cases = (
        ("Kerr oscillator", Diagonal(kerr_energies), (0, 0, kerr_n, 0, kerr_f2)),
        ("Harmonic", Diagonal(n), (0, 0, nbar, 0, 2nbar^2)),
        ("Displaced harmonic", b'*b,
            (displacement, displacement^2, d2+nbar,
             displacement*(d2+2nbar), d2^2+4d2*nbar+2nbar^2)),
    )

    results = NamedTuple[]
    @printf("βT = %.1f, Fock levels 0:%d, %d draws per representation\n", βT, K-1, N)
    @printf("Kerr ω = %.1f, χ = %.1f, omitted probability ≤ %.3g\n", ω, χ, kerr_tail_bound)
    for (name, H, targets) in cases
        (; state, ρ) = canonical_fock_state(H, βT)
        z = vec(sampleQ(rng, state, N))
        α, αplus = positiveP(rng, state, N)
        radius2 = abs2.(z)
        nplus = αplus .* α
        Q_values = (z, z.^2, radius2 .- 1, z .* (radius2 .- 2),
                    radius2.^2 .- 4radius2 .+ 2)
        P_values = (α, α.^2, nplus, α .* nplus, nplus.^2)
        moments = NamedTuple[]
        println("\n", name)
        println("Moment        analytic                     Q estimate                    +P estimate")
        for j in eachindex(operators)
            target = targets[j]
            projected = tr(ρ*operators[j])
            cutoff_error = abs(projected-target)
            cutoff_error <= 2e-8*max(1, abs(target)) ||
                error("Fock cutoff moment check failed for $name and $(labels[j])")
            Q = compare_thermal_moment(Q_values[j], projected)
            P = compare_thermal_moment(P_values[j], projected)
            Q.passed && P.passed ||
                error("Sample moment exceeded six standard errors for $name and $(labels[j])")
            @printf("%-9s %24s   %24s   %24s\n",
                labels[j], thermal_complex(target), thermal_complex(Q.estimate),
                thermal_complex(P.estimate))
            push!(moments, (; observable=labels[j], target, projected, cutoff_error, Q, P))
        end
        largest_z = maximum(m -> max(m.Q.largest_z, m.P.largest_z), moments)
        cutoff_error = maximum(m -> m.cutoff_error, moments)
        @printf("Largest component error %.3f standard errors, largest cutoff error %.3g\n",
            largest_z, cutoff_error)
        push!(results, (; name, state, ρ, z, α, αplus, moments, largest_z, cutoff_error))
    end
    return results
end

if abspath(PROGRAM_FILE) == @__FILE__
    check_canonical_thermal()
end
