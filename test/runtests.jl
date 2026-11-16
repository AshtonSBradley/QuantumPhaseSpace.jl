#!/usr/bin/env julia

using QuantumPhaseSpace, Test, Statistics, Aqua, Random

@testset "Aqua" begin
    Aqua.test_all(QuantumPhaseSpace)
end

# Compare sample means at six standard errors. Complex estimators must also
# pass the imaginary-component check, even when the target is real.
function checkmean(samples, target; nsigma=6, atol=1e-10)
    components = eltype(samples) <: Real ? (real,) : (real, imag)
    for component in components
        values = component.(samples)
        @test abs(mean(values) - component(target)) <=
              nsigma * std(values) / sqrt(length(values)) + atol
    end
end

include("hilbert.jl")
include("scientific.jl")
include("infrastructure.jl")
include("families.jl")
include("superpositions.jl")
include("envelopes.jl")
include("multimode.jl")
include("channels.jl")
include("gaussian.jl")
include("positivew.jl")
include("fermions.jl")
