@testset "Odd-cat envelope at the node and in the tail" begin
    # Compare the log implementation with odd Fock coefficients, which
    # retain the one-photon limit without subtracting nearly equal kets.
    for (α, z) in ((1e-6+2e-6im, 0.3+0.2im), (0.35im, 0.0im),
                   (0.35im, 0.0+3.0im), (0.35im, 4.0+0.5im))
        A, B = abs2(α), abs2(z)
        λ = 1-A/3
        v = coherentket(α, 48)
        v[1:2:end] .= 0
        normalize!(v)
        p = exp(-A*B/3+2real(QuantumPhaseSpace._logsinhc(conj(α)*z)))
        proposal = λ^2*B*exp(-λ*B)/π
        target = abs2(dot(coherentket(z, length(v)), v))/π
        @test 0 <= p <= 1
        @test proposal*p ≈ λ^2*sinh(A)/A*target rtol=1e-11
    end
end
