@testset "Independent odd-cat density and envelope" begin
    # Build normalised odd Fock coefficients directly by parity, avoiding any
    # subtraction of almost equal coherent kets or sampler overlap formulas.
    for α in (1e-6+2e-6im,0.3+0.2im,0.8im)
        A = abs2(α)
        K = 60
        v = coherentket(α,K)
        v[1:2:end] .= 0
        normalize!(v)
        rate = -expm1(-2A)*exp(-A)/(8A)
        for r in (0.01,0.3,1.0,2.0,4.0), θ in (0.1,0.8,2.0)
            z = r*cis(θ)
            t = conj(α)*z
            # Direct expression is safe on this modest verification grid.
            p = exp(-r^2/2-2A)*abs2(sinh(t)/t)
            proposal = r^2*exp(-r^2/2)/(4π)
            target = abs2(dot(coherentket(z,K),v))/π
            @test 0<=p<=1
            @test proposal*p ≈ rate*target rtol=1e-11
        end
    end
end
