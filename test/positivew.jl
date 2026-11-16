@testset "Exact doubled positive-W" begin
    rng = MersenneTwister(707)
    for (n, r) in ((1,0.25), (15,6.0))
        mixture = sum(pdf(Binomial(n,2/3),k)*pdf(Gamma(k+1,1.5),r^2)*2r for k in 0:n)
        @test mixture ≈ plaguerre(r,n) rtol=1e-12
    end
    @test isfinite(plaguerre_asymptotic(0.1,1))

    # Vacuum, a nontrivial photon number, and a large occupation beyond
    # the old approximation switch. Both paths are compared at n=3.
    N = 45_000
    for (n, method) in ((0,:direct), (3,:direct), (3,:convolution), (10_000,:direct))
        a, ap = positiveW(rng,Fock(n),N;method)
        R = a.*ap
        @test length(a) == length(ap) == N
        checkmean(a,0)
        checkmean(R,n+0.5)
        checkmean(R.^2,n^2+n+0.5)
        checkmean((real.(R) .- (n+0.5)).^2,2n+13/4)
    end
    a, ap = positiveW(rng,Coherent(0.4im),N)
    @test ap == conj.(a)
    checkmean(abs2.(a),0.16+0.5)
    @test all(isempty,positiveW(rng,SqueezedTwoMode(0.3,0.4),0))
    s = NOON(2)
    @test positiveW(MersenneTwister(8),s,20) == positiveW(MersenneTwister(8),s,20)
    @test_throws ArgumentError positiveW(Fock(1),3;method=:approximate)
end
