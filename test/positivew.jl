@testset "Exact doubled positive-W" begin
    rng = MersenneTwister(707)
    # Independent Laguerre/Gamma density identity in the safely evaluable range.
    for n in (0,1,4,15),r in (0.1,1.0,3.0,6.0)
        mixture = sum(pdf(Binomial(n,2/3),k)*pdf(Gamma(k+1,1.5),r^2)*2r for k in 0:n)
        @test mixture ≈ plaguerre(r,n) rtol=1e-12
    end
    for n in (0,1,40,319,320,321,1000,10000),method in (:direct,:convolution)
        a,ap = positiveW(rng,Fock(n),100_000;method)
        checkmean(a,0)
        checkmean(a.*ap,n+0.5)
        checkmean(a.^2 .* ap.^2,n^2+n+0.5)
        @test length(a)==length(ap)==100_000
    end
    for s in (EvenCat(0.7im),NOON(2),Thermal(0,0.5),SqueezedTwoMode(0.3,0.4))
        @test all(isempty,positiveW(rng,s,0))
        @test positiveW(MersenneTwister(8),s,20) == positiveW(MersenneTwister(8),s,20)
    end
    # Coherent +W convolution is also an ordinary Gaussian; this exceptional
    # conjugate case does not make doubled sampling an ordinary W interface.
    a,ap = positiveW(rng,Coherent(0.4im),100_000)
    @test ap==conj.(a)
    checkmean(a.*ap,0.16+0.5)
    @test_throws ArgumentError positiveW(Fock(1),3;method=:approximate)
end
