@testset "RNG, layout and validation" begin
    for state in (Coherent(1im),Fock(10),Thermal(0.2im,1),Squeezed(1,0.3im),
                  Crescent(1im,0.2,0.3),Bogoliubov(sqrt(2),1,4),SqueezedTwoMode(0.4,0.7))
        for sampler in (husimiQ,positiveP)
            @test sampler(MersenneTwister(14),state,20) == sampler(MersenneTwister(14),state,20)
            @test all(isempty,sampler(MersenneTwister(14),state,0))
            @test_throws ArgumentError sampler(state,-1)
            @test_throws ArgumentError sampler(state,1.5)
        end
        Z = sampleQ(MersenneTwister(4),state,20)
        @test size(Z) == (nmodes(state),20)
        dest = similar(Z)
        @test sampleQ!(MersenneTwister(4),dest,state) === dest
        @test Z == dest
        @test_throws DimensionMismatch sampleQ!(zeros(ComplexF64,nmodes(state)+1,2),state)
    end
    @test_throws ArgumentError Coherent(Inf)
    @test_throws ArgumentError Thermal(NaN,1)
    @test_throws ArgumentError Thermal(0,Inf)
    @test_throws ArgumentError Squeezed(0,NaN)
    @test_throws ArgumentError Crescent(0,0,Inf)
    @test_throws ArgumentError SqueezedTwoMode(0,Inf)
    @test_throws SamplingError reject(MersenneTwister(2),x->0.0,[-1,1],1,1;max_attempts=5)
    x,d = reject(MersenneTwister(2),x->1.0,[-1,1],20,1;diagnostics=true)
    @test d.proposals == d.accepted == 20
    @test acceptance_rate(d) == 1
    @test isnan(acceptance_rate(SamplingDiagnostics()))
    @test_throws ArgumentError reject(x->1.0,[-Inf,1],1,1)
    @test overlapcapability(Fock(2)) isa PureOverlap
    @test overlapcapability(Thermal(0,1)) isa SamplingOnly
    @test isfinite(real(logoverlap(Fock(10000),sqrt(10000))))
    @test logoverlap(Coherent(20),20) == 0
    a,ap = positiveP(MersenneTwister(32),Thermal(0.3im,2),100_000)
    @test ap == conj.(a)
    checkmean(ap.*a,2.09)
end
