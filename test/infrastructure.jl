@testset "RNG, return layout and validation" begin
    # Exercise the old single-mode and four-vector interfaces as well as
    # the matrix interface and a newer rejection sampler.
    for state in (Coherent(1im), SqueezedTwoMode(0.4,0.7), NOON(2))
        z = sampleQ(MersenneTwister(4), state, 20)
        @test size(z) == (nmodes(state), 20)
        dest = similar(z)
        @test sampleQ!(MersenneTwister(4), dest, state) === dest
        @test z == dest

        q = husimiQ(MersenneTwister(4), state, 20)
        @test length(q) == (state isa SqueezedTwoMode ? 4 : 2)
        @test q[2] == conj.(q[1])
        @test eltype(q[1]) === ComplexF64
        if state isa SqueezedTwoMode
            @test q[1] == z[1,:]
            @test q[3] == z[2,:]
            @test q[4] == conj.(q[3])
        end
        @test positiveP(MersenneTwister(14), state, 20) ==
              positiveP(MersenneTwister(14), state, 20)
    end
    @test all(isempty, positiveP(SqueezedTwoMode(0.4,0.7), 0))
    @test size(sampleQ(NOON(2), 0)) == (2, 0)
    @test_throws DimensionMismatch sampleQ!(zeros(ComplexF64, 3, 2), NOON(2))
    @test_throws ArgumentError sampleQ(Coherent(0), -1)
    @test_throws ArgumentError husimiQ(Fock(1), 1.5)
    @test_throws ArgumentError positiveP(NOON(2), -1)
    @test_throws MethodError positiveP("not a state", 1)

    @test_throws ArgumentError Coherent(Inf)
    @test_throws ArgumentError Fock(-1)
    @test_throws ArgumentError Thermal(0,-1)
    @test_throws ArgumentError Thermal(NaN,1)
    @test_throws ArgumentError Squeezed(0,NaN)
    @test_throws ArgumentError Crescent(0,0,Inf)
    @test_throws ArgumentError SqueezedTwoMode(-0.1,0)
    @test_throws ArgumentError SqueezedTwoMode(0,Inf)
    @test_throws ArgumentError Bogoliubov(0.1,0.1,1)
    @test_throws ArgumentError Bogoliubov(sqrt(2),1,-1)

    @test overlapcapability(Fock(2)) isa PureOverlap
    @test overlapcapability(Thermal(0,1)) isa SamplingOnly
    @test logoverlap(Coherent(20),20) == 0
    a, ap = positiveP(MersenneTwister(32), Thermal(0.3im,2), 20)
    @test ap == conj.(a)
    @test (a,ap) == glauberP(MersenneTwister(32), Thermal(0.3im,2), 20)
end

@testset "Bounded rejection helper" begin
    @test_throws SamplingError reject(MersenneTwister(2), x->0.0, [-1,1], 1, 1; max_attempts=5)
    x, d = reject(MersenneTwister(2), x->1.0, [-1,1], 20, 1; diagnostics=true)
    @test x isa Vector{Float64} && all(y -> -1 <= y <= 1, x)
    @test d.proposals == d.accepted == length(x) == 20
    @test acceptance_rate(d) == 1
    @test isnan(acceptance_rate(SamplingDiagnostics()))
    @test isempty(reject(x->1.0, [-1,1], 0, 1))
    @test_throws ArgumentError reject(x->1.0, [1,-1], 1, 1)
    @test_throws ArgumentError reject(x->1.0, [-Inf,1], 1, 1)
    @test_throws ArgumentError reject(x->1.0, [-1,1], 1, 0)
    for bad in (NaN, -1.0, 2.0)
        @test_throws ArgumentError reject(x->bad, [-1,1], 1, 1)
    end
end
