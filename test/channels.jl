# Kraus pure-loss channel in an independently truncated number basis.
function lossdensity(v,η)
    K = length(v)
    ρ = zeros(ComplexF64,K,K)
    for lost in 0:K-1
        w = zeros(ComplexF64,K)
        for n in lost:K-1
            w[n-lost+1] = sqrt(binomial(big(n),lost))*η^((n-lost)/2)*(1-η)^(lost/2)*v[n+1]
        end
        ρ .+= w*w'
    end
    return ρ
end
@testset "Loss and phase diffusion" begin
    rng = MersenneTwister(505)
    β = 0.6+0.4im
    K = 42
    a = annihilator(K)
    N = 45_000
    coherent = Coherent(β)
    @test sampleQ(MersenneTwister(8),Loss(coherent,1),20) ==
          sampleQ(MersenneTwister(8),coherent,20)
    z = vec(sampleQ(rng,Loss(coherent,0),N))
    checkmean(z,0)
    checkmean(abs2.(z),1)

    n, η = 4, 0.4
    p = [pdf(Binomial(n,η),k) for k in 0:n]
    for state in (Loss(Fock(n),η),NumberMixture(p))
        z = vec(sampleQ(rng,state,N))
        checkmean(abs2.(z),1+n*η)
        checkmean(abs2.(z).^2,2+4n*η+n*(n-1)*η^2)
    end

    phase = 0.8
    cat = Cat(β,phase)
    v = normalize!(coherentket(β,K)+cis(phase)*coherentket(-β,K))
    ρ = lossdensity(v,η)
    @test tr(ρ) ≈ 1 atol=1e-12
    z = vec(sampleQ(rng,Loss(cat,η),N))
    # First and mixed moments retain the interference-sensitive odd part.
    checkmean(z,tr(ρ*a))
    checkmean(conj.(z).*z.^2,tr(ρ*a^2*a'))
    checkmean(abs2.(z),tr(ρ*a*a'))

    σ = 0.6
    z = vec(sampleQ(rng,PhaseDiffused(coherent,σ),N))
    checkmean(z,β*exp(-σ^2/2))
    checkmean(z.^2,β^2*exp(-2σ^2))
    checkmean(abs2.(z),1+abs2(β))
    for common in (true,false)
        s = PhaseDiffused(ProductState(Coherent(β),Coherent(im*β)),σ;common)
        Z = sampleQ(rng,s,N)
        factor = common ? 1.0 : exp(-σ^2)
        checkmean(conj.(Z[1,:]).*Z[2,:],im*abs2(β)*factor)
    end
    s = ProductState(Fock(3),Fock(2))
    Z = sampleQ(rng,Loss(s,[0.2,0.8]),N)
    checkmean(abs2.(Z[1,:]),1.6)
    checkmean(abs2.(Z[2,:]),2.6)
    @test sampleQ(MersenneTwister(2),PhaseDiffused(cat,0),20) == sampleQ(MersenneTwister(2),cat,20)
    @test_throws ArgumentError Loss(cat,1.1)
    @test_throws ArgumentError PhaseDiffused(cat,-1)
    @test_throws ArgumentError PhaseDiffused(s,[0.1,0.2])
    @test_throws DimensionMismatch Loss(s,[0.2])
    @test overlapcapability(Loss(cat,0.4)) isa SamplingOnly
end
