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
    for η in (0.0,0.4,1.0)
        s = Loss(Coherent(β),η)
        z = vec(sampleQ(rng,s,90_000))
        checkmean(z,sqrt(η)*β)
        checkmean(abs2.(z).-1,η*abs2(β))
        checkmean(z.^2,η*β^2)
        if η==1
            @test sampleQ(MersenneTwister(8),s,20) == sampleQ(MersenneTwister(8),Coherent(β),20)
        end
        n = 4
        p = [pdf(Binomial(n,η),k) for k in 0:n]
        for state in (Loss(Fock(n),η),NumberMixture(p))
            z = vec(sampleQ(rng,state,90_000))
            checkmean(abs2.(z),1+n*η)
            checkmean(abs2.(z).^2,2+4n*η+n*(n-1)*η^2)
        end
    end
    cat = Cat(β,0.8)
    v = coherentket(β,K)+cis(0.8)*coherentket(-β,K); normalize!(v)
    for η in (0.3,0.8)
        ρ = lossdensity(v,η)
        @test tr(ρ) ≈ 1 atol=1e-12
        z = vec(sampleQ(rng,Loss(cat,η),100_000))
        for (p,q) in ((0,1),(0,2),(1,1),(2,2),(1,3))
            checkmean(conj.(z).^p .* z.^q,tr(ρ*a^q*adjoint(a)^p))
        end
    end
    σ = 0.6
    v = coherentket(β,K)
    for s in (PhaseDiffused(Coherent(β),σ),PhaseDiffused(Coherent(β),σ;common=false))
        z = vec(sampleQ(rng,s,100_000))
        for (p,q) in ((0,1),(0,2),(1,2),(0,3),(2,2))
            checkmean(conj.(z).^p .* z.^q,qmoment(v,a,p,q)*exp(-(p-q)^2*σ^2/2))
        end
    end
    for common in (true,false)
        s = PhaseDiffused(ProductState(Coherent(β),Coherent(im*β)),σ;common)
        Z = sampleQ(rng,s,100_000)
        factor = common ? 1.0 : exp(-σ^2)
        checkmean(conj.(Z[1,:]).*Z[2,:],im*abs2(β)*factor)
    end
    s = ProductState(Fock(3),Fock(2))
    Z = sampleQ(rng,Loss(s,[0.2,0.8]),100_000)
    checkmean(abs2.(Z[1,:]),1.6)
    checkmean(abs2.(Z[2,:]),2.6)
    @test sampleQ(MersenneTwister(2),PhaseDiffused(cat,0),20) == sampleQ(MersenneTwister(2),cat,20)
    @test_throws ArgumentError Loss(cat,1.1)
    @test_throws ArgumentError PhaseDiffused(cat,-1)
    @test_throws ArgumentError PhaseDiffused(s,[0.1,0.2])
    @test_throws DimensionMismatch Loss(s,[0.2])
    @test overlapcapability(Loss(cat,0.4)) isa SamplingOnly
end
