include("hilbert.jl")
@testset "Finite superpositions and cats" begin
    rng = MersenneTwister(301)
    K = 48
    a = annihilator(K)
    α,φ = 0.8+0.3im,0.7
    v = coherentket(α,K)+cis(φ)*coherentket(-α,K)
    normalize!(v)
    s = CoherentSuperposition([α,-α],[1,cis(φ)])
    checkqket(rng,s,v)
    for c in ([2im,2im*cis(φ)],[-1,-cis(φ)]), perm in ([1,2],[2,1])
        t = CoherentSuperposition([α,-α][perm],c[perm])
        @test t.expected_acceptance ≈ s.expected_acceptance
        for z in (0.2im,0.8-0.4im)
            @test abs2(exp(logoverlap(s,z))) ≈ abs2(exp(logoverlap(t,z)))
        end
    end
    @test_throws ArgumentError Superposition((Thermal(0,1),),[1])
    @test_throws ArgumentError Superposition((Fock(0),Fock(0)),[1,-1])
    @test_throws DimensionMismatch Superposition((Fock(0),ProductState(Fock(0),Fock(0))),[1,1])
    @test_throws SamplingError sampleQ(s,10;max_attempts=0)
    for t in (s,NOON(1),NOON(20),Compass(α),FockSuperposition([0,2,5],[1,im,0.3]))
        Z,d = sampleQ(rng,t,15_000;diagnostics=true)
        p = t.expected_acceptance
        @test abs(acceptance_rate(d)-p)<6sqrt(p*(1-p)/d.proposals)
        @test size(Z) == (nmodes(t),15_000)
    end
    # Sparse bosonic-code-like finite sum, with interference between |0> and |2>.
    v = zeros(ComplexF64,K); v[[1,3,6]] .= [1,im,0.3]; normalize!(v)
    checkqket(rng,FockSuperposition([0,2,5],[1,im,0.3]),v)
    β = 0.4+0.5im
    for m in (0,1,3)
        v = adjoint(a)^m*coherentket(β,K); normalize!(v)
        checkqket(rng,PhotonAddedCoherent(m,β),v)
    end
    @test_throws ArgumentError OddCat(0)
    for amplitude in (0.000001+0.000002im,0.3im,α), phase in (0,π,φ)
        t = Cat(amplitude,phase)
        Z,d = sampleQ(rng,t,40_000;diagnostics=true)
        z = vec(Z)
        A,δ,D = abs2(amplitude),exp(-2abs2(amplitude)),t.D
        μ = t.parity == 0 ? -im*amplitude*δ*sin(phase)/D : 0.0im
        checkmean(z,μ)
        checkmean(z.^2,amplitude^2)
        checkmean(abs2.(z),1+A*(1-δ*(t.parity == -1 ? -1 : cos(phase)))/D)
        p = maximum(PhaseSpaceTools._catrates(t))
        @test abs(acceptance_rate(d)-p)<6sqrt(p*(1-p)/d.proposals)
    end
    checkqket(rng,Cat(α,φ),(coherentket(α,K)+cis(φ)*coherentket(-α,K))/norm(coherentket(α,K)+cis(φ)*coherentket(-α,K)))
    # Complex displaced squeezed overlaps: compare with matrix exponentials.
    ϵ = 0.35cis(0.8)
    S = exp((conj(ϵ)*a^2-ϵ*adjoint(a)^2)/2)
    D = exp(β*adjoint(a)-conj(β)*a)
    v = D*S[:,1]
    for z in (0.2+0.1im,0.7-0.4im)
        @test exp(logoverlap(Squeezed(β,ϵ),z)) ≈ dot(coherentket(z,K),v) atol=1e-12
    end
    t = SqueezedCat(β,ϵ;phase=0.4)
    vm = exp(-β*adjoint(a)+conj(β)*a)*S[:,1]
    vs = v+cis(0.4)*vm; normalize!(vs)
    checkqket(rng,t,vs)
    @test t.normalisation ≈ norm(v+cis(0.4)*vm)^2/4 atol=1e-12
    # Displaced Fock Gram and phase, independently from exp(D).
    states = (Displaced(Fock(2),β),Displaced(Fock(3),-0.3im),Squeezed(β,ϵ))
    v1 = D[:,3]; v2 = exp(-0.3im*adjoint(a)-0.3im*a)[:,4]
    t = Superposition(states,[1,0.3im,-0.4])
    cols = hcat(v1,v2,v)
    @test t.gram ≈ cols'*cols atol=2e-12
    vs = cols*[1,0.3im,-0.4]; normalize!(vs)
    checkqket(rng,t,vs)
    # NOON phase-sensitive coherence at N=1 and N=2.
    for n in (1,2)
        Z = sampleQ(rng,NOON(n;phase=φ),100_000)
        checkmean(conj.(Z[1,:]).^n .* Z[2,:].^n,factorial(n)*cis(φ)/2)
    end
    @test real(logoverlap(Fock(10000),sqrt(10000))) > -10
    @test isfinite(real(logoverlap(CoherentSuperposition([40,-40],[1,1]),40)))
    @test size(sampleQ(rng,s,0)) == (1,0)
end
