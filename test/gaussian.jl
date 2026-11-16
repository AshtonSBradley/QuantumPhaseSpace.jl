@testset "Gaussian covariances and quadrature conventions" begin
    rng = MersenneTwister(606)
    # [x,p]=i and [X,Y]=2i describe the same a, with X=√2x,Y=√2p.
    aop = annihilator(8)
    xop,pop = (aop+aop')/sqrt(2),(aop-aop')/(im*sqrt(2))
    Xop,Yop = sqrt(2)*xop,sqrt(2)*pop
    @test (xop+im*pop)/sqrt(2) ≈ aop
    @test (Xop+im*Yop)/2 ≈ aop
    @test (xop*pop-pop*xop)[1:7,1:7] ≈ im*Matrix(I,7,7)
    @test (Xop*Yop-Yop*Xop)[1:7,1:7] ≈ 2im*Matrix(I,7,7)
    z,_ = wigner(rng,Gaussian(Coherent(0)),40_000)
    checkmean((sqrt(2)*real.(z)).^2,0.5)
    @test_throws ArgumentError Gaussian(0,[0.1 0;0 0.1])
    @test_throws ArgumentError Gaussian(0,[1.0 0.2;0 1.0])
    @test_throws ArgumentError Gaussian(0,[Inf 0;0 1])
    @test_throws DimensionMismatch Gaussian([0,0],Matrix(I,2,2))
    classical = Gaussian(0,[0.5 0;0 1.2])
    z,_ = glauberP(rng,classical,30_000)
    @test all(iszero,real.(z))
    checkmean(abs2.(z),0.35)
    @test_throws ArgumentError glauberP(Gaussian(Squeezed(0,0.3im)),1)
    r, φ = 0.4, 0.7
    N = 45_000
    s = TwoModeSqueezedThermal(r,φ,0.3,0.6)
    occupation1 = 0.3 + 1.9*sinh(r)^2
    occupation2 = 0.6 + 1.9*sinh(r)^2
    pair = 1.9*cis(2φ)*sinh(r)*cosh(r)
    for (sampler, offset) in ((wigner,0.5),(husimiQ,1),(positiveP,0))
        a, ap = sampler(rng,s,N)
        checkmean(ap[1,:].*a[1,:].-offset,occupation1)
        checkmean(ap[2,:].*a[2,:].-offset,occupation2)
        checkmean(a[1,:].*a[2,:],pair)
    end

    legacy = SqueezedTwoMode(r,φ)
    for (sampler, offset) in ((wigner,0.5),(husimiQ,1),(positiveP,0))
        a, ap, b, bp = sampler(rng,legacy,N)
        checkmean(ap.*a.-offset,sinh(r)^2)
        checkmean(a.*b,cis(2φ)*sinh(r)*cosh(r))
    end

    # An oriented single-mode covariance tests displacement and quadrature units.
    β, θ, nbar = 0.3im, 0.8, 0.7
    s = SqueezedThermal(β,r*cis(θ),nbar)
    z,_ = husimiQ(rng,s,N)
    checkmean(z,β)
    checkmean((z.-β).^2,-(2nbar+1)*cis(θ)*sinh(r)*cosh(r))
    checkmean(abs2.(z.-β),nbar+(2nbar+1)*sinh(r)^2+1)
    @test Gaussian(s.β,(2s.V)/2).V == s.V

    u,v = sqrt(2)*cis(0.2),cis(0.7)
    s = Gaussian(Bogoliubov(u,v,4))
    z,_ = wigner(rng,s,N)
    checkmean(abs2.(z),13.5)
    checkmean(z.^2,9u*conj(v))
    @test_throws ArgumentError GaussianBogoliubov(ones(2,2),zeros(2,2),0)
end

function added_gaussian_density(K,ϵ,nbar,m)
    a = annihilator(K)
    S = exp((conj(ϵ)*a^2-ϵ*adjoint(a)^2)/2)
    p = [nbar^n/(nbar+1)^(n+1) for n in 0:K-1]
    ρ = S*Diagonal(p)*S'
    ρ = adjoint(a)^m*ρ*a^m
    return ρ/tr(ρ),a
end
@testset "Photon-added centred Gaussian" begin
    rng = MersenneTwister(607)
    ϵ = 0.35cis(0.7)
    for (nbar,m) in ((0.0,0),(0.0,1),(0.4,2))
        g = SqueezedThermal(0,ϵ,nbar)
        s = PhotonAddedGaussian(g,m)
        ρ,a = added_gaussian_density(96,ϵ,nbar,m)
        ρsmall,asmall = added_gaussian_density(64,ϵ,nbar,m)
        @test abs(tr(ρ*a^2*adjoint(a)^2)-tr(ρsmall*asmall^2*adjoint(asmall)^2))<1e-7
        z = vec(sampleQ(rng,s,45_000))
        for (p,q) in ((0,2),(1,1),(2,2))
            checkmean(conj.(z).^p .* z.^q,tr(ρ*a^q*adjoint(a)^p))
        end
        if m==1
            v = eigvals(Symmetric(g.V+I/2))
            @test s.probabilities ≈ reverse(v)/sum(v)
        end
    end
    @test_throws ArgumentError PhotonAddedGaussian(Gaussian(Coherent(0.1)),1)
    @test_throws ArgumentError PhotonAddedGaussian(TwoModeSqueezedThermal(0.2,0),1)
    @test_throws ArgumentError PhotonAddedGaussian(Gaussian(Coherent(0)),-1)
end
