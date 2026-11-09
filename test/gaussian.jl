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
    z,_ = wigner(rng,Gaussian(Coherent(0)),100_000)
    checkmean((sqrt(2)*real.(z)).^2,0.5)
    checkmean((2real.(z)).^2,1)
    @test_throws ArgumentError Gaussian(0,[0.1 0;0 0.1])
    @test_throws ArgumentError Gaussian(0,[1.0 0.2;0 1.0])
    @test_throws ArgumentError Gaussian(0,[Inf 0;0 1])
    @test_throws DimensionMismatch Gaussian([0,0],Matrix(I,2,2))
    classical = Gaussian(0,[0.5 0;0 1.2])
    z,_ = glauberP(rng,classical,50_000)
    @test all(iszero,real.(z))
    checkmean(abs2.(z),0.35)
    @test_throws ArgumentError glauberP(Gaussian(Squeezed(0,0.3im)),1)
    states = (Gaussian(Coherent(0.2im)),SqueezedThermal(0.3im,0.4cis(0.8),0.7),
              TwoModeSqueezedThermal(0.4,0.7,0.3,0.6),Gaussian(Bogoliubov(sqrt(2)*cis(0.2),cis(0.7),4)))
    for s in states
        M = nmodes(s)
        for (sampler,order) in ((wigner,0.0),(husimiQ,0.5))
            result = sampler(rng,s,80_000)
            Z = reshape(first(result),M,:)
            R = zeros(2M,size(Z,2))
            R[1:2:end,:] .= sqrt(2)*real.(Z.-s.β)
            R[2:2:end,:] .= sqrt(2)*imag.(Z.-s.β)
            for i in 1:2M,j in 1:2M
                checkmean(R[i,:].*R[j,:],s.V[i,j]+order*(i==j))
            end
        end
        A,B = positiveP(rng,s,100_000)
        A,B = reshape(A,M,:),reshape(B,M,:)
        V = s.V
        for i in 1:M,j in 1:M
            n = (V[2i-1,2j-1]+V[2i,2j]+im*(V[2i-1,2j]-V[2i,2j-1]))/2-0.5*(i==j)
            m = (V[2i-1,2j-1]-V[2i,2j]+im*(V[2i-1,2j]+V[2i,2j-1]))/2
            checkmean((B[i,:].-conj(s.β[i])).*(A[j,:].-s.β[j]),n)
            checkmean((A[i,:].-s.β[i]).*(A[j,:].-s.β[j]),m)
        end
        # Passing the [X,Y]=2i covariance divided by two gives identical states.
        VXY = 2s.V
        @test Gaussian(s.β,VXY/2).V == s.V
        @test size(sampleQ(rng,s,0)) == (M,0)
    end
    r,φ = 0.4,0.7
    aa,ap,bb,bp = wigner(rng,SqueezedTwoMode(r,φ),100_000)
    checkmean(aa.*bb,cis(2φ)*sinh(r)*cosh(r))
    checkmean(ap.*aa.-0.5,sinh(r)^2)
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
    for nbar in (0.0,0.4),m in (0,1,2)
        g = SqueezedThermal(0,ϵ,nbar)
        s = PhotonAddedGaussian(g,m)
        ρ,a = added_gaussian_density(96,ϵ,nbar,m)
        ρsmall,asmall = added_gaussian_density(64,ϵ,nbar,m)
        @test abs(tr(ρ*a^2*adjoint(a)^2)-tr(ρsmall*asmall^2*adjoint(asmall)^2))<1e-7
        z = vec(sampleQ(rng,s,100_000))
        for (p,q) in ((0,1),(0,2),(1,1),(1,2),(2,2))
            checkmean(conj.(z).^p .* z.^q,tr(ρ*a^q*adjoint(a)^p))
        end
        if m==1
            v = eigvals(Symmetric(g.V+I/2))
            @test s.probabilities ≈ reverse(v)/sum(v)
        end
        @test size(sampleQ(rng,s,0)) == (1,0)
    end
    @test_throws ArgumentError PhotonAddedGaussian(Gaussian(Coherent(0.1)),1)
    @test_throws ArgumentError PhotonAddedGaussian(TwoModeSqueezedThermal(0.2,0),1)
    @test_throws ArgumentError PhotonAddedGaussian(Gaussian(Coherent(0)),-1)
end
