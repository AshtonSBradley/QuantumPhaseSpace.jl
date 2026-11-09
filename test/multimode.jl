@testset "Passive mixing and collective Fock states" begin
    rng = MersenneTwister(404)
    U = [1 im;im 1]/sqrt(2)
    input = ProductState(Fock(2),Fock(0))
    Z = sampleQ(rng,Passive(input,U),100_000)
    target = Matrix{ComplexF64}(I,2,2)+2U[:,1]*U[:,1]'
    for i in 1:2,j in 1:2
        checkmean(Z[i,:].*conj.(Z[j,:]),target[i,j])
    end
    for u in ([1.0+0im],[sqrt(0.3),sqrt(0.7)*cis(0.8)],[0.5,0.5im,sqrt(0.5)*cis(-0.4)]), n in (0,1,4)
        s = FixedNumberCondensate(n,u)
        M = length(u)
        Z = sampleQ(rng,s,100_000)
        for i in 1:M
            checkmean(Z[i,:],0)
            for j in 1:M
                checkmean(Z[i,:].*conj.(Z[j,:]),(i==j)+n*u[i]*conj(u[j]))
            end
            # E_Q |z_j|⁴ = <n_j(n_j-1)> + 4<n_j> + 2.
            checkmean(abs2.(Z[i,:]).^2,n*(n-1)*abs2(u[i])^2+4n*abs2(u[i])+2)
        end
        for method in (:specialised,:generic)
            a,ap = positiveP(rng,s,100_000;method)
            a,ap = reshape(a,M,:),reshape(ap,M,:)
            for i in 1:M,j in 1:M
                checkmean(ap[i,:].*a[j,:],n*conj(u[i])*u[j])
            end
        end
        @test size(sampleQ(rng,s,0)) == (M,0)
    end
    # Small two-mode Hilbert expansion: coefficients from the binomial theorem.
    n = 3; u = [sqrt(0.4),sqrt(0.6)*cis(0.7)]
    occupations = [n-k for k in 0:n]
    coeffs = [sqrt(binomial(n,k))*u[1]^(n-k)*u[2]^k for k in 0:n]
    expanded = FockSuperposition(vcat(transpose(occupations),transpose(collect(0:n))),coeffs)
    s = FixedNumberCondensate(n,u)
    for z in ([0.3im,0.2],[1.2+0.3im,-0.5im])
        @test exp(logoverlap(s,z)) ≈ exp(logoverlap(expanded,z)) atol=1e-13
    end
    @test_throws ArgumentError FixedNumberCondensate(1,[1,1])
    @test_throws ArgumentError FixedNumberCondensate(-1,[1])
    @test_throws ArgumentError Passive(input,[1 1;0 1])
    @test_throws DimensionMismatch Passive(input,ones(1,1))
    @test all(iszero,first(positiveP(FixedNumberCondensate(0,u),3)))
end
