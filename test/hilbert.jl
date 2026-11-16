# Independent finite-Fock calculations; these are validation tools, never sampler code.
function annihilator(K)
    a = zeros(ComplexF64,K,K)
    for n in 1:K-1
        a[n,n+1] = sqrt(n)
    end
    return a
end
function coherentket(β,K)
    v = zeros(ComplexF64,K)
    v[1] = exp(-abs2(β)/2)
    for n in 1:K-1
        v[n+1] = β*v[n]/sqrt(n)
    end
    return v
end
qmoment(v,a,p,q) = dot(v,a^q*adjoint(a)^p*v)
function checkqket(rng,s,v;N=45_000)
    a = annihilator(length(v))
    @test abs(norm(v)-1)<1e-10
    z = vec(sampleQ(rng,s,N))
    for (p,q) in ((0,2),(1,2),(2,2))
        checkmean(conj.(z).^p .* z.^q,qmoment(v,a,p,q))
    end
end
