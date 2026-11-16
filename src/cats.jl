"""
    Cat(α, phase=0)

Normalised `|α⟩+exp(i phase)|-α⟩`. Exact Gaussian rejection is used, except
small odd cats use a Gamma proposal with acceptance tending to one. `OddCat(0)` is
undefined; request `Fock(1)` explicitly for its limiting state. General phases
near destructive cancellation can be inefficient: inspect `sampleQ` diagnostics.
"""
struct Cat <: State
    α::ComplexF64
    phase::Float64
    parity::Int
    D::Float64
    function Cat(α::Number,phase::Real=0.0)
        all(isfinite,(α,phase)) || throw(ArgumentError("cat parameters must be finite."))
        α,phase = ComplexF64(α),rem2pi(Float64(phase),RoundNearest)
        parity = phase == 0 ? 1 : abs(phase) == Float64(π) ? -1 : 0
        A = abs2(α)
        D = parity == -1 ? -expm1(-2A) : parity == 1 ? 1+exp(-2A) :
            -expm1(-2A)+2exp(-2A)*cos(phase/2)^2
        D>0 && isfinite(D) || throw(ArgumentError("zero/undefined cat state; use Fock(1) explicitly for the odd α=0 limit."))
        new(α,phase,parity,D)
    end
end
"""Even-parity cat, using exact parity arithmetic."""
EvenCat(α) = Cat(α,0)
"""Odd-parity cat; α must be nonzero. Small finite amplitudes retain their exact cat state."""
OddCat(α) = Cat(α,π)
overlapcapability(::Cat) = PureOverlap()
function _logcoshsinh(t,sign)
    scale = abs(real(t))
    return scale+log(exp(t-scale)+sign*exp(-t-scale))-log(2)
end
function _logsinhc(t)
    if abs(t)<1e-3
        t2 = t*t
        return log1p(t2/6+t2^2/120+t2^3/5040)
    end
    return _logcoshsinh(t,-1)-log(t)
end
function logoverlap(s::Cat,z::Number)
    t = s.α*conj(z)
    h = s.parity == -1 ? (iszero(t) ? complex(-Inf) : log(t)+_logsinhc(t)) :
        _logcoshsinh(t-im*s.phase/2,1)+im*s.phase/2
    return -(abs2(z)+abs2(s.α))/2+log(2)+h-log(2s.D)/2
end
function _catrates(s)
    A = abs2(s.α)
    direct = s.D/2
    small = s.parity == -1 && 0<A<3 ? (1-A/3)^2 * (s.D/A) * exp(A)/2 : 0.0
    return direct,small
end
function _qfill!(rng,Z,s::Cat,d,limit)
    direct,small = _catrates(s)
    A = abs2(s.α)
    if small > direct
        return _qfill_odd_gamma!(rng,Z,s,d,limit,A/3,0.0)
    end
    return _qfill_cat_gaussian!(rng,Z,s,d,limit)
end

function _qfill_odd_gamma!(rng,Z,s::Cat,d,limit,broadening,offset)
    radial = Gamma(2,inv(1-broadening))
    while d.accepted<size(Z,2)
        _proposal!(d,limit)
        B = rand(rng,radial)
        z = sqrt(B)*cis(2π*rand(rng))
        t = conj(s.α)*z
        p = exp(-broadening*B-offset+2real(_logsinhc(t)))
        if rand(rng)<_acceptprob(p)
            d.accepted += 1
            Z[1,d.accepted] = z
        end
    end
    return Z
end

function _qfill_cat_gaussian!(rng,Z,s::Cat,d,limit)
    while d.accepted<size(Z,2)
        _proposal!(d,limit)
        z = (rand(rng,Bool) ? s.α : -s.α)+randnc(rng)
        t = conj(s.α)*z
        e = exp(-2abs(real(t)))
        sech = 2e/(1+e^2)
        interference = s.parity == 1 ? cos(2imag(t)) : s.parity == -1 ? -cos(2imag(t)) : cos(2imag(t)+s.phase)
        p = (1+interference*sech)/2
        if rand(rng)<_acceptprob(p)
            d.accepted += 1
            Z[1,d.accepted] = z
        end
    end
    return Z
end
_qpoint!(rng,z,s::Cat) = (_qfill!(rng,reshape(z,:,1),s,SamplingDiagnostics(),_attemptlimit(1)); z)
_catcomponents(s::Cat) = ((Coherent(s.α),Coherent(-s.α)),
    (1.0+0im,s.parity == -1 ? -1.0+0im : s.parity == 1 ? 1.0+0im : cis(s.phase)))
_preparation_primitives(s::Cat) = reduce(vcat,(_primitives(x) for x in first(_catcomponents(s))))
function _stateinner(s::Cat,t::State)
    states,c = _catcomponents(s)
    return sum(conj(cj)*_stateinner(x,t) for (cj,x) in zip(c,states))/sqrt(BigFloat(2)*s.D)
end
_stateinner(s::State,t::Cat) = conj(_stateinner(t,s))
function _stateinner(s::Cat,t::Cat)
    xs,c = _catcomponents(s)
    ys,b = _catcomponents(t)
    return sum(conj(cj)*bk*_stateinner(x,y) for (cj,x) in zip(c,xs),(bk,y) in zip(b,ys))/
        sqrt(BigFloat(4)*s.D*t.D)
end
_stateinner(s::Cat,t::Superposition) = conj(_stateinner(t,s))
function _stateinner(s::Superposition,t::Cat)
    return sum(conj(c)*_stateinner(x,t) for (c,x) in zip(s.coefficients,s.states))/sqrt(BigFloat(s.normalisation))
end
