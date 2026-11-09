module CatQShowcase

using PhaseSpaceTools, CairoMakie, Random, Statistics

# Independent analytic densities, used only for figure contours and validation.
coherent_overlap(z,α) = exp(-abs2(z-α)/2+im*imag(conj(z)*α))
function cat_density(z,α,parity)
    parity in (-1,1) || throw(ArgumentError("parity must be ±1."))
    amplitude = coherent_overlap(z,α)+parity*coherent_overlap(z,-α)
    return abs2(amplitude)/(2π*(1+parity*exp(-2abs2(α))))
end
mixture_density(z,α) = (exp(-abs2(z-α))+exp(-abs2(z+α)))/(2π)

function square_probability(density,halfwidth,bins)
    step = 2halfwidth/bins
    centres = range(-halfwidth+step/2;step,length=bins)
    return sum(density(complex(x,y)) for x in centres,y in centres)*step^2
end

function density_histogram(z;extent=3.6,bins=80)
    step = 2extent/bins
    centres = range(-extent+step/2;step,length=bins)
    counts = zeros(Int,bins,bins)
    captured = 0
    for point in z
        x,y = real(point),imag(point)
        if -extent<=x<extent && -extent<=y<extent
            i,j = floor(Int,(x+extent)/step)+1,floor(Int,(y+extent)/step)+1
            counts[i,j] += 1
            captured += 1
        end
    end
    # Normalise by ALL samples, retaining the mass outside the displayed area.
    density = counts./(length(z)*step^2)
    @assert isapprox(sum(density)*step^2,captured/length(z);atol=1e-14)
    return centres,density,captured/length(z)
end

function check_complex_mean(values,target)
    for part in (real,imag)
        discrepancy = abs(mean(part,values)-part(target))
        @assert discrepancy<=6std(part.(values))/sqrt(length(values))+1e-10
    end
end

function validate_samples(z,density,α,occupation,stats,expected_acceptance)
    N = length(z)
    check_complex_mean(z,0)
    check_complex_mean(z.^2,α^2)
    number = abs2.(z).-1
    check_complex_mean(number,occupation)
    # Numerical integration independently checks density normalisation and a
    # local probability in the interference region, not just global moments.
    mass = square_probability(density,7.0,400)
    @assert abs(mass-1)<1e-9
    central = square_probability(density,0.35,400)
    coarse = square_probability(density,0.35,200)
    @assert abs(central-coarse)<1e-6
    observed = count(t->abs(real(t))<0.35 && abs(imag(t))<0.35,z)/N
    central_se = sqrt(central*(1-central)/N)
    @assert abs(observed-central)<6central_se+1e-6
    rate = acceptance_rate(stats)
    rate_se = sqrt(expected_acceptance*(1-expected_acceptance)/stats.proposals)
    @assert abs(rate-expected_acceptance)<=6rate_se+1e-12
    return (;occupation=mean(number),occupation_target=occupation,
        occupation_se=std(number)/sqrt(N),central_probability=observed,
        central_target=central,central_se,acceptance=rate,
        acceptance_target=expected_acceptance,proposals=stats.proposals,
        analytic_mass=mass)
end

function make_comparison(α,N,seed,pathbase,amplitude_label)
    A = abs2(α)
    cases = (
        ("Even cat",EvenCat(α),z->cat_density(z,α,1),A*tanh(A),(1+exp(-2A))/2),
        ("Odd cat",OddCat(α),z->cat_density(z,α,-1),A/tanh(A),(1-exp(-2A))/2),
        ("Classical mixture",Mixture((Coherent(α),Coherent(-α)),[0.5,0.5]),
            z->mixture_density(z,α),A,1.0))
    rng = MersenneTwister(seed)
    checks = NamedTuple[]
    figure = with_theme(theme_latexfonts()) do
        fig = Figure(size=(1000,440),fontsize=20,figure_padding=18)
        sample_label = replace(string(N),r"(?<=\d)(?=(\d{3})+$)"=>",")
        Label(fig[0,1:4],"$amplitude_label   ·   $sample_label Q samples per state";
            fontsize=20,tellwidth=false)
        fine = range(-3.6,3.6;length=401)
        heat = nothing
        for (j,(name,state,density,occupation,expected_acceptance)) in enumerate(cases)
            Z,stats = sampleQ(rng,state,N;diagnostics=true)
            z = vec(Z)
            check = validate_samples(z,density,α,occupation,stats,expected_acceptance)
            centres,hist,captured = density_histogram(z)
            push!(checks,merge((;state=name,alpha_real=real(α),alpha_imag=imag(α),N),
                check,(;displayed_mass=captured)))
            ax = Axis(fig[1,j];title=name,xlabel="Re z",ylabel=j==1 ? "Im z" : "",
                aspect=DataAspect(),xticks=[-2,0,2],yticks=[-2,0,2],
                xticklabelsize=17,yticklabelsize=17,titlesize=21,
                xgridvisible=false,ygridvisible=false)
            heat = heatmap!(ax,centres,centres,hist;colormap=:Blues,colorrange=(0,0.2))
            exact = [density(complex(x,y)) for x in fine,y in fine]
            contour!(ax,fine,fine,exact;levels=[0.02,0.06,0.12,0.16],
                color=(:black,0.65),linewidth=1)
            scatter!(ax,[-real(α),real(α)],[-imag(α),imag(α)];
                marker=:xcross,markersize=14,color=:darkorange3)
            limits!(ax,-3.6,3.6,-3.6,3.6)
            hidespines!(ax,:t,:r)
        end
        Colorbar(fig[1,4],heat;label="Q(z)",width=14,ticks=[0,0.1,0.2],
            ticklabelsize=17,labelsize=20,spinewidth=0.5)
        Legend(fig[2,1:4],
            [LineElement(color=(:black,0.65),linewidth=1),
             MarkerElement(marker=:xcross,color=:darkorange3,markersize=14)],
            ["Exact Q contours","Coherent centres ±α"];
            orientation=:horizontal,framevisible=false,labelsize=18)
        colgap!(fig.layout,12)
        rowgap!(fig.layout,8)
        fig
    end
    save(pathbase*".png",figure;px_per_unit=2)
    save(pathbase*".pdf",figure)
    return checks
end

# A finite-Fock calculation is an independent check of the analytic line cuts.
# It is used only for validation, never for sampling or drawing the Q density.
function coherent_coefficients(β,K)
    v = zeros(ComplexF64,K)
    v[1] = exp(-abs2(β)/2)
    for n in 1:K-1
        v[n+1] = β*v[n]/sqrt(n)
    end
    return v
end

"""Show weak large-cat Q interference using logarithmic colour and a central cut."""
function make_large_cat_interference(;output=joinpath(@__DIR__,"..","docs","src","assets"))
    mkpath(output)
    α = 3.0
    A,d = α^2,exp(-2α^2)
    reference = exp(-A)/π
    ys = range(-2.5,2.5;length=2001)
    plus = [cat_density(im*y,α,1)/reference for y in ys]
    minus = [cat_density(im*y,α,-1)/reference for y in ys]
    mixed = [mixture_density(im*y,α)/reference for y in ys]
    plus_target = @. 2exp(-ys^2)*cos(α*ys)^2/(1+d)
    minus_target = @. 2exp(-ys^2)*sin(α*ys)^2/(1-d)
    @assert maximum(abs.(plus .- plus_target))<1e-12
    @assert maximum(abs.(minus .- minus_target))<1e-12
    @assert maximum(abs.(mixed .- exp.(-ys .^ 2)))<1e-12
    # Check the exact nodes and the fringe period of the actual positive Q law.
    @assert cat_density(0,α,-1)==0
    @assert cat_density(im*π/(2α),α,1)<1e-30
    @assert isapprox(cat_density(0,α,1)/reference,2/(1+d);rtol=1e-13)
    K = 80
    tail = ccdf(Poisson(A),K-1)
    @assert tail<1e-40
    v = coherent_coefficients(α,K)
    even,odd = copy(v),copy(v)
    even[2:2:end] .= 0
    odd[1:2:end] .= 0
    even ./= norm(even)
    odd ./= norm(odd)
    max_fock_error = 0.0
    for y in range(-2.5,2.5;length=101), (state,parity) in ((even,1),(odd,-1))
        coherent = coherent_coefficients(im*y,K)
        fock_q = abs2(dot(coherent,state))/π
        exact_q = cat_density(im*y,α,parity)
        @assert isapprox(fock_q,exact_q;rtol=1e-10,atol=1e-18)
        max_fock_error = max(max_fock_error,abs(fock_q-exact_q))
    end
    figure = with_theme(theme_latexfonts()) do
        fig = Figure(size=(1000,690),fontsize=20,figure_padding=20)
        Label(fig[0,1:4],"α = 3   ·   Analytic Husimi–Q densities";fontsize=21,tellwidth=false)
        xs,imaginary = range(-4.5,4.5;length=721),range(-2.3,2.3;length=369)
        cases = (("Even cat",z->cat_density(z,α,1)),
                 ("Odd cat",z->cat_density(z,α,-1)),
                 ("Classical mixture",z->mixture_density(z,α)))
        heat = nothing
        for (j,(name,density)) in enumerate(cases)
            ax = Axis(fig[1,j];title=name,xlabel="Re z",ylabel=j==1 ? "Im z" : "",
                aspect=DataAspect(),xticks=[-3,0,3],yticks=[-2,0,2],
                xgridvisible=false,ygridvisible=false,xticklabelsize=17,yticklabelsize=17)
            values = [density(complex(x,y)) for x in xs,y in imaginary]
            @assert minimum(values)>=0
            heat = heatmap!(ax,xs,imaginary,max.(values,1e-9);
                colormap=:Blues,colorscale=log10,colorrange=(1e-9,0.2),rasterize=2)
            scatter!(ax,[-α,α],[0,0];marker=:xcross,markersize=12,color=:darkorange3)
            limits!(ax,-4.5,4.5,-2.3,2.3)
            hidespines!(ax,:t,:r)
        end
        Colorbar(fig[1,4],heat;label="Q(z), log scale",width=14,
            ticks=([1e-8,1e-5,1e-2,1e-1],["10⁻⁸","10⁻⁵","10⁻²","10⁻¹"]),
            ticklabelsize=17,labelsize=18,spinewidth=0.5)
        cut = Axis(fig[2,1:4];title="Central overlap cut: Re z = 0",xlabel="Im z",
            ylabel="Q(i y) / Qmix(0)",xticks=-2:1:2,yticks=[0,1,2],
            xgridvisible=false,ygridvisible=false,xticklabelsize=17,yticklabelsize=17)
        lines!(cut,ys,plus;color=:dodgerblue3,linewidth=2.5,label="Even cat")
        lines!(cut,ys,minus;color=:darkorange3,linewidth=2.5,label="Odd cat")
        lines!(cut,ys,mixed;color=:gray40,linewidth=2,linestyle=:dash,label="Classical mixture")
        limits!(cut,-2.5,2.5,0,2.13)
        hidespines!(cut,:t,:r)
        Legend(fig[3,1:4],cut;orientation=:horizontal,framevisible=false,labelsize=18)
        Label(fig[4,1:4],"Common line-cut scale: Qmix(0) = exp(−9)/π ≈ 3.93 × 10⁻⁵";
            fontsize=18,tellwidth=false)
        colgap!(fig.layout,12)
        rowgap!(fig.layout,10)
        fig
    end
    save(joinpath(output,"cat-q-interference.png"),figure;px_per_unit=2)
    save(joinpath(output,"cat-q-interference.pdf"),figure)
    open(joinpath(output,"cat-q-interference-cut.tsv"),"w") do io
        println(io,"imag_z\tQ_even\tQ_odd\tQ_mixture")
        for y in range(-2.5,2.5;length=401)
            println(io,join((y,cat_density(im*y,α,1),cat_density(im*y,α,-1),mixture_density(im*y,α)),'\t'))
        end
    end
    checks = (;alpha=α,fringe_period=π/α,mixture_origin=reference,
        even_origin=cat_density(0,α,1),fock_cutoff=K,poisson_tail=tail,max_fock_error)
    println("Large-cat Q validation: ",checks)
    return checks
end

"""
    make_cat_q_figures(; N=300_000, output=...)

Regenerate the sampled comparisons, large-cat interference view and validation data. No
finite-window sampling is used. The plotted histograms keep their full-batch
normalisation, and the record reports the fraction inside the display window.
"""
function make_cat_q_figures(;N::Integer=300_000,
        output=joinpath(@__DIR__,"..","docs","src","assets"))
    N>=2 || throw(ArgumentError("figure validation needs N≥2."))
    mkpath(output)
    main = make_comparison(1.2+0im,N,4242,joinpath(output,"cat-q"),"α = 1.2")
    rotated = make_comparison(1.2cis(π/6),N,4243,
        joinpath(output,"cat-q-rotated"),"α = 1.2 exp(iπ/6)")
    checks = vcat(main,rotated)
    record = joinpath(output,"cat-q-validation.tsv")
    open(record,"w") do io
        println(io,join(string.(keys(first(checks))),'\t'))
        for row in checks
            println(io,join(values(row),'\t'))
        end
    end
    make_large_cat_interference(;output)
    println("Saved cat Q comparisons and validation: ",record)
    return checks
end

end # module

if abspath(PROGRAM_FILE)==@__FILE__
    CatQShowcase.make_cat_q_figures()
end
