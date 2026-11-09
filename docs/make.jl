using Documenter, PhaseSpaceTools

makedocs(sitename="PhaseSpaceTools.jl",
pages = [
"Overview" => "index.md",
"Introduction to Phase Space" => "phasespace.md",
"Quantum States" => "states.md",
"Distributions" => "distributions.md",
"Exact Q Sampling" => "exact-sampling.md",
"Solving SDES" => "sdes.md"
],
format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true")
)

deploydocs(
    repo = "github.com/AshtonSBradley/PhaseSpaceTools.jl.git",
    devbranch = "master"                  # match your default branch
)

