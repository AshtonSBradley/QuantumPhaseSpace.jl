using Documenter, QuantumPhaseSpace

makedocs(sitename="QuantumPhaseSpace.jl",
pages = [
"Overview" => "index.md",
"Introduction to Phase Space" => "phasespace.md",
"Quantum States" => "states.md",
"Distributions" => "distributions.md",
"Exact Q Sampling" => "exact-sampling.md",
"Example: Cat States" => "cat-q.md",
"Examples: Fock-Basis States" => "fock-basis.md",
"Fermionic Truncated Wigner" => "fermions.md",
"Support and Examples" => "support.md",
"Solving SDES" => "sdes.md"
],
format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true")
)

deploydocs(
    repo = "github.com/AshtonSBradley/QuantumPhaseSpace.jl.git",
    devbranch = "master",
    versions = ["stable" => "v^", "v#.#", "dev" => "dev"],
)
