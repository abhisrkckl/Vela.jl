using Documenter, Vela

makedocs(
    sitename = "Vela.jl",
    authors = "Abhimanyu Susobhanan",
    pages = [
        "Home" => "index.md",
        "installation.md",
        "pyvela.md",
        "Explanation & API Reference" => [
            "precision.md",
            "quantities.md",
            "toas.md",
            "timing-model.md",
            "residuals.md",
            "likelihood.md",
            "parameters.md",
            "priors.md",
            "red-noise.md",
        ],
    ],
    # format = Documenter.LaTeX(platform = "none")
)

deploydocs(repo = "github.com/abhisrkckl/Vela.jl.git")
