using Documenter, Vela

makedocs(
    sitename = "Vela.jl",
    authors = "Abhimanyu Susobhanan",
    pages = [
        "Home" => "index.md",
        "getting-started.md",
        "tutorial.md",
        "Explanation & API Reference" =>
            ["precision.md", "quantities.md", "toas.md", "timing-model.md", "residuals.md"],
        # "api-reference.md",
    ],
    # format = Documenter.LaTeX(platform = "none")
)
