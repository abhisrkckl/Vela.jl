using Documenter, Vela

makedocs(
    sitename = "Vela.jl",
    authors = "Abhimanyu Susobhanan",
    pages = [
        "Home" => "index.md",
        "getting-started.md",
        "tutorial.md",
        "explanation.md",
        # "api-reference.md",
    ],
    # format = Documenter.LaTeX(platform = "none")
)
