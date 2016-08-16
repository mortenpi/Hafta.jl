using Documenter, Hafta

makedocs(
    modules = [Hafta],
    format  = Documenter.Formats.HTML,
    sitename = "Hafta.jl",
    pages = [
        "Introduction" => "index.md",
        "api.md",
    ]
)

deploydocs(
    repo = "github.com/JuliaDocs/Documenter.jl.git",
    target = "build",
    make = () -> nothing,
    deps = () -> nothing,
    julia = "0.5",
)
