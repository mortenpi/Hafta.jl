using Documenter, Hafta

makedocs(
    modules = [Hafta],
    format  = :html,
    sitename = "Hafta.jl",
    pages = [
        "Introduction" => "index.md",
        "api.md",
    ],
    html_prettyurls = true,
)

deploydocs(
    repo = "github.com/mortenpi/Hafta.jl.git",
    target = "build",
    julia = "0.6",
    make = nothing,
    deps = nothing,
)
