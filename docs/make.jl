using GeneralAttractors
using Documenter

DocMeta.setdocmeta!(
    GeneralAttractors,
    :DocTestSetup,
    :(using GeneralAttractors);
    recursive = true,
)

makedocs(;
    modules = [GeneralAttractors],
    authors = "Federico Claudi, Sarthak Chandran",
    repo = "https://github.com/FedeClaudi/GeneralAttractors.jl/blob/{commit}{path}#{line}",
    sitename = "GeneralAttractors.jl",
    format = Documenter.HTML(;
        prettyurls = get(ENV, "CI", "false") == "true",
        canonical = "https://FedeClaudi.github.io/GeneralAttractors.jl",
        edit_link = "main",
        assets = String[],
    ),
    pages = ["Home" => "index.md"],
)

deploydocs(; repo = "github.com/FedeClaudi/GeneralAttractors.jl", devbranch = "main")
