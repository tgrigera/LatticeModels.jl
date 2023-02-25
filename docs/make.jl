using LatticeModels
using Documenter

DocMeta.setdocmeta!(LatticeModels, :DocTestSetup, :(using LatticeModels); recursive=true)

makedocs(;
    modules=[LatticeModels],
    authors="Tomas S. Grigera <tgrigera@iflysib.unlp.edu.ar> and contributors",
    repo="https://github.com/tgrigera/LatticeModels.jl/blob/{commit}{path}#{line}",
    sitename="LatticeModels.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://tgrigera.github.io/LatticeModels.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/tgrigera/LatticeModels.jl",
    devbranch="main",
)
