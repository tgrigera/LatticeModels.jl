using LatticeModels
using Documenter

DocMeta.setdocmeta!(LatticeModels, :DocTestSetup, :(using LatticeModels); recursive=true)

makedocs(;
    sitename="LatticeModels.jl",
    modules=[LatticeModels],
    authors="Tomás S. Grigera",
    #repo="https://github.com/tgrigera/LatticeModels.jl/blob/{commit}{path}#{line}",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://tgrigera.github.io/LatticeModels.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Graphs" => "graphs.md",
        "Models" => "models.md"
    ],
)

if get(ENV, "CI", nothing) == "true" 
    deploydocs( repo = "github.com/tgrigera/LatticeModels.jl.git",
                    devbranch="main")
end
