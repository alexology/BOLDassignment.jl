using Documenter
using BOLDassignment

DocMeta.setdocmeta!(BOLDassignment, :DocTestSetup, :(using BOLDassignment); recursive=true)

makedocs(;
    modules=[BOLDassignment],
    sitename="BOLDassignment.jl",
    authors="Alex Laini",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://alexology.github.io/BOLDassignment.jl",
        edit_link="main",
    ),
    pages=[
        "Home" => "index.md",
        "The LCA assignment" => "lca.md",
        "Top hits as in BOLDigger3" => "boldigger3.md",
        "The BOLD service" => "service.md",
        "API" => "api.md",
    ],
    checkdocs=:exports,
)

deploydocs(; repo="github.com/alexology/BOLDassignment.jl", devbranch="main")
