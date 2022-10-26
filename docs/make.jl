using NormalForms
using Documenter

DocMeta.setdocmeta!(NormalForms, :DocTestSetup, :(using NormalForms); recursive=true)

makedocs(;
    modules=[NormalForms],
    authors="Brandon Flores",
    repo="https://github.com/brainandforce/NormalForms.jl/blob/{commit}{path}#{line}",
    sitename="NormalForms.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://brainandforce.github.io/NormalForms.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/brainandforce/NormalForms.jl",
    devbranch="main",
)
