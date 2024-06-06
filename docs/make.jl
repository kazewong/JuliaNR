using JuliaNR
using Documenter

DocMeta.setdocmeta!(JuliaNR, :DocTestSetup, :(using JuliaNR); recursive=true)

makedocs(;
    modules=[JuliaNR],
    authors="Kaze Wong <kazewong.physics@gmail.com> and contributors",
    sitename="JuliaNR.jl",
    format=Documenter.HTML(;
        canonical="https://kazewong.github.io/JuliaNR.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/kazewong/JuliaNR.jl",
    devbranch="main",
)
