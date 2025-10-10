using TimestepperTestCases
using Documenter

DocMeta.setdocmeta!(TimestepperTestCases, :DocTestSetup, :(using TimestepperTestCases); recursive=true)

makedocs(;
    modules=[TimestepperTestCases],
    authors="Simone Silvestri <silvestri.simone0@gmail.com> and contributors",
    sitename="TimestepperTestCases.jl",
    format=Documenter.HTML(;
        canonical="https://simone-silvestri.github.io/TimestepperTestCases.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/simone-silvestri/TimestepperTestCases.jl",
    devbranch="main",
)
