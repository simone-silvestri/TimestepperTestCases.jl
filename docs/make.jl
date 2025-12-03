using TimestepperTestCases
using Documenter

DocMeta.setdocmeta!(TimestepperTestCases, :DocTestSetup, :(using TimestepperTestCases); recursive=true)

experiment_pages = [
    "Internal Tide" => "experiments/internal_tide.md",
    "Idealized Coast" => "experiments/idealized_coast.md",
    "Channel Flow" => "experiments/channel.md",
]

pages = [
    "Home" => "index.md",
    "Getting Started" => "getting_started.md",
    "Experiments" => experiment_pages,
    "Diagnostics" => "diagnostics.md",
    "Stability Analysis" => "stability.md",
    "Notebooks" => "notebooks.md",
    "API Reference" => "api_reference.md",
]

makedocs(;
    modules=[TimestepperTestCases],
    authors="Simone Silvestri <silvestri.simone0@gmail.com> and contributors",
    sitename="TimestepperTestCases.jl",
    format=Documenter.HTML(;
        canonical="https://simone-silvestri.github.io/TimestepperTestCases.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=pages,
)

deploydocs(;
    repo="github.com/simone-silvestri/TimestepperTestCases.jl",
    devbranch="main",
)
