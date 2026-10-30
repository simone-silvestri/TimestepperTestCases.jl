using TimestepperTestCases
using TimestepperTheory
using Documenter
using Literate

DocMeta.setdocmeta!(TimestepperTestCases, :DocTestSetup, :(using TimestepperTestCases); recursive=true)

DocMeta.setdocmeta!(TimestepperTheory, :DocTestSetup, :(using TimestepperTheory); recursive=true)

example = joinpath(@__DIR__, "..", "lib", "TimestepperTheory", "examples", "time_discretization.jl")
results = joinpath(dirname(example), "..", "output")
generated = joinpath(@__DIR__, "src", "generated")
mkpath(generated)

if !("--reuse-figures" in ARGS)
    run(`$(Base.julia_cmd()) --project=$(@__DIR__) $example`)
end

page = Literate.markdown(example, generated; execute=false, flavor=Literate.CommonMarkFlavor(), credit=false)
for image in eachmatch(r"!\[.*?\]\(([^)]+)\)", read(page, String))
    cp(joinpath(results, image.captures[1]), joinpath(generated, image.captures[1]); force=true)
end

experiment_pages = [
    "Internal Tide" => "experiments/internal_tide.md",
    "Channel Flow" => "experiments/channel.md",
]

pages = [
    "Home" => "index.md",
    "Getting Started" => "getting_started.md",
    "Experiments" => experiment_pages,
    "Diagnostics" => "diagnostics.md",
    "Time-discretization theory" => [
        "Overview" => "theory.md",
        "Example guide" => "theory_example.md",
        "Manuscript example" => "generated/time_discretization.md",
        "Theory API" => "theory_api.md",
    ],
    "Stability Analysis" => "stability.md",
    "Notebooks" => "notebooks.md",
    "API Reference" => "api_reference.md",
]

makedocs(;
    modules=[TimestepperTestCases, TimestepperTestCases.BuoyancyVarianceDissipationComputations, TimestepperTheory],
    authors="Simone Silvestri <silvestri.simone0@gmail.com> and contributors",
    sitename="TimestepperTestCases.jl",
    format=Documenter.HTML(;
        canonical="https://simone-silvestri.github.io/TimestepperTestCases.jl",
        edit_link="main",
        prettyurls=get(ENV, "CI", "false") == "true",
        assets=String[],
    ),
    pages=pages,
    checkdocs=:exports,
    warnonly=[:missing_docs, :cross_references],
)
