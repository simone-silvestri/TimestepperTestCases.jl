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

figure_directories = (joinpath(@__DIR__, "..", "examples", "output"),
                      joinpath(@__DIR__, "..", "internal_tide", "figures"), joinpath(@__DIR__, "figures"))

for name in ("internal_tide", "internal_tide_figures", "channel", "global_ocean", "global_ocean_results", "global_surface_averages")
    source = joinpath(@__DIR__, "..", "examples", "figures", name * ".jl")
    local page = Literate.markdown(source, generated; execute=false, flavor=Literate.CommonMarkFlavor(), credit=false)
    content = read(page, String)
    for image in eachmatch(r"!\[.*?\]\(([^)]+)\)", content)
        paths = [joinpath(directory, image.captures[1]) for directory in figure_directories]
        found = findfirst(isfile, paths)
        if isnothing(found)
            content = replace(content, image.match => "Run this example with local or published data to generate the figure.")
        else
            cp(paths[found], joinpath(generated, image.captures[1]); force=true)
        end
    end
    write(page, content)
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
    "Data-driven examples" => [
        "Data and execution" => "data_examples.md",
        "Internal tide" => "generated/internal_tide.md",
        "Internal tide: campaign figures" => "generated/internal_tide_figures.md",
        "Channel" => "generated/channel.md",
        "Global cost" => "generated/global_ocean.md",
        "Global diffusivity" => "generated/global_ocean_results.md",
        "Global surface averages" => "generated/global_surface_averages.md",
    ],
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
