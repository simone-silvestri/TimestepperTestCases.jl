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

weights_example = joinpath(dirname(example), "reconstruction_weights.jl")
Literate.markdown(weights_example, generated; execute=true, credit=false)

figure_directories = (joinpath(@__DIR__, "..", "examples", "output"), joinpath(@__DIR__, "figures"))
on_ci = get(ENV, "CI", "false") == "true"

for name in ("internal_tide_figures", "channel", "global_ocean_results", "global_surface_results", "global_surface_averages", "global_ocean")
    source = joinpath(@__DIR__, "..", "examples", "figures", name * ".jl")
    local page = Literate.markdown(source, generated; execute=false, flavor=Literate.CommonMarkFlavor(), credit=false)
    content = read(page, String)
    for image in eachmatch(r"!\[.*?\]\(([^)]+)\)", content)
        paths = [joinpath(directory, image.captures[1]) for directory in figure_directories]
        found = findfirst(isfile, paths)
        if isnothing(found)
            on_ci && error("$(image.captures[1]), shown in examples/figures/$name.jl, is missing from docs/figures.")
            @warn "Figure missing; the page shows a placeholder" figure=image.captures[1] example=name
            content = replace(content, image.match => "Run this example with local or published data to generate the figure.")
        else
            cp(paths[found], joinpath(generated, image.captures[1]); force=true)
        end
    end
    write(page, content)
end

pages = [
    "Home" => "index.md",
    "Getting started" => "getting_started.md",
    "Experiments" => [
        "Internal tide" => [
            "Configuration and run" => "experiments/internal_tide.md",
            "Manuscript figures" => "generated/internal_tide_figures.md",
        ],
        "Re-entrant channel" => [
            "Configuration and run" => "experiments/channel.md",
            "Manuscript figures" => "generated/channel.md",
        ],
        "Quarter-degree global ocean" => [
            "Configuration and run" => "experiments/global_ocean.md",
            "Numerical diffusivity" => "generated/global_ocean_results.md",
            "Surface velocity" => "generated/global_surface_results.md",
            "Surface tracers" => "generated/global_surface_averages.md",
            "Computational cost" => "generated/global_ocean.md",
        ],
    ],
    "Reproducing the paper" => "data_examples.md",
    "Variance-budget diagnostics" => "diagnostics.md",
    "Time-discretization theory" => [
        "Overview" => "theory.md",
        "Manuscript example" => "generated/time_discretization.md",
        "Reconstruction weights" => "generated/reconstruction_weights.md",
        "Theory API" => "theory_api.md",
    ],
    "API reference" => "api_reference.md",
]

makedocs(;
    modules=[TimestepperTestCases, TimestepperTestCases.BuoyancyVarianceDissipationComputations, TimestepperTheory],
    authors="Simone Silvestri <silvestri.simone0@gmail.com> and contributors",
    sitename="TimestepperTestCases.jl",
    format=Documenter.HTML(;
        canonical="https://simone-silvestri.github.io/TimestepperTestCases.jl",
        edit_link="main",
        prettyurls=on_ci,
        assets=String[],
    ),
    pages=pages,
    checkdocs=:exports,
    warnonly=[:missing_docs, :cross_references],
)
