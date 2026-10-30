function save_results(directory, name, figure, rows)
    mkpath(directory)
    save(joinpath(directory, name * ".pdf"), figure)
    save(joinpath(directory, name * ".png"), figure; px_per_unit=220/72)
    open(joinpath(directory, name * ".csv"), "w") do stream
        println(stream, join(string.(keys(first(rows))), ','))
        for row in rows
            println(stream, join(values(row), ','))
        end
    end
    @info "Saved $name"
    return nothing
end
