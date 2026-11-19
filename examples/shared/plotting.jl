using CairoMakie
using LaTeXStrings

manuscript_font = joinpath(pkgdir(Makie.MathTeXEngine), "assets", "fonts", "ComputerModern", "cmr10.ttf")
set_theme!(Theme(fontsize=16, fonts=(; regular=manuscript_font, bold=manuscript_font), Legend=(; framevisible=true)))

function save_example(figure, name; pt_per_unit=0.75)
    directory = joinpath(@__DIR__, "..", "output")
    mkpath(directory)
    for extension in ("png", "pdf")
        path = joinpath(directory, name * "." * extension)
        extension == "pdf" ? save(path, figure; pt_per_unit) : save(path, figure)
    end
    return figure
end

publication_label(label) = replace(label, "WRK3-SE-SM05" => "WRK3-SE-F1")

function signed_power_label(value)
    iszero(value) && return L"0"
    sign = value < 0 ? "-" : ""
    exponent = round(Int, log10(abs(value)))
    return LaTeXString(sign * "10^{" * string(exponent) * "}")
end
