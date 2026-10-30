using Literate
using Markdown

output_directory = joinpath(@__DIR__, "..", "output")
source = joinpath(@__DIR__, "time_discretization.jl")
markdown_path = Literate.markdown(source, output_directory; execute=false, flavor=Literate.CommonMarkFlavor(), credit=false)
content = sprint(show, MIME("text/html"), Markdown.parse(read(markdown_path, String)))

html = """
<!doctype html>
<html lang="en">
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Time discretization — TimestepperTheory</title>
<style>
body { max-width: 980px; margin: 48px auto; padding: 0 28px; color: #20252a; background: #fff; font: 17px/1.65 Georgia, serif; }
h1, h2 { line-height: 1.25; font-weight: normal; } h1 { margin-bottom: 28px; } h2 { margin-top: 48px; border-top: 1px solid #ddd; padding-top: 24px; }
pre { padding: 18px 22px; background: #f6f7f8; border: 1px solid #e2e5e8; border-radius: 4px; overflow: auto; }
code { font: 13px/1.65 Menlo, monospace; } img { display: block; max-width: 100%; height: auto; margin: 28px auto; }
a { color: #136a99; } @media (max-width: 600px) { body { margin-top: 24px; padding: 0 16px; } }
</style>
<main>$content</main>
</html>
"""

write(joinpath(output_directory, "time_discretization.html"), html)
@info "Rendered time discretization" path=joinpath(output_directory, "time_discretization.html")
