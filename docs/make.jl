using Documenter
using DocumenterCitations
using GeometricIntegrators
using RungeKutta
using Weave


ENV["GKSwstype"] = "100"

cp(normpath(@__FILE__, "../../AUTHORS.md"), normpath(@__FILE__, "../src/authors.md"); force=true)

bib = CitationBibliography(joinpath(@__DIR__, "src", "GeometricIntegrators.bib"))

DocMeta.setdocmeta!(GeometricIntegrators, :DocTestSetup, :(using GeometricIntegrators); recursive=true)

weave(joinpath(@__DIR__, "src", "methods.jmd"), out_path = joinpath(@__DIR__, "src"), doctype = "github")

makedocs(bib,
    sitename = "GeometricIntegrators.jl",
    format = Documenter.HTML(
               prettyurls = get(ENV, "CI", nothing) == "true",
               mathengine = MathJax3(Dict(
                 :tex => Dict(
                   :macros => Dict(
                     :bigtimes => "\\mathop{\\vcenter{\\huge\\times}}",
                     :ext => "\\mathsf{d}",
                     :id => "\\mathop{id}",
                     :rank => "\\mathop{rank}",
                     :spn => "\\mathop{span}",
                     :identity => "\\mathbb{I}",
                     :abs => ["\\left \\vert #1 \\right \\vert", 1],
                     :mf => ["\\mathcal{#1}", 1],
                     :tb => ["\\mathsf{T}_{#1} #2", 2, ""],
                     :cb => ["\\mathsf{T}_{#1}^{*} #2", 2, ""],
                   ),
                   :inlineMath => [["\$","\$"], ["\\(","\\)"]],
                   :tags => "ams",
                 ),
                 :options => Dict(
                   :ignoreHtmlClass => "tex2jax_ignore",
                   :processHtmlClass => "tex2jax_process",
                 ),
               ), true)
             ),
    pages = ["Home" => "index.md",
             "Tutorial" => "tutorial.md",
             "Problems" => "problems.md",
             "Methods"    => "methods.md",
             "Integrators" => [
                "Runge-Kutta" => "integrators/rk.md",
                "Splitting"   => "integrators/splitting.md",
                "Variational" => "integrators/variational.md",
                "VPRK"        => "integrators/vprk.md",
                "SPARK"       => "integrators/spark.md",
                "DVI"         => "integrators/dvi.md",
                # "CGVI"        => "integrators/cgvi.md",
               #  "DGVI"        => "integrators/dgvi.md",
               #  "HPI"         => "integrators/hpi.md",
               #  "HPG"         => "integrators/hpg.md",
               ],
             "Modules" => [
               # "Discontinuities"     => "modules/discontinuities.md",
                 "Methods"             => "modules/methods.md",
                 "Integrators"         => "modules/integrators.md",
                 "Problems"            => "modules/problems.md",
                 "Equations"           => "modules/equations.md",
               # "Simulations"         => "modules/simulations.md",
                 "Solutions"           => "modules/solutions.md",
                 "Runge-Kutta Tableaus" => "modules/rungekutta.md",
                 "Partitioned Runge-Kutta Tableaus" => "modules/rungekutta_partitioned.md",
                 "SPARK Methods"       => "modules/spark.md",
                 ],
             "Developer Docs" =>[
                "Integrators"         => "developer/integrators.md",
                "Projections"         => "developer/projections.md",
                "Code Integration"    => "developer/code_integration.md",
                "Custom Integrators"  => "developer/custom_integrators.md",
                "Adaptive Time Stepping" => "developer/adaptive_time_stepping.md",
               ],
             "Release Notes" => "releasenotes.md",
             "Bibliography" => "bibliography.md",
             "Authors" => "authors.md",
             "License" => "LICENSE.md",
             ],
    modules = [GeometricIntegrators,
               GeometricBase,
               GeometricEquations,
               GeometricSolutions,
               RungeKutta]
)

deploydocs(
    repo   = "github.com/JuliaGNI/GeometricIntegrators.jl",
    devurl = "latest",
    devbranch = "main",
    push_preview = true,
)
