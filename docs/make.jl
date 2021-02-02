using Documenter
using DocumenterCitations
using Weave
using GeometricIntegrators

bib = CitationBibliography("GeometricIntegrators.bib")

makedocs(bib,
    sitename = "GeometricIntegrators.jl",
    format = Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true"),
    pages = ["Home" => "index.md",
             "Tutorial" => "tutorial/tutorial.md",
             "Integrators" => [
                "Overview"    => "integrators.md",
                "Runge-Kutta" => "integrators/rk.md",
                "Splitting"   => "integrators/splitting.md",
                "VPRK"        => "integrators/vprk.md",
                "SPARK"       => "integrators/spark.md",
               #  "CGVI"        => "integrators/cgvi.md",
               #  "DGVI"        => "integrators/dgvi.md",
               #  "HPG"         => "integrators/hpg.md",
               ],
             "Modules" => [
                "Discontinuities"   => "modules/discontinuities.md",
                "Equations"         => "modules/equations.md",
                "Integrators"       => "modules/integrators.md",
                "Interpolation"     => "modules/interpolation.md",
                "Simulations"       => "modules/simulations.md",
                "Solutions"         => "modules/solutions.md",
               ],
             "Tableaus" => [
                "Runge-Kutta Methods"  => "tableaus/rungekutta.md",
                "Partitioned Runge-Kutta Methods"  => "tableaus/prk.md",
                "Splitting Methods"  => "tableaus/splitting.md",
                "VPRK Methods"  => "tableaus/vprk.md",
                "SPARK Methods"  => "tableaus/spark.md",
               ],
             "Release Notes" => "releasenotes.md",
             "Bibliography" => "bibliography.md",
             ]
)

deploydocs(
    repo   = "github.com/JuliaGNI/GeometricIntegrators.jl",
)
