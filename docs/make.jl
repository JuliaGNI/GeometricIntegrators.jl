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
                "VPRK"        => "integrators/vprk.md",
                "SPARK"       => "integrators/spark.md",
               #  "CGVI"        => "integrators/cgvi.md",
               #  "DGVI"        => "integrators/dgvi.md",
               #  "HPG"         => "integrators/hpg.md",
                "Splitting"   => "integrators/splitting.md",
                "Stochastic"  => "integrators/stochastic.md"],
             "Modules"  => [
                "Basis Functions"   => "modules/basis_functions.md",
                "Discontinuities"   => "modules/discontinuities.md",
                "Equations"         => "modules/equations.md",
                "Integrators"       => "modules/integrators.md",
                "Interpolation"     => "modules/interpolation.md",
                "Linear Solvers"    => "modules/solvers_linear.md",
                "Nonlinear Solvers" => "modules/solvers_nonlinear.md",
                "Quadrature Rules"  => "modules/quadratures.md",
                "Simulations"       => "modules/simulations.md",
                "Solutions"         => "modules/solutions.md",
                "Tableaus"          => "modules/tableaus.md"],
             "Release Notes" => "releasenotes.md",
             "Bibliography" => "bibliography.md",
             ]
)

deploydocs(
    repo   = "github.com/JuliaGNI/GeometricIntegrators.jl",
)
