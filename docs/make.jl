using Documenter
using GeometricIntegrators

makedocs(
    sitename = "GeometricIntegrators.jl",
    format = Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true"),
    pages = ["Home" => "index.md",
             "Tutorial" => "tutorial.md",
             "Integrators" => [
                "Overview"    => "integrators.md",
                "Splitting"   => "integrators/splitting.md",
                "Runge-Kutta" => "integrators/rk.md",
                "VPRK"        => "integrators/vprk.md",
                "SPARK"       => "integrators/spark.md",
                "CGVI"        => "integrators/cgvi.md",
                "DGVI"        => "integrators/dgvi.md",
                "HPG"         => "integrators/hpg.md",
                "Stochastic"  => "integrators/stochastic.md"],
             "Modules"  => [
                "Basis Functions"   => "modules/basis_functions.md",
                "Discontinuities"   => "modules/discontinuities.md",
                "Equations"         => "modules/equations.md",
                "Integrators"       => "modules/integrators.md",
                "Stochastic Integrators"  => "modules/integrators_stochastic.md",
                "SPARK Integrators"       => "modules/integrators_spark.md",
                "VPRK Integrators"        => "modules/integrators_vprk.md",
                "Interpolation"     => "modules/interpolation.md",
                "Linear Solvers"    => "modules/solvers_linear.md",
                "Nonlinear Solvers" => "modules/solvers_nonlinear.md",
                "Quadrature Rules"  => "modules/quadratures.md",
                "Simulations"       => "modules/simulations.md",
                "Solutions"         => "modules/solutions.md",
                "Tableaus"          => "modules/tableaus.md"],
             ]
)

deploydocs(
    repo   = "github.com/DDMGNI/GeometricIntegrators.jl.git",
)
