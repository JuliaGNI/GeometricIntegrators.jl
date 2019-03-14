using Documenter, GeometricIntegrators

makedocs(
    sitename = "GeometricIntegrators.jl",
    format = Documenter.HTML(),
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
                "HPG"         => "integrators/hpg.md"],
             "Modules"  => [
                "Basis Functions"   => "modules/basis_functions.md",
                "Discontinuities"   => "modules/discontinuities.md",
                "Equations"         => "modules/equations.md",
                "Integrators"       => "modules/integrators.md",
                "Interpolation"     => "modules/interpolation.md",
                "Linear Solvers"    => "modules/solvers_linear.md",
                "Nonlinear Solvers" => "modules/solvers_nonlinear.md",
                "Quadrature Rules"  => "modules/quadratures.md",
                "Solutions"         => "modules/solutions.md",
                "Tableaus"          => "modules/tableaus.md"],
             ]
)

deploydocs(
    repo   = "github.com/DDMGNI/GeometricIntegrators.jl.git",
)
