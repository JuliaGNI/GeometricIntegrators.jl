using Documenter, GeometricIntegrators

makedocs(
    sitename = "GeometricIntegrators.jl",
    format = :html,
    pages = ["Home" => "index.md",
             "Tutorial" => "tutorial.md",
             "Modules"  => [
                "Equations"         => "modules/equations.md",
                "Integrators"       => "modules/integrators.md",
                "Linear Solvers"    => "modules/solvers_linear.md",
                "Nonlinear Solvers" => "modules/solvers_nonlinear.md",
                "Tableaus"          => "modules/tableaus.md"]
             ]
)
