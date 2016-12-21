using Documenter, GeometricIntegrators

makedocs(
    sitename = "GeometricIntegrators.jl",
    format = :html,
    pages = ["Home" => "index.md",
             "Tutorial" => "tutorial.md",
             "Modules"  => [
                "Equations"         => "modules/equations.md",
                "Solutions"         => "modules/solutions.md",
                "Integrators"       => "modules/integrators.md",
                "Interpolation"     => "modules/interpolation.md",
                "Linear Solvers"    => "modules/solvers_linear.md",
                "Nonlinear Solvers" => "modules/solvers_nonlinear.md",
                "Tableaus"          => "modules/tableaus.md"]
             ]
)

deploydocs(
    repo   = "github.com/DDMGNI/GeometricIntegrators.jl.git",
    target = "build",
    deps   = nothing,
    make   = nothing
)
