using Documenter, GeometricIntegrators

makedocs(
    sitename = "GeometricIntegrators.jl",
    format = :html,
    pages = ["Home" => "index.md",
             "Tutorial" => "tutorial.md",
             "Modules"  => [
                "Basis Functions"   => "modules/basis_functions.md",
                "Equations"         => "modules/equations.md",
                "Integrators"       => "modules/integrators.md",
                "Interpolation"     => "modules/interpolation.md",
                "Linear Solvers"    => "modules/solvers_linear.md",
                "Nonlinear Solvers" => "modules/solvers_nonlinear.md",
                "Solutions"         => "modules/solutions.md",
                "Tableaus"          => "modules/tableaus.md",
                "Example Problems"  => "modules/problems.md"]
             ]
)

deploydocs(
    repo   = "github.com/DDMGNI/GeometricIntegrators.jl.git",
    target = "build",
    deps   = nothing,
    make   = nothing
)
