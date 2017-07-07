__precompile__()

module Solvers

    using ..Config
    using ..Utils

    export LinearSolver, LUSolver, LUSolverLAPACK,
           factorize!, solve!

    include("solvers/linear/linear_solvers.jl")
    include("solvers/linear/lu_solver.jl")
    include("solvers/linear/lu_solver_lapack.jl")

    export NonlinearSolver, AbstractNewtonSolver, NewtonSolver, QuasiNewtonSolver,
           NonlinearFunctionParameters,
           residual_initial!, residual_absolute!, residual_relative!,
           printSolverStatus, solverConverged, solverStatusOK,
           solve!, function_stages!

    include("solvers/nonlinear/nonlinear_solvers.jl")
    include("solvers/nonlinear/jacobian.jl")
    include("solvers/nonlinear/abstract_fixed_point_solver.jl")
    include("solvers/nonlinear/abstract_newton_solver.jl")
    include("solvers/nonlinear/fixed_point_solver.jl")
    include("solvers/nonlinear/newton_solver.jl")
    include("solvers/nonlinear/quasi_newton_solver.jl")


    function __init__()
        default_params = (
            (:ls_solver, :lapack),
            (:nls_atol,  2eps()),
            (:nls_rtol,  2eps()),
            (:nls_stol,  2eps()),
            (:nls_nmax,  10000),
            (:nls_nmin,  0),
            (:nls_nwarn, 100),
            (:nls_solver, QuasiNewtonSolver),
            (:jacobian_autodiff, true),
            (:jacobian_fd_ϵ, 8sqrt(eps())),
            (:quasi_newton_refactorize, 5),
            (:linesearch_nmax, 50),
            (:linesearch_armijo_λ₀, 1.0),
            (:linesearch_armijo_σ₀, 0.1),
            (:linesearch_armijo_σ₁, 0.5),
            (:linesearch_armijo_ϵ,  0.5),
        )

        for param in default_params
            add_config(param...)
        end
    end

end
