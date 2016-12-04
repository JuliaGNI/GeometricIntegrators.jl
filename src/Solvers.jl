__precompile__()

module Solvers

    using ..Utils

    include("utils/macro_utils.jl")

    export LinearSolver, LUSolver, LUSolverLAPACK,
           factorize!, solve!

    include("solvers/linear/linear_solvers.jl")
    include("solvers/linear/lu_solver.jl")
    include("solvers/linear/lu_solver_lapack.jl")

    export NonlinearSolver, AbstractNewtonSolver, NewtonSolver, QuasiNewtonSolver,
           NonlinearFunctionParameters,
           solve!, function_stages!

    include("solvers/nonlinear/nonlinear_solvers.jl")
    include("solvers/nonlinear/jacobian.jl")
    include("solvers/nonlinear/abstract_fixed_point_solver.jl")
    include("solvers/nonlinear/abstract_newton_solver.jl")
    include("solvers/nonlinear/fixed_point_solver.jl")
    include("solvers/nonlinear/newton_solver.jl")
    include("solvers/nonlinear/quasi_newton_solver.jl")

end
