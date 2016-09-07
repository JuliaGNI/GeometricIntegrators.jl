__precompile__()

module GeomDAE

include("utils/macro_utils.jl")
include("utils/matrix_utils.jl")

include("solvers/linear/linear_solvers.jl")

export LinearSolver, factorize!, solve!

include("solvers/linear/lu_solver_lapack.jl")

export LUSolver, factorize!, solve!

include("solvers/nonlinear/nonlinear_solvers.jl")

export NonlinearSolver, solve!

include("solvers/nonlinear/jacobian.jl")

include("solvers/nonlinear/abstract_newton_solver.jl")

export AbstractNewtonSolver

include("solvers/nonlinear/newton_solver.jl")

export NewtonSolver, solve!

include("solvers/nonlinear/quasi_newton_solver.jl")

export QuasiNewtonSolver, solve!

include("equations/equations.jl")

export Equation, ODE, PODE, DAE, PDAE

include("integrators/tableaus.jl")

export Tableau, TableauRK, TableauERK, TableauDIRK, TableauFIRK, TableauSIRK,
       TableauPRK, TableauSARK, TableauSPARK, TableauGLM,
       showTableau, writeTableauToFile, readTableauERKFromFile

include("integrators/tableaus_erk.jl")

export getTableauExplicitEuler, getTableauExplicitMidpoint, getTableauHeun,
       getTableauKutta, getTableauERK4, getTableauERK438

include("integrators/tableaus_dirk.jl")

export getTableauCrouzeix

include("integrators/tableaus_firk.jl")

export getTableauImplicitEuler, getTableauImplicitMidpoint,
       getTableauGLRK1, getTableauGLRK2, getTableauGLRK3

include("integrators/tableaus_prk.jl")

export getTableauSymplecticEulerA, getTableauSymplecticEulerB

include("integrators/solutions.jl")

export Solution, SolutionODE, SolutionPODE, SolutionDAE, SolutionPDAE,
       reset

include("integrators/integrators.jl")

export Integrator, IntegratorERK, IntegratorDIRK, IntegratorFIRK,
       IntegratorPRK, IntegratorSARK, IntegratorSPARK,
       solve, solve!

include("utils/hdf5_utils.jl")

export createHDF5, writeSolutionToHDF5

end
