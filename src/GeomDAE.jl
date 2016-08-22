__precompile__()

module GeomDAE

include("utils/matrix_utils.jl")

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

include("integrators/solutions.jl")

export Solution, SolutionODE, SolutionPODE, SolutionDAE, SolutionPDAE,
       reset

include("integrators/integrators.jl")

export Integrator, IntegratorERK, IntegratorDIRK, IntegratorFIRK,
       IntegratorPRK, IntegratorSARK, IntegratorSPARK,
       solve, solve!

end
