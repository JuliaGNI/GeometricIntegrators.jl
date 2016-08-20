__precompile__()

module GeomDAE

include("utils/matrix_utils.jl")

include("equations/equations.jl")

export Equation, ODE, PODE, DAE, PDAE

include("integrators/tableaus.jl")

export Tableau, TableauRK, TableauERK, TableauIRK, TableauNLIRK, TableauPRK,
       TableauSARK, TableauSPARK, TableauGLM,
       showTableau, writeTableauToFile, readTableauERKFromFile

include("integrators/tableaus_erk.jl")

export getTableauExplicitEuler, getTableauExplicitMidpoint, getTableauHeun,
       getTableauKutta, getTableauERK4, getTableauERK438

include("integrators/tableaus_irk.jl")

export getTableauCrouzeix

include("integrators/tableaus_nlirk.jl")

export getTableauImplicitEuler, getTableauImplicitMidpoint,
       getTableauGLRK1, getTableauGLRK2, getTableauGLRK3

include("integrators/solutions.jl")

export Solution, SolutionODE, SolutionPODE, SolutionDAE, SolutionPDAE,
       reset

include("integrators/integrators.jl")

export Integrator, IntegratorERK, IntegratorIRK, IntegratorNLIRK,
       IntegratorPRK, IntegratorSARK, IntegratorSPARK

end
