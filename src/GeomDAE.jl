__precompile__()

module GeomDAE

include("utils/matrix_utils.jl")

include("equations/equations.jl")

export Equation, ODE, PODE, DAE, PDAE

include("integrators/tableaus.jl")

export Tableau, TableauRK, TableauERK, TableauIRK, TableauNLIRK, TableauPRK,
       TableauSPARK, TableauGLM,
       showTableau, writeTableauToFile, readTableauERKFromFile

include("integrators/tableaus_erk.jl")
include("integrators/tableaus_irk.jl")
include("integrators/tableaus_nlirk.jl")

export getTableauExplicitMidpoint, getTableauHeun
export getTableauCrouzeix
export getTableauImplicitMidpoint, getTableauGLRK1, getTableauGLRK2

include("integrators/solutions.jl")

export Solution, SolutionODE, SolutionPODE, SolutionDAE, SolutionPDAE,
       reset

include("integrators/integrators.jl")


end
