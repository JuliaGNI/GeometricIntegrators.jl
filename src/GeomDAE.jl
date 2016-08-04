__precompile__()

module GeomDAE

include("utils/matrix_utils.jl")
include("integrators/tableaus.jl")

export Tableau, TableauRK, TableauERK, TableauIRK, TableauNLIRK,
       TableauPRK, TableauSPARK, TableauGLM

end
