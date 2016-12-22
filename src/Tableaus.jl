__precompile__()

module Tableaus

    using ..Utils
    
    include("utils/macro_utils.jl")

    export CoefficientsRK, CoefficientsARK, show_coefficients

    include("tableaus/coefficients.jl")

    export AbstractTableau, showTableau

    include("tableaus/tableaus.jl")

    export AbstractTableauRK, AbstractTableauIRK, AbstractTableauPRK,
           writeTableauToFile

    include("tableaus/abstract_tableau_rk.jl")

    export TableauERK, readTableauERKFromFile

    include("tableaus/tableau_erk.jl")

    export TableauDIRK, TableauFIRK, TableauSIRK

    include("tableaus/tableau_dirk.jl")
    include("tableaus/tableau_firk.jl")
    include("tableaus/tableau_sirk.jl")

    export TableauEPRK, TableauIPRK, TableauVPRK

    include("tableaus/tableau_eprk.jl")
    include("tableaus/tableau_iprk.jl")
    include("tableaus/tableau_vprk.jl")

    export TableauIPARK, TableauSARK, TableauSPARK

    include("tableaus/tableau_ipark.jl")
    include("tableaus/tableau_sark.jl")
    include("tableaus/tableau_spark.jl")

    export TableauGLM

    include("tableaus/tableau_glm.jl")


    include("tableaus/tableaus_erk.jl")

    export getTableauExplicitEuler, getTableauExplicitMidpoint, getTableauHeun,
           getTableauKutta, getTableauERK4, getTableauERK438

    include("tableaus/tableaus_dirk.jl")

    export getTableauCrouzeix

    include("tableaus/tableaus_firk.jl")

    export getTableauImplicitEuler, getTableauImplicitMidpoint,
           getTableauGLRK1, getTableauGLRK2, getTableauGLRK3

    include("tableaus/tableaus_sirk.jl")

    include("tableaus/tableaus_eprk.jl")

    export getTableauSymplecticEulerA, getTableauSymplecticEulerB

    include("tableaus/tableaus_vprk.jl")

    export getTableauLobIIIAB2

    include("tableaus/tableaus_ipark.jl")

    export getTableauSymplecticProjection

    include("tableaus/tableaus_spark.jl")

end
