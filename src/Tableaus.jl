__precompile__()

module Tableaus

    using ..BasisFunctions
    using ..Utils

    export CoefficientsRK, CoefficientsARK, CoefficientsPRK, CoefficientsMRK

    include("tableaus/coefficients.jl")

    export AbstractTableau

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

    export TableauARK, TableauSARK

    include("tableaus/tableau_ark.jl")
    include("tableaus/tableau_sark.jl")

    export TableauPARK, TableauSPARK, TableauVPARK, TableauVSPARK

    include("tableaus/tableau_park.jl")
    include("tableaus/tableau_spark.jl")
    include("tableaus/tableau_vpark.jl")
    include("tableaus/tableau_vspark.jl")

    export TableauGLM

    include("tableaus/tableau_glm.jl")


    export getTableauGLRK

    include("tableaus/tableaus_glrk.jl")


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

    include("tableaus/tableaus_spark.jl")

    export getTableauLobIIIAB2

    include("tableaus/tableaus_vprk.jl")

    export getTableauSymplecticProjection,
           getTableauLobIIIAB2p

    include("tableaus/tableaus_vpark.jl")

    export getTableauSymmetricSymplecticProjection,
           getTableauLobIIIAB2sp

    include("tableaus/tableaus_vspark.jl")

end
