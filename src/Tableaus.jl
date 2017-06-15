__precompile__()

module Tableaus

    using ..CommonFunctions
    using ..BasisFunctions
    using ..Quadratures
    using ..Utils


    export CoefficientsRK, CoefficientsARK, CoefficientsPRK, CoefficientsMRK

    include("tableaus/coefficients.jl")

    export get_symplectic_conjugate_coefficients, check_symplecticity, check_symmetry,
           check_order_conditions_B, check_order_conditions_C, check_order_conditions_D,
           compute_symplecticity_error

    include("tableaus/coefficients_symplectic.jl")

    export getCoefficientsGLRK,
           getCoefficientsGLRK1, getCoefficientsGLRK2, getCoefficientsGLRK3,
           getCoefficientsGLRK4, getCoefficientsGLRK5, getCoefficientsGLRK6

    include("tableaus/coefficients_glrk.jl")

    export getCoefficientsLobIII2,  getCoefficientsLobIII3,  getCoefficientsLobIII4,
           getCoefficientsLobIIIA2, getCoefficientsLobIIIA3, getCoefficientsLobIIIA4,
           getCoefficientsLobIIIB2, getCoefficientsLobIIIB3, getCoefficientsLobIIIB4,
           getCoefficientsLobIIIC2, getCoefficientsLobIIIC3, getCoefficientsLobIIIC4,
           getCoefficientsLobIIID2, getCoefficientsLobIIID3, getCoefficientsLobIIID4,
           getCoefficientsLobIIIE2, getCoefficientsLobIIIE3, getCoefficientsLobIIIE4,
           getCoefficientsLobIIIF2, getCoefficientsLobIIIF3, getCoefficientsLobIIIF4

    include("tableaus/coefficients_lob.jl")

    export getCoefficientsSRK3

    include("tableaus/coefficients_srk.jl")


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


    include("tableaus/tableaus_erk.jl")

    export getTableauExplicitEuler, getTableauExplicitMidpoint, getTableauHeun,
           getTableauKutta, getTableauERK4, getTableauERK438

    include("tableaus/tableaus_dirk.jl")

    export getTableauCrouzeix

    include("tableaus/tableaus_firk.jl")

    export getTableauImplicitEuler, getTableauImplicitMidpoint,
           getTableauGLRK1, getTableauGLRK2, getTableauGLRK3,
           getTableauGLRK, getTableauSRK3

    include("tableaus/tableaus_sirk.jl")

    include("tableaus/tableaus_eprk.jl")

    export getTableauSymplecticEulerA, getTableauSymplecticEulerB

    include("tableaus/tableaus_spark.jl")

    export getTableauLobIIIAIIIB2, getTableauLobIIIAIIIB3, getTableauLobIIIAIIIB4,
           getTableauLobIIIC2, getTableauLobIIIC3, getTableauLobIIIC4,
           getTableauLobIIID2, getTableauLobIIID3, getTableauLobIIID4,
           getTableauLobIIIE2, getTableauLobIIIE3, getTableauLobIIIE4,
           getTableauLobIIIF2, getTableauLobIIIF3, getTableauLobIIIF4,
           getTableauVPGLRK,
           getTableauVPSRK3


    include("tableaus/tableaus_vprk.jl")

    export getTableauSymplecticProjection,
           getTableauLobIIIAIIIB2pSymplectic, getTableauLobIIIAIIIB3pSymplectic,
           getTableauGLRKpSymplectic

    include("tableaus/tableaus_vpark.jl")

    export getTableauSymmetricProjection,
           getTableauLobIIIAIIIB2pSymmetric, getTableauLobIIIAIIIB3pSymmetric,
           getTableauGLRKpSymmetric

    include("tableaus/tableaus_vspark.jl")

end
