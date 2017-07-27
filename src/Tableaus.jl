__precompile__()

module Tableaus

    using ..Config
    using ..CommonFunctions
    using ..BasisFunctions
    using ..Quadratures
    using ..Utils


    export CoefficientsRK, CoefficientsARK, CoefficientsPRK, CoefficientsMRK

    include("tableaus/coefficients.jl")

    export get_symplectic_conjugate_coefficients, symplecticize,
           check_symplecticity, check_symmetry, compute_symplecticity_error,
           check_order_conditions_B, check_order_conditions_C, check_order_conditions_D


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
           getCoefficientsLobIIIF2, getCoefficientsLobIIIF3, getCoefficientsLobIIIF4,
           getCoefficientsLobIIIG2, getCoefficientsLobIIIG3, getCoefficientsLobIIIG4

    include("tableaus/coefficients_lob.jl")

    export getCoefficientsRadIIA2, getCoefficientsRadIIA3

    include("tableaus/coefficients_rad.jl")

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

    export AbstractTableauSplitting, TableauSplittingGS, TableauSplittingNS, TableauSplittingSS

    include("tableaus/tableau_splitting.jl")


    include("tableaus/tableaus_erk.jl")

    export getTableauExplicitEuler, getTableauExplicitMidpoint, getTableauHeun,
           getTableauKutta, getTableauERK4, getTableauERK438

    include("tableaus/tableaus_dirk.jl")

    export getTableauCrouzeix

    include("tableaus/tableaus_firk.jl")

    export getTableauImplicitEuler, getTableauImplicitMidpoint,
           getTableauGLRK,
           getTableauLobIIIA2, getTableauLobIIIA3, getTableauLobIIIA4,
           getTableauLobIIIB2, getTableauLobIIIB3, getTableauLobIIIB4,
           getTableauLobIIIC2, getTableauLobIIIC3, getTableauLobIIIC4,
           getTableauLobIIID2, getTableauLobIIID3, getTableauLobIIID4,
           getTableauLobIIIE2, getTableauLobIIIE3, getTableauLobIIIE4,
           getTableauLobIIIF2, getTableauLobIIIF3, getTableauLobIIIF4,
           getTableauLobIIIG2, getTableauLobIIIG3, getTableauLobIIIG4,
           getTableauRadIIA2,  getTableauRadIIA3,
           getTableauSRK3

    include("tableaus/tableaus_sirk.jl")

    include("tableaus/tableaus_eprk.jl")

    export getTableauSymplecticEulerA, getTableauSymplecticEulerB

    include("tableaus/tableaus_spark.jl")

    export getTableauVPGLRK,
           getTableauVPLobIIIAIIIB2, getTableauVPLobIIIAIIIB3, getTableauVPLobIIIAIIIB4,
           getTableauVPLobIIIBIIIA2, getTableauVPLobIIIBIIIA3, getTableauVPLobIIIBIIIA4,
           getTableauVPLobIIIC2, getTableauVPLobIIIC3, getTableauVPLobIIIC4,
           getTableauVPLobIIID2, getTableauVPLobIIID3, getTableauVPLobIIID4,
           getTableauVPLobIIIE2, getTableauVPLobIIIE3, getTableauVPLobIIIE4,
           getTableauVPLobIIIF2, getTableauVPLobIIIF3, getTableauVPLobIIIF4,
           getTableauVPLobIIIG2, getTableauVPLobIIIG3, getTableauVPLobIIIG4,
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

    export getTableauLieA, getTableauLieB, getTableauStrang,
           getTableauMcLachlan2, getTableauMcLachlan4,
           getTableauTripleJump, getTableauSuzukiFractal

    include("tableaus/tableaus_splitting.jl")


    function __init__()
        add_config(:tab_compensated_summation, true)
    end

end
