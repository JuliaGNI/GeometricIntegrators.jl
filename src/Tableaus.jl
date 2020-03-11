module Tableaus

    using ..Config
    using ..CommonFunctions
    using ..BasisFunctions
    using ..Quadratures
    using ..Integrators
    using ..Integrators.SPARK
    using ..Integrators.Stochastic
    using ..Utils


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

    export getCoefficientsLobIII,  getCoefficientsLobIIIA, getCoefficientsLobIIIB,
           getCoefficientsLobIIIC, getCoefficientsLobIIID, getCoefficientsLobIIIE

    include("tableaus/coefficients_lob.jl")

    export getCoefficientsRadIIA2, getCoefficientsRadIIA3

    include("tableaus/coefficients_rad.jl")

    export getCoefficientsPGLRK, getTableauPGLRK

    include("tableaus/coefficients_pglrk.jl")

    export getCoefficientsSRK3

    include("tableaus/coefficients_srk.jl")


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
    include("tableaus/tableaus_siprk.jl")
    include("tableaus/tableaus_sisprk.jl")

    export  getTableauStochasticGLRK, getTableauStochasticDIRK
    export  getTableauStochasticStormerVerlet, getTableauStochasticSymplecticEuler
    export  getTableauStochasticLobIIIABD2, getTableauModifiedStochasticStormerVerlet

    include("tableaus/tableaus_serk.jl")

    export getTableauPlaten, getTableauBurrageR2, getTableauBurrageCL
    export getTableauBurrageE1, getTableauBurrageG5, getTableauStochasticHeun
    export getTableauStochasticEuler

    include("tableaus/tableaus_werk.jl")

    export getTableauRosslerRS1, getTableauRosslerRS2

    include("tableaus/tableaus_wirk.jl")

    export getTableauSRKw1, getTableauSRKw2

    include("tableaus/tableaus_eprk.jl")

    export getTableauSymplecticEulerA, getTableauSymplecticEulerB,
           getTableauLobattoIIIAIIIB2, getTableauLobattoIIIBIIIA2

    export getTableauSPARKGLRK

    include("tableaus/tableaus_spark.jl")

    export getTableauVPGLRK,
           getTableauVPLobIIIA2, getTableauVPLobIIIA3, getTableauVPLobIIIA4,
           getTableauVPLobIIIB2, getTableauVPLobIIIB3, getTableauVPLobIIIB4,
           getTableauVPLobIIIC2, getTableauVPLobIIIC3, getTableauVPLobIIIC4,
           getTableauVPLobIIID2, getTableauVPLobIIID3, getTableauVPLobIIID4,
           getTableauVPLobIIIE2, getTableauVPLobIIIE3, getTableauVPLobIIIE4,
           getTableauVPLobIIIF2, getTableauVPLobIIIF3, getTableauVPLobIIIF4,
           getTableauVPLobIIIG2, getTableauVPLobIIIG3, getTableauVPLobIIIG4,
           getTableauVPSRK3,
           getTableauVPLobIIIAIIIA2, getTableauVPLobIIIAIIIA3, getTableauVPLobIIIAIIIA4,
           getTableauVPRadIIAIIA2, getTableauVPRadIIAIIA3

    include("tableaus/tableaus_vprk.jl")


    export getTableauSymplecticProjection,
           getTableauLobIIIAIIIB2pSymplectic,
           getTableauLobIIIAIIIB3pSymplectic,
           getTableauLobIIIAIIIB4pSymplectic,
           getTableauGLRKpSymplectic

    include("tableaus/tableaus_vpark.jl")

    export getTableauVSPARKMidpointProjection,
           getTableauVSPARKSymmetricProjection,
           getTableauVSPARKSymplecticProjection,
           getTableauVSPARKGLRKpMidpoint,
           getTableauVSPARKGLRKpSymmetric,
           getTableauVSPARKGLRKpSymplectic,
           getTableauVSPARKLobIIIAIIIB2pSymmetric,
           getTableauVSPARKLobIIIAIIIB3pSymmetric,
           getTableauVSPARKLobIIIAIIIB4pSymmetric

    include("tableaus/tableaus_vspark.jl")

    export getTableauVSPARKLobIIIAB,
           getTableauVSPARKLobIIIC,
           getTableauVSPARKLobIIID,
           getTableauVSPARKLobIIIE,
           getTableauVSPARKGLRKLobIIIAB,
           getTableauVSPARKGLRKLobIIIC,
           getTableauVSPARKGLRKLobIIID,
           getTableauVSPARKGLRKLobIIIE

    include("tableaus/tableaus_vspark_psecondary.jl")

    export getTableauHPARK, getTableauHPARKGLRK,
           getTableauHPARKLobIIIAIIIB2,
           getTableauHPARKLobIIIAIIIB3,
           getTableauHPARKLobIIIAIIIB4

    include("tableaus/tableaus_hpark.jl")

    export getTableauHSPARKSymmetricProjection,
           getTableauHSPARKGLRKpSymmetric,
           getTableauHSPARKLobIIIAIIIB2pSymmetric,
           getTableauHSPARKLobIIIAIIIB3pSymmetric,
           getTableauHSPARKLobIIIAIIIB4pSymmetric

    include("tableaus/tableaus_hspark.jl")

    export getTableauSLRKLobIIIAB,
           getTableauSLRKLobIIIC,
           getTableauSLRKLobIIID,
           getTableauSLRKLobIIIE

    include("tableaus/tableaus_slrk.jl")

    export getTableauLieA, getTableauLieB, getTableauStrang,
           getTableauMcLachlan2, getTableauMcLachlan4,
           getTableauTripleJump, getTableauSuzukiFractal

    include("tableaus/tableaus_splitting.jl")


    function __init__()
        add_config(:tab_compensated_summation, true)
    end

end
