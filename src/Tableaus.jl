module Tableaus

    using LinearAlgebra: mul!

    using ..Config
    using ..CommonFunctions
    using ..BasisFunctions
    using ..Quadratures
    using ..Integrators
    using ..Integrators.Stochastic
    using ..Integrators.SPARK
    using ..Integrators.VPRK
    using ..Utils

    using GeometricIntegrators.Integrators.SPARK: get_ã_vspark_primary,
                                                  get_α_vspark_primary,
                                                  compute_ã_vspark_primary,
                                                  compute_α_vspark_primary

    export getCoefficientsGLRK,
           getCoefficientsGLRK1, getCoefficientsGLRK2, getCoefficientsGLRK3,
           getCoefficientsGLRK4, getCoefficientsGLRK5, getCoefficientsGLRK6

    include("tableaus/coefficients_glrk.jl")


    export getCoefficientsLobIII2,  getCoefficientsLobIII3,  getCoefficientsLobIII4,  getCoefficientsLobIII5,
           getCoefficientsLobIIIA2, getCoefficientsLobIIIA3, getCoefficientsLobIIIA4, getCoefficientsLobIIIA5,
           getCoefficientsLobIIIB2, getCoefficientsLobIIIB3, getCoefficientsLobIIIB4, getCoefficientsLobIIIB5,
           getCoefficientsLobIIIC2, getCoefficientsLobIIIC3, getCoefficientsLobIIIC4, getCoefficientsLobIIIC5,
           getCoefficientsLobIIID2, getCoefficientsLobIIID3, getCoefficientsLobIIID4, getCoefficientsLobIIID5,
           getCoefficientsLobIIIE2, getCoefficientsLobIIIE3, getCoefficientsLobIIIE4, getCoefficientsLobIIIE5,
           getCoefficientsLobIIIF2, getCoefficientsLobIIIF3, getCoefficientsLobIIIF4,
           getCoefficientsLobIIIG2, getCoefficientsLobIIIG3, getCoefficientsLobIIIG4

    export getCoefficientsLobIII,  getCoefficientsLobIIIA, getCoefficientsLobIIIB,
           getCoefficientsLobIIIC, getCoefficientsLobIIID, getCoefficientsLobIIIE,
           getCoefficientsLobIIIF, getCoefficientsLobIIIG

    include("tableaus/coefficients_lob.jl")

    export getCoefficientsRadIIA2, getCoefficientsRadIIA3

    include("tableaus/coefficients_rad.jl")

    export getCoefficientsPGLRK #, getTableauPGLRK

    include("tableaus/coefficients_pglrk.jl")

    export getCoefficientsSRK3

    include("tableaus/coefficients_srk.jl")


    include("tableaus/tableaus_erk.jl")

    export getTableauExplicitEuler, getTableauExplicitMidpoint, getTableauHeun,
           getTableauKutta, getTableauRunge, getTableauERK4, getTableauERK438

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

    include("tableaus/tableaus_iprk.jl")

    export getTableauIPGLRK           
    
    include("tableaus/tableaus_spark.jl")

    export TableauSPARKGLRK,
           TableauSPARKLobIIIAIIIB,
           TableauSPARKLobIIIBIIIA,
           TableauSPARKGLRKLobIIIAIIIB,
           TableauSPARKGLRKLobIIIBIIIA,
           TableauSPARKLobatto,
           TableauSPARKLobABC,
           TableauSPARKLobABD,
           TableauSPARKVPRK,
           TableauSPARKGLVPRK

    include("tableaus/tableaus_vprk.jl")

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

    export TableauSymplecticProjection,
           TableauLobIIIAIIIBpSymplectic,
           TableauLobIIIBIIIApSymplectic,
           TableauGLRKpSymplectic

    include("tableaus/tableaus_vpark.jl")

    export TableauVSPARKLobIIIAIIIBProjection,
           TableauVSPARKLobIIIBIIIAProjection,
           TableauVSPARKModifiedLobIIIAIIIBProjection,
           TableauVSPARKModifiedLobIIIBIIIAProjection,
           TableauVSPARKInternalProjection,
           TableauVSPARKModifiedInternalProjection,
           TableauVSPARKMidpointProjection,
           TableauVSPARKModifiedMidpointProjection,
           TableauVSPARKSymmetricProjection,
           TableauVSPARKSymplecticProjection,
           TableauVSPARKGLRKpLobIIIAIIIB,
           TableauVSPARKGLRKpLobIIIBIIIA,
           TableauVSPARKGLRKpModifiedLobIIIAIIIB,
           TableauVSPARKGLRKpModifiedLobIIIBIIIA,
           TableauVSPARKGLRKpInternal,
           TableauVSPARKGLRKpModifiedInternal,
           TableauVSPARKGLRKpMidpoint,
           TableauVSPARKGLRKpModifiedMidpoint,
           TableauVSPARKGLRKpSymmetric,
           TableauVSPARKGLRKpSymplectic,
           TableauVSPARKLobIIIAIIIBpLobIIIAIIIB,
           TableauVSPARKLobIIIBIIIApLobIIIAIIIB,
           TableauVSPARKLobIIIAIIIBpLobIIIBIIIA,
           TableauVSPARKLobIIIBIIIApLobIIIBIIIA,
           TableauVSPARKLobIIIAIIIBpModifiedLobIIIAIIIB,
           TableauVSPARKLobIIIAIIIBpModifiedLobIIIBIIIA,
           TableauVSPARKLobIIIBIIIApModifiedLobIIIAIIIB,
           TableauVSPARKLobIIIBIIIApModifiedLobIIIBIIIA,
           TableauVSPARKLobIIIAIIIBpMidpoint,
           TableauVSPARKLobIIIBIIIApMidpoint,
           TableauVSPARKLobIIIAIIIBpModifiedMidpoint,
           TableauVSPARKLobIIIBIIIApModifiedMidpoint,
           TableauVSPARKLobIIIAIIIBpSymmetric,
           TableauVSPARKLobIIIBIIIApSymmetric,
           TableauVSPARKLobABCCD,
           TableauVSPARKLobABCCE,
           TableauVSPARKLobABDE,
           TableauVSPARKLobABED,
           TableauVSPARKLobABD,
           TableauVSPARKLobABE,
           TableauVSPARKLobDE,
           TableauVSPARKLobED

    include("tableaus/tableaus_vspark_primary.jl")

    export getTableauVSPARKLobIIIAB,
           getTableauVSPARKLobIIIC,
           getTableauVSPARKLobIIID,
           getTableauVSPARKLobIIIE,
           getTableauVSPARKGLRKLobIIIAB,
           getTableauVSPARKGLRKLobIIIC,
           getTableauVSPARKGLRKLobIIID,
           getTableauVSPARKGLRKLobIIIE

    include("tableaus/tableaus_vspark_secondary.jl")

    export getTableauHPARK,
           TableauHPARKGLRK,
           TableauHPARKLobIIIAIIIB,
           TableauHPARKLobIIIBIIIA

    include("tableaus/tableaus_hpark.jl")

    export TableauHSPARKSymmetricProjection,
           TableauHSPARKGLRKpSymmetric,
           TableauHSPARKLobIIIAIIIBpSymmetric,
           TableauHSPARKLobIIIBIIIApSymmetric

    include("tableaus/tableaus_hspark_primary.jl")

    export TableauHSPARKLobIIIAB,
           TableauHSPARKLobIIIBA,
           TableauHSPARKLobIIIC,
           TableauHSPARKLobIIID,
           TableauHSPARKLobIIIE,
           TableauHSPARKGLRKLobIIIAB,
           TableauHSPARKGLRKLobIIIBA,
           TableauHSPARKGLRKLobIIIC,
           TableauHSPARKGLRKLobIIID,
           TableauHSPARKGLRKLobIIIE

    include("tableaus/tableaus_hspark_secondary.jl")

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
