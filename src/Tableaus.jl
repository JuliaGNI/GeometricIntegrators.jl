module Tableaus

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

    include("tableaus/coefficients_sa.jl")

    export CoefficientsGLRK

    include("tableaus/coefficients_glrk.jl")


    export CoefficientsLobattoIIIA, CoefficientsLobattoIIIB,
           CoefficientsLobattoIIIC, CoefficientsLobattoIIIC̄,
           CoefficientsLobattoIIID, CoefficientsLobattoIIIE,
           CoefficientsLobattoIIIF, CoefficientsLobattoIIIG

    include("tableaus/coefficients_lob.jl")

    export CoefficientsRadauIA, CoefficientsRadauIIA

    include("tableaus/coefficients_rad.jl")

    export CoefficientsSRK3

    include("tableaus/coefficients_srk.jl")


    include("tableaus/tableaus_erk.jl")

    export TableauExplicitEuler, TableauExplicitMidpoint,
           TableauHeun, TableauRunge, TableauKutta,
           TableauERK4, TableauERK416, TableauERK438

    include("tableaus/tableaus_dirk.jl")

    export TableauCrouzeix

    include("tableaus/tableaus_firk.jl")

    export TableauImplicitEuler, TableauImplicitMidpoint,
           TableauGLRK, TableauRadauIA,  TableauRadauIIA, TableauSRK3,
           TableauLobattoIIIA, TableauLobattoIIIB, TableauLobattoIIIC, TableauLobattoIIIC̄,
           TableauLobattoIIID, TableauLobattoIIIE, TableauLobattoIIIF, TableauLobattoIIIG
           

    include("tableaus/tableaus_sirk.jl")
    include("tableaus/tableaus_siprk.jl")
    include("tableaus/tableaus_sisprk.jl")

    export  TableauStochasticGLRK, TableauStochasticDIRK
    export  TableauStochasticStormerVerlet, TableauStochasticSymplecticEuler
    export  TableauStochasticLobattoIIIABD2, TableauModifiedStochasticStormerVerlet

    include("tableaus/tableaus_serk.jl")

    export TableauPlaten, TableauBurrageR2, TableauBurrageCL
    export TableauBurrageE1, TableauBurrageG5, TableauStochasticHeun
    export TableauStochasticEuler

    include("tableaus/tableaus_werk.jl")

    export TableauRosslerRS1, TableauRosslerRS2

    include("tableaus/tableaus_wirk.jl")

    export TableauSRKw1, TableauSRKw2

    include("tableaus/tableaus_eprk.jl")

    export TableauSymplecticEulerA, TableauSymplecticEulerB,
           TableauLobattoIIIAIIIB2, TableauLobattoIIIBIIIA2

    include("tableaus/tableaus_iprk.jl")

    export TableauIPGLRK           
    
    include("tableaus/tableaus_spark.jl")

    export TableauSPARKGLRK,
           TableauSPARKLobattoIIIAIIIB,
           TableauSPARKLobattoIIIBIIIA,
           TableauSPARKGLRKLobattoIIIAIIIB,
           TableauSPARKGLRKLobattoIIIBIIIA,
           TableauSPARKLobatto,
           TableauSPARKLobABC,
           TableauSPARKLobABD,
           TableauSPARKVPRK,
           TableauSPARKGLVPRK

    include("tableaus/tableaus_vprk.jl")

    export TableauVPGLRK,
           TableauVPLobattoIIIA,
           TableauVPLobattoIIIB,
           TableauVPLobattoIIIC,
           TableauVPLobattoIIID,
           TableauVPLobattoIIIE,
           TableauVPLobattoIIIF,
           TableauVPLobattoIIIG,
           TableauVPSRK3,
           TableauVPLobattoIIIAIIIA,
           TableauVPLobattoIIIBIIIB,
           TableauVPRadauIIAIIA

    export TableauSymplecticProjection,
           TableauLobattoIIIAIIIBpSymplectic,
           TableauLobattoIIIBIIIApSymplectic,
           TableauGLRKpSymplectic

    include("tableaus/tableaus_vpark.jl")

    export TableauVSPARKLobattoIIIAIIIBProjection,
           TableauVSPARKLobattoIIIBIIIAProjection,
           TableauVSPARKModifiedLobattoIIIAIIIBProjection,
           TableauVSPARKModifiedLobattoIIIBIIIAProjection,
           TableauVSPARKInternalProjection,
           TableauVSPARKModifiedInternalProjection,
           TableauVSPARKMidpointProjection,
           TableauVSPARKModifiedMidpointProjection,
           TableauVSPARKSymmetricProjection,
           TableauVSPARKSymplecticProjection,
           TableauVSPARKGLRKpLobattoIIIAIIIB,
           TableauVSPARKGLRKpLobattoIIIBIIIA,
           TableauVSPARKGLRKpModifiedLobattoIIIAIIIB,
           TableauVSPARKGLRKpModifiedLobattoIIIBIIIA,
           TableauVSPARKGLRKpInternal,
           TableauVSPARKGLRKpModifiedInternal,
           TableauVSPARKGLRKpMidpoint,
           TableauVSPARKGLRKpModifiedMidpoint,
           TableauVSPARKGLRKpSymmetric,
           TableauVSPARKGLRKpSymplectic,
           TableauVSPARKLobattoIIIAIIIBpLobattoIIIAIIIB,
           TableauVSPARKLobattoIIIBIIIApLobattoIIIAIIIB,
           TableauVSPARKLobattoIIIAIIIBpLobattoIIIBIIIA,
           TableauVSPARKLobattoIIIBIIIApLobattoIIIBIIIA,
           TableauVSPARKLobattoIIIAIIIBpModifiedLobattoIIIAIIIB,
           TableauVSPARKLobattoIIIAIIIBpModifiedLobattoIIIBIIIA,
           TableauVSPARKLobattoIIIBIIIApModifiedLobattoIIIAIIIB,
           TableauVSPARKLobattoIIIBIIIApModifiedLobattoIIIBIIIA,
           TableauVSPARKLobattoIIIAIIIBpMidpoint,
           TableauVSPARKLobattoIIIBIIIApMidpoint,
           TableauVSPARKLobattoIIIAIIIBpModifiedMidpoint,
           TableauVSPARKLobattoIIIBIIIApModifiedMidpoint,
           TableauVSPARKLobattoIIIAIIIBpSymmetric,
           TableauVSPARKLobattoIIIBIIIApSymmetric,
           TableauVSPARKLobABCCD,
           TableauVSPARKLobABCCE,
           TableauVSPARKLobABDE,
           TableauVSPARKLobABED,
           TableauVSPARKLobABD,
           TableauVSPARKLobABE,
           TableauVSPARKLobDE,
           TableauVSPARKLobED

    include("tableaus/tableaus_vspark_primary.jl")

    export TableauVSPARKLobattoIIIAB,
           TableauVSPARKLobattoIIIBA,
           TableauVSPARKLobattoIIICC̄,
           TableauVSPARKLobattoIIIC̄C,
           TableauVSPARKLobattoIIID,
           TableauVSPARKLobattoIIIE,
           TableauVSPARKGLRKLobattoIIIAB,
           TableauVSPARKGLRKLobattoIIIBA,
           TableauVSPARKGLRKLobattoIIICC̄,
           TableauVSPARKGLRKLobattoIIIC̄C,
           TableauVSPARKGLRKLobattoIIID,
           TableauVSPARKGLRKLobattoIIIE

    include("tableaus/tableaus_vspark_secondary.jl")

    export getTableauHPARK,
           TableauHPARKGLRK,
           TableauHPARKLobattoIIIAIIIB,
           TableauHPARKLobattoIIIBIIIA

    include("tableaus/tableaus_hpark.jl")

    export TableauHSPARKSymmetricProjection,
           TableauHSPARKGLRKpSymmetric,
           TableauHSPARKLobattoIIIAIIIBpSymmetric,
           TableauHSPARKLobattoIIIBIIIApSymmetric

    include("tableaus/tableaus_hspark_primary.jl")

    export TableauHSPARKLobattoIIIAB,
           TableauHSPARKLobattoIIIBA,
           TableauHSPARKLobattoIIICC̄,
           TableauHSPARKLobattoIIIC̄C,
           TableauHSPARKLobattoIIID,
           TableauHSPARKLobattoIIIE,
           TableauHSPARKGLRKLobattoIIIAB,
           TableauHSPARKGLRKLobattoIIIBA,
           TableauHSPARKGLRKLobattoIIICC̄,
           TableauHSPARKGLRKLobattoIIIC̄C,
           TableauHSPARKGLRKLobattoIIID,
           TableauHSPARKGLRKLobattoIIIE

    include("tableaus/tableaus_hspark_secondary.jl")

    export TableauSLRKLobattoIIIAB,
           TableauSLRKLobattoIIIBA,
           TableauSLRKLobattoIIICC̄,
           TableauSLRKLobattoIIIC̄C,
           TableauSLRKLobattoIIID,
           TableauSLRKLobattoIIIE

    include("tableaus/tableaus_slrk.jl")

    export TableauLieA, TableauLieB, TableauStrang,
           TableauMcLachlan2, TableauMcLachlan4,
           TableauTripleJump, TableauSuzukiFractal

    include("tableaus/tableaus_splitting.jl")


    function __init__()
        add_config(:tab_compensated_summation, true)
    end

end
