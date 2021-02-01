module Tableaus

    import GenericLinearAlgebra
    
    using CompactBasisFunctions
    using QuadratureRules
    using Reexport

    @reexport using RungeKutta


    using ..Config
    using ..Common
    using ..Integrators
    using ..Integrators.SPARK
    using ..Integrators.VPRK
    using ..Utils

    using GeometricIntegrators.Integrators.SPARK: get_ã_vspark_primary,
                                                  get_α_vspark_primary,
                                                  compute_ã_vspark_primary,
                                                  compute_α_vspark_primary

    export TableauGauss

    
    include("tableaus/coefficients_glrk.jl")
    include("tableaus/coefficients_lob.jl")


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
           TableauVPLobattoIIIC̄,
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
           TableauGausspSymplectic

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

    export TableauLieA,
           TableauLieB,
           TableauMarchuk,
           TableauMcLachlan2,
           TableauMcLachlan4,
           TableauStrang,
           TableauStrangA,
           TableauStrangB,
           TableauSuzukiFractal,
           TableauTripleJump

    include("tableaus/tableaus_splitting.jl")


    function __init__()
        add_config(:tab_compensated_summation, true)
    end

end
