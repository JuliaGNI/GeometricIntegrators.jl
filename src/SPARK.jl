module SPARK

    using CompactBasisFunctions
    using Documenter: @doc
    using GenericLinearAlgebra
    using LinearAlgebra: Diagonal
    using QuadratureRules
    using RungeKutta
    using RungeKutta.Tableaus
    using SimpleSolvers

    using GeometricBase
    using GeometricEquations

    import GeometricBase: tableau
    import GeometricBase.Utils: @define
    import RungeKutta: AbstractTableau, Tableau, nstages, eachstage


    import ..Integrators

    import ..Integrators: GeometricIntegrator, Newton
    import ..Integrators: HDAEMethod, IDAEMethod, LDAEMethod, PDAEMethod
    import ..Integrators: SolutionStepPDAE
    import ..Integrators: StageVector
    import ..Integrators: Extrapolation, HermiteExtrapolation
    import ..Integrators: initial_guess!, integrate_step!, residual!, solutionstep!
    import ..Integrators: CacheDict, Cache, CacheType, IDAEIntegratorCache
    import ..Integrators: AbstractCoefficients,
                          @CoefficientsRK, @HeaderCoefficientsRK
    import ..Integrators: create_internal_stage_vector,
                          update!, update_vector_fields!, update_multiplier!,
                          initialize!, initsolver, internal, nlsolution,
                          internal_variables
    import ..Integrators: current, equation, equations, timestep, eachstage, nstages
    import ..Integrators: cache, caches, iguess, method, problem, solver


    export CoefficientsARK, CoefficientsPRK, CoefficientsMRK, CoefficientsIRK,
           CoefficientsSPARK

    export AbstractIntegratorSPARK
    export AbstractTableauSPARK, TableauSPARK

    export SPARKMethod, IntegratorSPARK, TableauSPARK

    export HPARK, IntegratorHPARK, TableauHPARK
    export VPARK, IntegratorVPARK, TableauVPARK

    export HSPARK, IntegratorHSPARK, TableauHSPARK
    export HSPARKprimary, IntegratorHSPARKprimary, TableauHSPARKprimary
    export HSPARKsecondary, IntegratorHSPARKsecondary, TableauHSPARKsecondary

    export VSPARK, IntegratorVSPARK, TableauVSPARK
    export VSPARKprimary, IntegratorVSPARKprimary, TableauVSPARKprimary
    export VSPARKsecondary, IntegratorVSPARKsecondary, TableauVSPARKsecondary

    export IntegratorSLRK, SLRK


    include("spark/abstract.jl")
    include("spark/coefficients.jl")
    include("spark/cache.jl")

    include("spark/integrators_spark_common.jl")
    include("spark/integrators_spark_tableau.jl")
    include("spark/integrators_spark.jl")
    
    include("spark/integrators_vpark.jl")
    include("spark/integrators_hpark.jl")

    include("spark/integrators_vspark_common.jl")
    include("spark/integrators_vspark.jl")
    include("spark/integrators_vspark_primary.jl")
    include("spark/integrators_vspark_secondary.jl")

    include("spark/integrators_hspark_common.jl")
    include("spark/integrators_hspark.jl")
    include("spark/integrators_hspark_primary.jl")
    include("spark/integrators_hspark_secondary.jl")

    include("spark/integrators_slrk.jl")



    include("spark/coefficients_glrk.jl")
    include("spark/coefficients_lob.jl")


    include("spark/tableaus_spark.jl")

    export SPARKGLRK,
           SPARKLobattoIIIAIIIB,
           SPARKLobattoIIIBIIIA,
           SPARKGLRKLobattoIIIAIIIB,
           SPARKGLRKLobattoIIIBIIIA,
           SPARKLobatto,
           SPARKLobABC,
           SPARKLobABD,
           SPARKVPRK,
           SPARKGLVPRK

    export TableauSymplecticProjection,
           TableauLobattoIIIAIIIBpSymplectic,
           TableauLobattoIIIBIIIApSymplectic,
           TableauGausspSymplectic

    include("spark/tableaus_vpark.jl")

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

    include("spark/tableaus_vspark_primary.jl")

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

    include("spark/tableaus_vspark_secondary.jl")

    export getTableauHPARK,
           TableauHPARKGLRK,
           TableauHPARKLobattoIIIAIIIB,
           TableauHPARKLobattoIIIBIIIA

    include("spark/tableaus_hpark.jl")

    export TableauHSPARKSymmetricProjection,
           TableauHSPARKGLRKpSymmetric,
           TableauHSPARKLobattoIIIAIIIBpSymmetric,
           TableauHSPARKLobattoIIIBIIIApSymmetric

    include("spark/tableaus_hspark_primary.jl")

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

    include("spark/tableaus_hspark_secondary.jl")

    export SLRKLobattoIIIAB,
           SLRKLobattoIIIBA,
           SLRKLobattoIIICC̄,
           SLRKLobattoIIIC̄C,
           SLRKLobattoIIID,
           SLRKLobattoIIIE

    include("spark/tableaus_slrk.jl")

end
