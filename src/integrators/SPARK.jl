module SPARK

    using CompactBasisFunctions
    using Documenter: @doc
    using LinearAlgebra
    using QuadratureRules
    using RungeKutta
    using RungeKutta: AbstractTableau, @TableauHeader
    using SimpleSolvers

    using ..GeometricBase
    using ..GeometricEquations
    using ..Config

    import ..Solutions: SolutionStepODE, SolutionStepPDAE, SolutionVector, update!

    import ..Integrators

    import ..Integrators: AbstractIntegratorRK, PDAEIntegrator
    import ..Integrators: InitialGuessIODE, InitialGuessODE, InitialGuessPODE, Parameters
    import ..Integrators: IntegratorCache, IDAEIntegratorCache, ODEIntegratorCache, CacheDict, CacheType
    import ..Integrators: AbstractTableau,
                          @CoefficientsRK, @HeaderCoefficientsRK
    import ..Integrators: create_internal_stage_vector, create_nonlinear_solver,
                          update_vector_fields!, update_solution!, update_multiplier!,
                          initialize!
    import ..Integrators: equation, equations, tableau, timestep,
                          eachdim, eachstage, nstages

    import ..Utils: @define


    export CoefficientsARK, CoefficientsPRK, CoefficientsMRK, CoefficientsIRK,
           CoefficientsSPARK

    export AbstractIntegratorSPARK
    export AbstractTableauSPARK, TableauSPARK

    export IntegratorSPARK, TableauSPARK

    export IntegratorHPARK, TableauHPARK
    export IntegratorVPARK, TableauVPARK

    export IntegratorHSPARK, TableauHSPARK
    export IntegratorHSPARKprimary, TableauHSPARKprimary
    export IntegratorHSPARKsecondary, TableauHSPARKsecondary

    export IntegratorVSPARK, TableauVSPARK
    export IntegratorVSPARKprimary, TableauVSPARKprimary
    export IntegratorVSPARKsecondary, TableauVSPARKsecondary

    export IntegratorSLRK, TableauSLRK
    export IntegratorPGLRK, CoefficientsPGLRK


    include("spark/abstract_integrator_spark.jl")
    include("spark/coefficients.jl")

    include("spark/integrators_spark_cache.jl")
    include("spark/integrators_spark_common.jl")
    include("spark/integrators_spark_tableau.jl")
    include("spark/integrators_spark_parameters.jl")

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

    include("spark/integrators_pglrk.jl")
    include("spark/integrators_slrk.jl")

end
