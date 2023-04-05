module SPARK

    using Documenter: @doc
    using LinearAlgebra: Diagonal
    using RungeKutta
    using SimpleSolvers

    import RungeKutta: AbstractTableau, nstages, eachstage

    using ..GeometricBase
    using ..GeometricEquations
    using ..Config

    using ..Methods

    import ..Methods: tableau

    import ..Solutions: SolutionStepPDAE, SolutionVector, update!

    import ..Integrators

    import ..Integrators: Integrator, PDAEIntegrator, Parameters, Newton
    import ..Integrators: InitialGuess, InitialGuessIODE, InitialGuessPODE, Extrapolation, HermiteExtrapolation
    import ..Integrators: initialguess!, initial_guess!, integrate_step!, function_stages!
    import ..Integrators: CacheDict, Cache, CacheType, IDAEIntegratorCache
    import ..Integrators: AbstractCoefficients,
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


    include("spark/abstract_integrator_spark.jl")
    include("spark/coefficients.jl")

    include("spark/integrators_spark_cache.jl")
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

end
