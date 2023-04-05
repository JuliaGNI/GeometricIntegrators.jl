module VPRK

    using Documenter: @doc
    using LinearAlgebra: mul!
    using RungeKutta
    using SimpleSolvers

    import RungeKutta: AbstractTableau
    
    using ..GeometricBase
    using ..GeometricEquations
    using ..Config
    using ..Utils

    using ..Methods

    import ..Solutions: SolutionStepPODE, SolutionStepPDAE, SolutionPDAE, SolutionVector

    import ..Integrators

    import ..Integrators: Integrator, PDAEIntegrator, Parameters, Newton
    import ..Integrators: InitialGuess, InitialGuessIODE, InitialGuessPODE, Extrapolation, HermiteExtrapolation
    import ..Integrators: initialguess!, initial_guess!, integrate_step!, function_stages!
    import ..Integrators: IODEIntegrator, IODEIntegratorCache,
                          AbstractIntegratorIRK, AbstractIntegratorPRK
    import ..Integrators: CacheDict, Cache, CacheType, nlsolution
    import ..Integrators: IntegratorCache, OldCacheDict, Parameters
    import ..Integrators: AbstractCoefficients, CoefficientsPGLRK,
                          @CoefficientsRK, @HeaderTableau, @HeaderCoefficientsRK
    import ..Integrators: create_internal_stage_vector, create_nonlinear_solver,
                          update_vector_fields!, update_multiplier!,
                          initialize!, update!
    import ..Integrators: equation, equations, tableau, timestep,
                          eachdim
    import ..Integrators: solver


    # export IntegratorVPRK, IntegratorVPRKpNone

    # export IntegratorVPRKdegenerate

    # export IntegratorVPRKpStandard, IntegratorVPRKpSymplectic,
    # export IntegratorVPRKpInternal, IntegratorVPRKpMidpoint,
    #        IntegratorVPRKpSymmetric, IntegratorVPRKpTableau,
        #    IntegratorVPRKpSecondary, IntegratorVPRKpVariational#,
        #    IntegratorVPRKpVariationalQ, IntegratorVPRKpVariationalP

    # export IntegratorVPRKpLegendre#, TableauVPRKpLegendre

    # include("vprk/integrators_vprk_abstract.jl")
    # include("vprk/integrators_vprk_cache.jl")
    # include("vprk/integrators_vprk_parameters.jl")
    # include("vprk/integrators_vprk_common.jl")
    # include("vprk/integrators_vprk.jl")
    # include("vprk/integrators_vprk_degenerate.jl")
    # include("vprk/integrators_vprk_pinternal.jl")
    # include("vprk/integrators_vprk_pmidpoint.jl")
    # include("vprk/integrators_vprk_pstandard.jl")
    # include("vprk/integrators_vprk_psecondary.jl")
    # include("vprk/integrators_vprk_psymmetric.jl")
    # include("vprk/integrators_vprk_pvariational.jl")
    # include("vprk/integrators_vprk_plegendre.jl")
    # include("vprk/integrators_vprk_ptableau.jl")

end
