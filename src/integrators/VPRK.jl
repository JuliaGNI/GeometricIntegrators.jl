module VPRK

    using Documenter: @doc
    using LinearAlgebra: mul!
    using RungeKutta
    using SimpleSolvers

    import RungeKutta: AbstractTableau
    
    using ..GeometricBase
    using ..GeometricEquations

    using ..Methods

    import ..Solutions: SolutionStepPODE, SolutionStepPDAE, SolutionPDAE

    import ..Integrators

    import ..Integrators: Integrator, PDAEIntegrator, Newton
    import ..Integrators: InitialGuess, Extrapolation, HermiteExtrapolation
    import ..Integrators: initialguess!, initial_guess!, integrate_step!, residual!
    import ..Integrators: IODEIntegrator, IODEIntegratorCache,
                          AbstractIntegratorIRK, AbstractIntegratorPRK
    import ..Integrators: CacheDict, Cache, CacheType, nlsolution
    import ..Integrators: IntegratorCache
    import ..Integrators: AbstractCoefficients, CoefficientsPGLRK
    import ..Integrators: create_internal_stage_vector, 
                          update_vector_fields!, update_multiplier!,
                          initialize!, update!
    import ..Integrators: equation, equations, tableau, timestep
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
    # include("vprk/integrators_vprk_psecondary.jl")
    # include("vprk/integrators_vprk_pvariational.jl")
    # include("vprk/integrators_vprk_plegendre.jl")
    # include("vprk/integrators_vprk_ptableau.jl")

end
