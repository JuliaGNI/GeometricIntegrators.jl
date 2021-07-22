module VPRK

    using Documenter: @doc
    using LinearAlgebra: mul!
    using RungeKutta
    using SimpleSolvers

    import GeometricEquations: IODE, LODE, LDAE, get_functions, hassecondary

    using ..GeometricBase
    using ..Config
    using ..Utils

    import ..Solutions: AtomicSolutionPODE, AtomicSolutionPDAE, SolutionPDAE, SolutionVector
    import ..Solutions: update!

    import ..Integrators

    import ..Integrators: IODEIntegrator, IODEIntegratorCache, InitialGuessIODE,
                          AbstractIntegratorIRK, AbstractIntegratorPRK
    import ..Integrators: IntegratorCache, CacheDict, CacheType, Parameters
    import ..Integrators: AbstractTableau, AbstractCoefficients,
                          CoefficientsPGLRK,
                          @CoefficientsRK, @HeaderTableau, @HeaderCoefficientsRK,
                          get_symplectic_conjugate_coefficients
    import ..Integrators: create_internal_stage_vector, create_nonlinear_solver,
                          update_vector_fields!, update_solution!, update_multiplier!,
                          initialize!
    import ..Integrators: equation, equations, tableau, timestep,
                          eachdim, eachstage, nstages


    export IntegratorVPRK, IntegratorVPRKpNone, TableauVPRK

    export IntegratorVPRKdegenerate

    export IntegratorVPRKpStandard, IntegratorVPRKpSymplectic,
         IntegratorVPRKpInternal, IntegratorVPRKpMidpoint,
         IntegratorVPRKpSymmetric, IntegratorVPRKpTableau,
         IntegratorVPRKpSecondary, IntegratorVPRKpVariational,
         IntegratorVPRKpVariationalQ, IntegratorVPRKpVariationalP

    export IntegratorVPRKpLegendre#, TableauVPRKpLegendre

    include("vprk/integrators_vprk_abstract.jl")
    include("vprk/integrators_vprk_cache.jl")
    include("vprk/integrators_vprk_tableau.jl")
    include("vprk/integrators_vprk_parameters.jl")
    include("vprk/integrators_vprk_common.jl")
    include("vprk/integrators_vprk.jl")
    include("vprk/integrators_vprk_degenerate.jl")
    include("vprk/integrators_vprk_pinternal.jl")
    include("vprk/integrators_vprk_pmidpoint.jl")
    include("vprk/integrators_vprk_pstandard.jl")
    include("vprk/integrators_vprk_psecondary.jl")
    include("vprk/integrators_vprk_psymmetric.jl")
    include("vprk/integrators_vprk_ptableau.jl")
    include("vprk/integrators_vprk_pvariational.jl")
    include("vprk/integrators_vprk_plegendre.jl")

end
