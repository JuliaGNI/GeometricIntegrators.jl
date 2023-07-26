module Integrators

    using Reexport

    @reexport using GeometricBase
    
    using Documenter: @doc
    using ForwardDiff
    using GeometricBase
    # using GeometricBase.Config
    using GeometricEquations
    using GeometricSolutions
    using LinearAlgebra
    using OffsetArrays
    using QuadratureRules
    using SimpleSolvers

    using ..Discontinuities
    using ..Extrapolators
    using ..Methods
    using ..Solutions

    import Base: Callable

    import CompactBasisFunctions
    import CompactBasisFunctions: Basis
    import CompactBasisFunctions: nbasis

    import GeometricBase: description, reference, nconstraints, tableau, reset!
    import GeometricBase.Utils: @big, @define, compensated_summation
    
    import RungeKutta
    import RungeKutta: eachstage, nstages

    import ..Methods: hasnullvector, initmethod, implicit_update, nullvector

    import ..Solutions: current, cut_periodic_solution!, update!


    # compat workaroung
    Base.ndims(prob::GeometricProblem) = length(vec(prob.ics.q))


    export InitialGuess, NoInitialGuess
    export initialguess!

    include("initial_guess/initial_guess.jl")
    include("initial_guess/hermite.jl")
    include("initial_guess/midpoint.jl")


    export AbstractCoefficients

    include("integrators/abstract_coefficients.jl")


    export AbstractIntegrator, DeterministicIntegrator, StochasticIntegrator, Integrator
    export ODEIntegrator, DAEIntegrator, SDEIntegrator,
           PODEIntegrator, PDAEIntegrator, PSDEIntegrator,
           IODEIntegrator, IDAEIntegrator,
           HODEIntegrator, HDAEIntegrator,
           LODEIntegrator, LDAEIntegrator,
           SPSDEIntegrator

    export IntegratorCache, IntegratorConstructor
    export equation, timestep

    include("integrators/integrator_cache.jl")
    include("integrators/abstract_integrator.jl")


    export Integrator
    export integrate, integrate!, integrate_step!
    
    export NoSolver, NoProjection
    export default_solver, default_iguess, default_projection

    include("integrators/integrator.jl")


    export ExplicitEuler, ImplicitEuler
    
    include("integrators/various/integrators_explicit_euler.jl")
    include("integrators/various/integrators_implicit_euler.jl")


    export AbstractIntegratorRK, AbstractIntegratorIRK, AbstractIntegratorPRK, IntegratorRK

    export get_symplectic_conjugate_coefficients, symplecticize,
           check_symplecticity, symplecticity_conditions, 
           check_symmetry, compute_symplecticity_error

    include("integrators/rk/abstract.jl")
    include("integrators/rk/common.jl")
    include("integrators/rk/updates.jl")
    include("integrators/rk/tableaus.jl")


    export IntegratorERK
    export IntegratorIRK
    export IntegratorIRKimplicit
    export IntegratorDIRK
    #export IntegratorMidpointImplicit, IntegratorSRKimplicit

    include("integrators/rk/integrators_erk.jl")
    include("integrators/rk/integrators_irk.jl")
    include("integrators/rk/integrators_irk_implicit.jl")
    include("integrators/rk/integrators_dirk.jl")
    # include("integrators/rk/integrators_midpoint_implicit.jl")
    # include("integrators/rk/integrators_srk_implicit.jl")


    export IntegratorEPRK
    export IntegratorIPRK
    export IntegratorIPRKimplicit
    # export IntegratorFLRK
    # export IntegratorPGLRK, CoefficientsPGLRK

    include("integrators/rk/integrators_eprk.jl")
    include("integrators/rk/integrators_iprk.jl")
    include("integrators/rk/integrators_iprk_implicit.jl")
    # include("integrators/rk/integrators_flrk.jl")
    include("integrators/rk/pglrk_coefficients.jl")
    # include("integrators/rk/pglrk_integrators.jl")


    export ExactSolution,
           IntegratorExactODE,
           IntegratorSplitting,
           IntegratorComposition,
           AbstractTableauSplitting,
           Composition,
           Splitting,
           SplittingCoefficientsGeneral,
           SplittingCoefficientsNonSymmetric,
           SplittingCoefficientsGS,
           SplittingCoefficientsSS

    export LieA,
           LieB,
           Strang,
           Marchuk,
           StrangA,
           StrangB,
           McLachlan2,
           McLachlan4,
           TripleJump,
           SuzukiFractal
           
    include("integrators/splitting/exact_solution.jl")
    include("integrators/splitting/splitting_coefficients.jl")
    include("integrators/splitting/splitting_methods.jl")
    include("integrators/splitting/splitting_integrator.jl")
    include("integrators/splitting/composition_methods.jl")
    include("integrators/splitting/composition_integrator.jl")


    export IntegratorVPRK

    include("integrators/vi/integrators_vprk_cache.jl")
    include("integrators/vi/integrators_vprk.jl")


    # export IntegratorCGVI, IntegratorDGVI, IntegratorDGVIEXP,
    #        IntegratorDGVIPI, IntegratorDGVIP0, IntegratorDGVIP1

    # include("integrators/cgvi/integrators_cgvi.jl")
    # include("integrators/dgvi/integrators_dgvi.jl")
    # include("integrators/dgvi/integrators_dgvi_experimental.jl")
    # include("integrators/dgvi/integrators_dgvi_path_integral.jl")
    # include("integrators/dgvi/integrators_dgvi_projection_initial.jl")
    # include("integrators/dgvi/integrators_dgvi_projection_final.jl")


    export IntegratorDVIA, IntegratorDVIB,
           IntegratorCMDVI, IntegratorCTDVI,
           IntegratorDVRK
    
    include("integrators/dvi/integrators_dvi_a.jl")
    include("integrators/dvi/integrators_dvi_b.jl")
    include("integrators/dvi/integrators_cmdvi.jl")
    include("integrators/dvi/integrators_ctdvi.jl")
    include("integrators/dvi/integrators_dvrk.jl")


    include("integrators/integrators.jl")
    include("integrators/methods.jl")


    include("projections/cache.jl")
    include("projections/common.jl")
    include("projections/midpoint_projection.jl")
    include("projections/standard_projection.jl")
    include("projections/symmetric_projection.jl")
    

    # function __init__()
    #     default_params = (
    #         (:ig_extrapolation, HermiteExtrapolation),
    #         (:ig_extrapolation_stages, 5),
    #         (:int_show_progress_nmin,  1000),
    #         (:tab_compensated_summation, true),
    #     )

    #     for param in default_params
    #         add_config(param...)
    #     end
    # end

end
