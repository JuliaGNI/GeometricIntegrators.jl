module Integrators

    using Reexport

    @reexport using GeometricBase
    
    using Documenter: @doc
    using ForwardDiff
    using GeometricBase
    using LinearAlgebra
    using OffsetArrays
    using QuadratureRules
    using SimpleSolvers

    using GeometricBase.Config
    using GeometricBase.Utils
    using GeometricEquations
    using GeometricSolutions

    using ..Discontinuities
    using ..Extrapolators
    using ..Methods
    using ..Solutions

    import Base: Callable

    import CompactBasisFunctions
    import CompactBasisFunctions: Basis
    import CompactBasisFunctions: nbasis

    import GeometricEquations: nconstraints
    
    import RungeKutta
    import RungeKutta: eachstage, nstages

    import ..Extrapolators: extrapolate_ode!, extrapolate_iode!, extrapolate_pode!

    import ..Methods: hasnullvector, initmethod, implicit_update, nullvector, tableau

    import ..Solutions: current, cut_periodic_solution!, reset!


    # compat workaroung
    Base.ndims(prob::GeometricProblem) = length(vec(prob.ics.q))


    export InitialGuess, NoInitialGuess
    # export InitialGuessODE, InitialGuessIODE, InitialGuessPODE
    export initialguess!

    include("integrators/initial_guess/initial_guess.jl")
    include("integrators/initial_guess/hermite.jl")
    include("integrators/initial_guess/midpoint.jl")
    
    # include("integrators/initial_guess/initial_guess_ode.jl")
    # include("integrators/initial_guess/initial_guess_iode.jl")
    # include("integrators/initial_guess/initial_guess_pode.jl")


    export AbstractCoefficients

    include("integrators/abstract_coefficients.jl")
    include("integrators/abstract_tableau.jl")


    export AbstractIntegrator, DeterministicIntegrator, StochasticIntegrator, Integrator
    export ODEIntegrator, DAEIntegrator, SDEIntegrator,
           PODEIntegrator, PDAEIntegrator, PSDEIntegrator,
           IODEIntegrator, IDAEIntegrator,
           HODEIntegrator, HDAEIntegrator,
           LODEIntegrator, LDAEIntegrator,
           SPSDEIntegrator

    export IntegratorCache, IntegratorConstructor
    export equation, timestep
    # export function_stages!, NonlinearFunctionParameters

    include("integrators/abstract_integrator.jl")
    include("integrators/integrator_cache.jl")
    include("integrators/integrators_common.jl")


    export Integrator
    export integrate, integrate!, integrate_step!
    
    export NoSolver, NoProjection
    export default_solver, default_iguess, default_projection

    include("integrators/integrator.jl")


    export IntegratorExplicitEuler
    
    include("integrators/various/integrators_explicit_euler.jl")


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


    export IntegratorSplitting,
           IntegratorComposition,
           IntegratorExactODE,
           AbstractTableauSplitting,
           Composition,
           Splitting,
           SplittingCoefficientsGeneral,
           SplittingCoefficientsNonSymmetric,
           SplittingCoefficientsGS,
           SplittingCoefficientsSS

    export coefficients

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
           
    include("integrators/splitting/integrators_exact_ode.jl")
    include("integrators/splitting/splitting_coefficients.jl")
    include("integrators/splitting/splitting_methods.jl")
    include("integrators/splitting/integrators_composition.jl")
    include("integrators/splitting/integrators_splitting.jl")


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


    include("integrators/SPARK.jl")
    include("integrators/VPRK.jl")


    include("integrators/integrators.jl")
    include("integrators/methods.jl")


    include("projections/cache.jl")
    include("projections/common.jl")
    include("projections/midpoint_projection.jl")
    include("projections/standard_projection.jl")
    include("projections/symmetric_projection.jl")
    

    function __init__()
        default_params = (
            (:ig_extrapolation, HermiteExtrapolation),
            (:ig_extrapolation_stages, 5),
            (:int_show_progress_nmin,  1000),
        )

        for param in default_params
            add_config(param...)
        end
    end

end
