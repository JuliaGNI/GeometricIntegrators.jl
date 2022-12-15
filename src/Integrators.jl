module Integrators

    using Reexport

    @reexport using GeometricBase
    
    using CompactBasisFunctions
    using Documenter: @doc
    using LinearAlgebra
    using OffsetArrays
    using QuadratureRules
    using RungeKutta
    using SimpleSolvers

    using GeometricBase.Config
    using GeometricBase.Utils
    using GeometricEquations
    using GeometricEquations: _get_v, _get_f, _get_v̄, _get_f̄
    using GeometricSolutions

    using ..Discontinuities
    using ..Solutions


    import GeometricBase: timestep

    import RungeKutta: AbstractTableau, nstages


    # compat workaroung
    Base.ndims(prob::GeometricProblem) = length(vec(prob.ics.q))


    export Extrapolation,
           EulerExtrapolation,
           MidpointExtrapolation,
           HermiteExtrapolation
    export extrapolate!

    include("integrators/extrapolation/extrapolation.jl")
    include("integrators/extrapolation/aitken_neville.jl")
    include("integrators/extrapolation/euler.jl")
    include("integrators/extrapolation/hermite.jl")
    include("integrators/extrapolation/midpoint.jl")


    export InitialGuess, InitialGuessODE, InitialGuessIODE, InitialGuessPODE,
           initialize!

    include("integrators/initial_guess/initial_guess_ode.jl")
    include("integrators/initial_guess/initial_guess_iode.jl")
    include("integrators/initial_guess/initial_guess_pode.jl")


    export AbstractCoefficients

    include("integrators/abstract_coefficients.jl")
    include("integrators/abstract_tableau.jl")


    export Integrator, DeterministicIntegrator, StochasticIntegrator
    export ODEIntegrator, DAEIntegrator, SDEIntegrator,
           PODEIntegrator, PDAEIntegrator, PSDEIntegrator,
           IODEIntegrator, IDAEIntegrator,
           HODEIntegrator, HDAEIntegrator,
           LODEIntegrator, LDAEIntegrator,
           SPSDEIntegrator

    export IntegratorCache, IntegratorConstructor
    export integrate, integrate!, integrate_step!, equation, timestep
    export function_stages!#, NonlinearFunctionParameters

    include("integrators/abstract_integrator.jl")
    include("integrators/integrator_cache.jl")
    include("integrators/integrators_common.jl")


    export IntegratorExplicitEuler
    
    include("integrators/various/integrators_explicit_euler.jl")


    export AbstractIntegratorRK, AbstractIntegratorIRK, AbstractIntegratorPRK, IntegratorRK

    export get_symplectic_conjugate_coefficients, symplecticize,
           check_symplecticity, symplecticity_conditions, 
           check_symmetry, compute_symplecticity_error

    include("integrators/rk/abstract_integrator_rk.jl")
    include("integrators/rk/tableaus.jl")


    export IntegratorERK
    export IntegratorDIRK
    export IntegratorFIRK
    export IntegratorFIRKimplicit, IntegratorMidpointImplicit, IntegratorSRKimplicit

    include("integrators/rk/integrators_erk.jl")
    include("integrators/rk/integrators_dirk.jl")
    include("integrators/rk/integrators_firk.jl")
    include("integrators/rk/integrators_firk_implicit.jl")
    include("integrators/rk/integrators_midpoint_implicit.jl")
    include("integrators/rk/integrators_srk_implicit.jl")


    export IntegratorEPRK
    export IntegratorIPRK
    export IntegratorPRKimplicit
    export IntegratorFLRK

    include("integrators/rk/integrators_eprk.jl")
    include("integrators/rk/integrators_iprk.jl")
    include("integrators/rk/integrators_prk_implicit.jl")
    include("integrators/rk/integrators_flrk.jl")


    export IntegratorPGLRK, CoefficientsPGLRK

    include("integrators/rk/integrators_pglrk.jl")


    export IntegratorSplitting,
           IntegratorComposition,
           IntegratorExactODE,
           AbstractTableauSplitting,
           TableauSplitting,
           TableauSplittingGS,
           TableauSplittingNS,
           TableauSplittingSS
           
    include("integrators/splitting/integrators_exact_ode.jl")
    include("integrators/splitting/splitting_tableau.jl")
    include("integrators/splitting/integrators_composition.jl")
    include("integrators/splitting/integrators_splitting.jl")


    export IntegratorCGVI, IntegratorDGVI, IntegratorDGVIEXP,
           IntegratorDGVIPI, IntegratorDGVIP0, IntegratorDGVIP1

    include("integrators/cgvi/integrators_cgvi.jl")
    include("integrators/dgvi/integrators_dgvi.jl")
    include("integrators/dgvi/integrators_dgvi_experimental.jl")
    include("integrators/dgvi/integrators_dgvi_path_integral.jl")
    include("integrators/dgvi/integrators_dgvi_projection_initial.jl")
    include("integrators/dgvi/integrators_dgvi_projection_final.jl")


    export IntegratorDVIA, IntegratorDVIB,
           IntegratorCMDVI, IntegratorCTDVI
    
    include("integrators/dvi/integrators_dvi_a.jl")
    include("integrators/dvi/integrators_dvi_b.jl")
    include("integrators/dvi/integrators_cmdvi.jl")
    include("integrators/dvi/integrators_ctdvi.jl")


    include("integrators/SPARK.jl")
    include("integrators/VPRK.jl")


    include("integrators/integrators.jl")


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
