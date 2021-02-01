module Integrators

    using CompactBasisFunctions
    using DelimitedFiles
    using Documenter: @doc
    using LinearAlgebra
    using OffsetArrays
    using QuadratureRules
    using SimpleSolvers

    using ..Common
    using ..Config
    using ..Interpolation
    using ..Utils

    using ..Discontinuities
    using ..Equations
    using ..Solutions


    export InitialGuess, InitialGuessODE, InitialGuessIODE, InitialGuessPODE,
           initialize!


    import ..Equations: _get_v̄, _get_f̄


    include("integrators/initial_guess/extrapolation.jl")
    include("integrators/initial_guess/initial_guess_ode.jl")
    include("integrators/initial_guess/initial_guess_iode.jl")
    include("integrators/initial_guess/initial_guess_pode.jl")


    export AbstractCoefficients, AbstractTableau

    include("integrators/abstract_coefficients.jl")
    include("integrators/abstract_tableau.jl")


    export Integrator, DeterministicIntegrator, StochasticIntegrator
    export ODEIntegrator, DAEIntegrator, SDEIntegrator,
           PODEIntegrator, PDAEIntegrator, PSDEIntegrator,
           IODEIntegrator, IDAEIntegrator,
           HODEIntegrator, HDAEIntegrator,
           VODEIntegrator, VDAEIntegrator,
           SPSDEIntegrator

    export IntegratorCache, IntegratorConstructor
    export integrate, integrate!, integrate_step!, equation, timestep
    export function_stages!#, NonlinearFunctionParameters

    include("integrators/abstract_integrator.jl")
    include("integrators/integrator_cache.jl")
    include("integrators/integrators_common.jl")


    export IntegratorExplicitEuler
    
    include("integrators/various/integrators_explicit_euler.jl")


    export CoefficientsRK, HeaderCoefficientsRK,
           AbstractTableauRK, AbstractTableauIRK, AbstractTableauPRK,
           IntegratorRK, writeTableauToFile

    export get_symplectic_conjugate_coefficients, symplecticize,
           check_symplecticity, symplecticity_conditions, 
           check_symmetry, compute_symplecticity_error,
           check_order_conditions_B, check_order_conditions_C, check_order_conditions_D

    include("integrators/rk/abstract_integrator_rk.jl")
    include("integrators/rk/coefficients.jl")
    include("integrators/rk/tableaus.jl")


    export IntegratorERK, TableauERK, readTableauERKFromFile

    include("integrators/rk/integrators_erk.jl")


    export IntegratorDIRK, TableauDIRK
    export IntegratorFIRK, TableauFIRK
    export IntegratorFIRKimplicit, IntegratorMidpointImplicit, IntegratorSRKimplicit

    include("integrators/rk/integrators_dirk.jl")
    include("integrators/rk/integrators_firk.jl")
    include("integrators/rk/integrators_firk_implicit.jl")
    include("integrators/rk/integrators_midpoint_implicit.jl")
    include("integrators/rk/integrators_srk_implicit.jl")


    export IntegratorEPRK, TableauEPRK
    export IntegratorIPRK, TableauIPRK
    export IntegratorFLRK

    include("integrators/rk/integrators_eprk.jl")
    include("integrators/rk/integrators_iprk.jl")
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


    include("integrators/SPARK.jl")
    include("integrators/VPRK.jl")


    include("integrators/integrators.jl")


    function __init__()
        default_params = (
            (:ig_interpolation, HermiteInterpolation),
            (:ig_extrapolation_stages, 5),
            (:int_show_progress_nmin,  1000),
        )

        for param in default_params
            add_config(param...)
        end
    end

end
