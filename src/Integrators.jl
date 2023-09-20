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
    using PrettyTables
    using QuadratureRules
    using RungeKutta.Tableaus
    using RungeKutta.PartitionedTableaus
    using SimpleSolvers

    using ..Discontinuities
    using ..Extrapolators

    import Base: Callable

    import CompactBasisFunctions
    import CompactBasisFunctions: Basis
    import CompactBasisFunctions: nbasis

    import GeometricBase: description, reference, tableau, order
    import GeometricBase: equations, nconstraints, parameters, timestep
    import GeometricBase: integrate, integrate!
    import GeometricBase: reset!
    import GeometricBase.Utils: @big, @define, compensated_summation
    
    import RungeKutta
    import RungeKutta: eachstage, nstages
    import RungeKutta: AbstractTableau, Tableau, PartitionedTableau, SymplecticTableau, SymplecticPartitionedTableau

    import SimpleSolvers: SolverMethod


    # compat workaround
    Base.ndims(prob::EquationProblem) = length(vec(prob.ics.q))


    const DEFAULT_NSAVE = 1
    const DEFAULT_NWRITE = 0


    export AbstractCoefficients

    include("integrators/abstract_coefficients.jl")


    export GeometricMethod
    export ODEMethod, PODEMethod, HODEMethod, IODEMethod, LODEMethod, SODEMethod
    export DAEMethod, PDAEMethod, HDAEMethod, IDAEMethod, LDAEMethod
    export RKMethod, ERKMethod, IRKMethod, DIRKMethod
    export PRKMethod, EPRKMethod, IPRKMethod, VPRKMethod, DVIMethod
    export AbstractSplittingMethod

    export RK, PRK

    export default_solver, default_iguess, default_projection
    export initmethod, implicit_update
    export internal_variables

    export AbstractTableau, Tableau, PartitionedTableau, SymplecticTableau, SymplecticPartitionedTableau
    export name, order, description, reference
    export nstages, eachstage, coefficients, weights, nodes
    export hasnullvector, nullvector
    export isexplicit, isimplicit, issymmetric, issymplectic, isenergypreserving, isstifflyaccurate

    include("integrators/methods.jl")


    export DataSeries, TimeSeries, Solution, AbstractSolution, DeterministicSolution
    export current, previous, history

    export SolutionODE, SolutionPODE
    export SolutionDAE, SolutionPDAE

    export SolutionStep,
        SolutionStepODE, SolutionStepPODE,
        SolutionStepDAE, SolutionStepPDAE

    export reset!, update!, update_vector_fields!, cut_periodic_solution!

    include("solutions/solution_step.jl")
    include("solutions/solution.jl")
    include("solutions/solution_step_ode.jl")
    include("solutions/solution_step_pode.jl")
    include("solutions/solution_step_dae.jl")
    include("solutions/solution_step_pdae.jl")
    include("solutions/solution_step_constructors.jl")


    export InitialGuess, NoInitialGuess
    export initialguess!

    include("initial_guess/initial_guess.jl")
    include("initial_guess/hermite.jl")
    include("initial_guess/midpoint.jl")


    export IntegratorCache, IntegratorConstructor
    export equation, timestep

    include("integrators/integrator_cache.jl")


    export NoSolver

    include("integrators/solvers.jl")


    export GeometricIntegrator
    export integrate, integrate!, integrate_step!

    include("integrators/integrator.jl")


    export ExplicitEuler, ImplicitEuler
    
    include("integrators/euler/explicit_euler.jl")
    include("integrators/euler/implicit_euler.jl")


    export AbstractIntegratorRK, AbstractIntegratorIRK, AbstractIntegratorPRK, IntegratorRK
    
    export ERKIntegrator
    export IntegratorIRK
    export IntegratorIRKimplicit
    export IntegratorDIRK
    #export IntegratorMidpointImplicit, IntegratorSRKimplicit

    export IntegratorEPRK
    export IPRKIntegrator
    export IntegratorIPRKimplicit
    # export IntegratorFLRK
    # export IntegratorPGLRK, CoefficientsPGLRK

    export get_symplectic_conjugate_coefficients, symplecticize,
           check_symplecticity, symplecticity_conditions, 
           check_symmetry, compute_symplecticity_error

    include("integrators/rk/abstract.jl")
    include("integrators/rk/common.jl")
    include("integrators/rk/updates.jl")
    include("integrators/rk/tableaus.jl")

    include("integrators/rk/integrators_erk.jl")
    include("integrators/rk/integrators_irk.jl")
    include("integrators/rk/integrators_irk_implicit.jl")
    include("integrators/rk/integrators_dirk.jl")
    # include("integrators/rk/integrators_midpoint_implicit.jl")
    # include("integrators/rk/integrators_srk_implicit.jl")

    include("integrators/rk/integrators_eprk.jl")
    include("integrators/rk/integrators_iprk.jl")
    include("integrators/rk/integrators_iprk_implicit.jl")
    # include("integrators/rk/integrators_flrk.jl")
    include("integrators/rk/pglrk_coefficients.jl")
    # include("integrators/rk/pglrk_integrators.jl")

    include("integrators/rk/methods.jl")


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

    include("integrators/splitting/exact_solution.jl")
    include("integrators/splitting/splitting_coefficients.jl")
    include("integrators/splitting/splitting_methods.jl")
    include("integrators/splitting/splitting_integrator.jl")
    include("integrators/splitting/composition_methods.jl")
    include("integrators/splitting/composition_integrator.jl")


    export IntegratorVPRK

    include("integrators/vi/vi_methods.jl")
    include("integrators/vi/position_momentum_common.jl")
    include("integrators/vi/position_momentum_cache.jl")
    include("integrators/vi/position_momentum_midpoint.jl")
    include("integrators/vi/position_momentum_trapezoidal.jl")
    include("integrators/vi/vprk_methods.jl")
    include("integrators/vi/vprk_cache.jl")
    include("integrators/vi/vprk_integrator.jl")


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
    
    include("integrators/dvi/dvi_common.jl")
    include("integrators/dvi/dvi_cache.jl")
    include("integrators/dvi/dvi_euler.jl")
    include("integrators/dvi/dvi_midpoint.jl")
    include("integrators/dvi/dvi_trapezoidal.jl")
    include("integrators/dvi/dvrk.jl")


    export HPImidpointIntegrator, HPItrapezoidalIntegrator

    include("integrators/hpi/hpi_common.jl")
    include("integrators/hpi/hpi_cache.jl")
    include("integrators/hpi/hpi_midpoint.jl")
    include("integrators/hpi/hpi_trapezoidal.jl")


    export NoProjection, projection

    include("projections/projection.jl")
    include("projections/cache.jl")
    include("projections/common.jl")
    include("projections/midpoint_projection.jl")
    include("projections/standard_projection.jl")
    include("projections/symmetric_projection.jl")

    include("projections/methods.jl")
    include("integrators/vi/vprk_projected.jl")
    

    include("integrators/method_list.jl")


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
