module Integrators

using Reexport

@reexport using GeometricBase
@reexport using GeometricIntegratorsBase

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
using StaticArrays

using ..Discontinuities
# using ..Extrapolators

import Base: Callable

import CompactBasisFunctions
import CompactBasisFunctions: Basis
import CompactBasisFunctions: nbasis

import GeometricBase: description, reference, tableau, order
import GeometricBase: equations, initialguess, nconstraints, parameters, timestep
import GeometricBase: integrate, integrate!, reset!, solutionstep!, update!
import GeometricBase.Utils: @big, @define, compensated_summation

import GeometricIntegratorsBase: ODEMethod, PODEMethod, HODEMethod, IODEMethod, LODEMethod, SODEMethod, DELEMethod
import GeometricIntegratorsBase: DAEMethod, PDAEMethod, HDAEMethod, IDAEMethod, LDAEMethod
import GeometricIntegratorsBase: ODEIntegratorCache, IODEIntegratorCache, PODEIntegratorCache, DELEIntegratorCache
import GeometricIntegratorsBase: DAEIntegratorCache, IDAEIntegratorCache, PDAEIntegratorCache

import GeometricIntegratorsBase: Cache, CacheType, GeometricIntegrator, NoProjection, ProjectionMethod
import GeometricIntegratorsBase: cache, hasnullvector, iguess, internal, method, nlsolution, nullvector, problem, projection, solver
import GeometricIntegratorsBase: components!, initialize!, initial_guess!, integrate_step!, residual!
import GeometricIntegratorsBase: default_options, initialize!, initmethod, initsolver, internal_variables, copy_internal_variables!
import GeometricIntegratorsBase: Newton

import GeometricIntegratorsBase: isexplicit, isimplicit, issymmetric, issymplectic, isenergypreserving, isstifflyaccurate, implicit_update
import GeometricIntegratorsBase: isodemethod, ispodemethod, ishodemethod, isiodemethod, islodemethod, issodemethod
import GeometricIntegratorsBase: isdaemethod, ispdaemethod, ishdaemethod, isidaemethod, isldaemethod
import GeometricIntegratorsBase: default_solver, default_linesearch, default_iguess, solversize

import RungeKutta
import RungeKutta: eachstage, nstages
import RungeKutta: AbstractTableau, Tableau, PartitionedTableau, SymplecticTableau, SymplecticPartitionedTableau

import SimpleSolvers: SolverMethod


# compat workaround
# Base.ndims(prob::EquationProblem) = length(vec(prob.ics.q))


const DEFAULT_NSAVE = 1
const DEFAULT_NWRITE = 0


export AbstractCoefficients

include("integrators/abstract_coefficients.jl")


export GeometricMethod
export ODEMethod, PODEMethod, HODEMethod, IODEMethod, LODEMethod, SODEMethod
export DAEMethod, PDAEMethod, HDAEMethod, IDAEMethod, LDAEMethod
export DELEMethod, Midpoint, Trapezoidal
export RKMethod, ERKMethod, IRKMethod, DIRKMethod
export PRKMethod, EPRKMethod, IPRKMethod, VPRKMethod, DVIMethod
export AbstractSplittingMethod

export RK, PRK

export default_solver, default_iguess, default_projection
export initmethod, implicit_update
export internal_variables

export AbstractTableau, Tableau, PartitionedTableau, SymplecticTableau, SymplecticPartitionedTableau
export name, order, description, reference
# export nstages, eachstage, coefficients, weights, nodes
export hasnullvector, nullvector
export isexplicit, isimplicit, issymmetric, issymplectic, isenergypreserving, isstifflyaccurate
export isodemethod, ispodemethod, ishodemethod, isiodemethod, islodemethod, issodemethod
export isdaemethod, ispdaemethod, ishdaemethod, isidaemethod, isldaemethod

# include("integrators/methods.jl")


export DataSeries, TimeSeries, Solution, AbstractSolution, DeterministicSolution
export current, previous, history

export SolutionODE, SolutionPODE
export SolutionDAE, SolutionPDAE

export SolutionStep,
    SolutionStepODE, SolutionStepPODE,
    SolutionStepDAE, SolutionStepPDAE

export reset!, update!, update_vector_fields!, cut_periodic_solution!

# include("solutions/solution_step.jl")
# include("solutions/solution.jl")
# include("solutions/solution_step_ode.jl")
# include("solutions/solution_step_pode.jl")
# include("solutions/solution_step_dae.jl")
# include("solutions/solution_step_pdae.jl")
# include("solutions/solution_step_constructors.jl")


export equation, timestep


export NoSolver

export GeometricIntegrator
export integrate, integrate!, solutionstep!

include("integrators/integrator.jl")



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


include("integrators/vi/vi_methods.jl")
include("integrators/vi/deleqs.jl")
include("integrators/vi/deleqs_methods.jl")
include("integrators/vi/position_momentum_common.jl")
include("integrators/vi/position_momentum_cache.jl")
include("integrators/vi/position_momentum_midpoint.jl")
include("integrators/vi/position_momentum_trapezoidal.jl")
include("integrators/vi/vprk_methods.jl")
include("integrators/vi/vprk_cache.jl")
include("integrators/vi/vprk_integrator.jl")


include("integrators/cgvi/integrators_cgvi.jl")
# include("integrators/dgvi/integrators_dgvi.jl")
# include("integrators/dgvi/integrators_dgvi_experimental.jl")
# include("integrators/dgvi/integrators_dgvi_path_integral.jl")
# include("integrators/dgvi/integrators_dgvi_projection_initial.jl")
# include("integrators/dgvi/integrators_dgvi_projection_final.jl")


include("integrators/dvi/dvi_common.jl")
include("integrators/dvi/dvi_cache.jl")
include("integrators/dvi/dvi_euler.jl")
include("integrators/dvi/dvi_midpoint.jl")
include("integrators/dvi/dvi_trapezoidal.jl")
include("integrators/dvi/dvrk.jl")


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
