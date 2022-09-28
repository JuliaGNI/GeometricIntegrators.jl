module Solutions

using Reexport

using HDF5
using OffsetArrays
using SharedArrays

using GeometricBase
using GeometricBase.Config
using GeometricBase.Utils
using GeometricEquations

@reexport using GeometricSolutions

export DEFAULT_NSAVE, DEFAULT_NWRITE

const DEFAULT_NSAVE = 1
const DEFAULT_NWRITE = 0


export DataSeries, TimeSeries, Solution, AbstractSolution, DeterministicSolution
export current, previous

export SolutionODE, SolutionPODE
export SolutionDAE, SolutionPDAE

export SolutionStep,
       SolutionStepODE, SolutionStepPODE,
       SolutionStepDAE, SolutionStepPDAE

export update!, cut_periodic_solution!

include("solutions/solution_step.jl")
include("solutions/solution.jl")
include("solutions/solution_step_ode.jl")
include("solutions/solution_step_pode.jl")
include("solutions/solution_step_dae.jl")
include("solutions/solution_step_pdae.jl")
include("solutions/solution_step_constructors.jl")

end
