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

export AtomicSolution,
    AtomicSolutionODE, AtomicSolutionPODE,
    AtomicSolutionDAE, AtomicSolutionPDAE
export update!, cut_periodic_solution!

include("solutions/atomic_solution.jl")
include("solutions/atomic_solution_ode.jl")
include("solutions/atomic_solution_pode.jl")
include("solutions/atomic_solution_dae.jl")
include("solutions/atomic_solution_pdae.jl")



export SolutionODE, SolutionPODE
export SolutionDAE, SolutionPDAE

include("solutions/solution.jl")

include("solutions/atomic_solution_constructors.jl")

end
