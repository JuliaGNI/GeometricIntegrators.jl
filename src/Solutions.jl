module Solutions

    using Base: TwicePrecision

    using HDF5
    using OffsetArrays
    using SharedArrays

    using GeometricBase
    using GeometricBase.Config
    using GeometricBase.Utils
    using GeometricEquations

    export DEFAULT_NSAVE, DEFAULT_NWRITE

    const DEFAULT_NSAVE = 1
    const DEFAULT_NWRITE = 0

    export SolutionVector

    SolutionVector{DT} = Union{Vector{DT}, Vector{TwicePrecision{DT}}}


    export get_data!, set_data!
    export DataSeries, PDataSeries, SDataSeries

    include("solutions/dataseries.jl")

    export TimeSeries, compute_timeseries!

    include("solutions/timeseries.jl")

    export Solution, ParallelSolution, DeterministicSolution
    export nsave, nsamples, timesteps, counter, offset, lastentry, hdf5
    export get_solution, get_solution!, set_solution!

    include("solutions/solution.jl")

    export AtomicSolution,
           AtomicSolutionODE, AtomicSolutionPODE,
           AtomicSolutionDAE, AtomicSolutionPDAE
    export update!, cut_periodic_solution!

    include("solutions/atomic_solution.jl")
    include("solutions/atomic_solution_ode.jl")
    include("solutions/atomic_solution_pode.jl")
    include("solutions/atomic_solution_dae.jl")
    include("solutions/atomic_solution_pdae.jl")

    export SolutionODE, SSolutionODE, PSolutionODE, SolutionPODE, SSolutionPODE, PSolutionPODE
    export SolutionDAE, SSolutionDAE, PSolutionDAE, SolutionPDAE, SSolutionPDAE, PSolutionPDAE
    export get_initial_conditions, get_initial_conditions!, set_initial_conditions!,
           create_hdf5, create_hdf5!

    include("solutions/solution_ode.jl")
    include("solutions/solution_pode.jl")
    include("solutions/solution_dae.jl")
    include("solutions/solution_pdae.jl")

    include("solutions/solution_hdf5.jl")

    include("solutions/atomic_solution_constructors.jl")
    
end
