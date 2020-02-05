module Solutions

    using HDF5
    using SharedArrays

    using Base: TwicePrecision
    using Random

    using ..CommonFunctions
    using ..Config
    using ..Equations
    using ..Utils


    export DEFAULT_NSAVE, DEFAULT_NWRITE

    const DEFAULT_NSAVE = 1
    const DEFAULT_NWRITE = 0
    const DEFAULT_SCONV = :strong

    export SolutionVector

    SolutionVector{DT} = Union{Vector{DT}, Vector{TwicePrecision{DT}}}


    export get_data!, set_data!
    export DataSeries, PDataSeries, SDataSeries

    include("solutions/dataseries.jl")

    export TimeSeries, compute_timeseries!

    include("solutions/timeseries.jl")

    export SemiMartingale
    export WienerProcess, generate_wienerprocess!

    include("solutions/wienerprocess.jl")

    export Solution, StochasticSolution
    export nsave, ntime, timesteps, offset, conv, hdf5
    export get_solution, get_solution!, set_solution!

    include("solutions/solution.jl")

    export AtomicSolution,
           AtomicSolutionODE, AtomicSolutionPODE,
           AtomicSolutionDAE, AtomicSolutionPDAE,
           AtomicSolutionSDE, AtomicSolutionPSDE
    export update!, cut_periodic_solution!,
           get_increments, get_increments!, set_increments!

    include("solutions/atomic_solution.jl")
    include("solutions/atomic_solution_ode.jl")
    include("solutions/atomic_solution_pode.jl")
    include("solutions/atomic_solution_dae.jl")
    include("solutions/atomic_solution_pdae.jl")
    include("solutions/atomic_solution_sde.jl")
    include("solutions/atomic_solution_psde.jl")

    export SolutionODE, SSolutionODE, PSolutionODE, SolutionPODE, SSolutionPODE, PSolutionPODE
    export SolutionDAE, SSolutionDAE, PSolutionDAE, SolutionPDAE, SSolutionPDAE, PSolutionPDAE
    export SolutionSDE, SSolutionSDE, PSolutionSDE, SolutionPSDE, SSolutionPSDE, PSolutionPSDE
    export get_initial_conditions, get_initial_conditions!, set_initial_conditions!,
           create_hdf5, create_hdf5!

    include("solutions/solution_ode.jl")
    include("solutions/solution_pode.jl")
    include("solutions/solution_dae.jl")
    include("solutions/solution_pdae.jl")
    include("solutions/solution_sde.jl")
    include("solutions/solution_psde.jl")

    include("solutions/solution_hdf5.jl")

end
