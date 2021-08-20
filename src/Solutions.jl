module Solutions

    using HDF5
    using OffsetArrays
    using SharedArrays

    using GeometricBase
    using GeometricBase.Config
    using GeometricBase.Utils
    using GeometricEquations

    using GeometricBase: fromarray

    export DEFAULT_NSAVE, DEFAULT_NWRITE

    const DEFAULT_NSAVE = 1
    const DEFAULT_NWRITE = 0


    export DataSeries

    include("solutions/dataseries.jl")

    export Solution, ParallelSolution, DeterministicSolution
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
    export get_initial_conditions, get_initial_conditions!, set_initial_conditions!
    export create_hdf5, create_hdf5!, write_to_hdf5, hdf5

    include("solutions/solution_ode.jl")
    include("solutions/solution_pode.jl")
    include("solutions/solution_dae.jl")
    include("solutions/solution_pdae.jl")

    include("solutions/solution_hdf5.jl")

    include("solutions/atomic_solution_constructors.jl")
    
end
