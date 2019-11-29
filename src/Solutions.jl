module Solutions

    using HDF5
    using SharedArrays

    using Base: TwicePrecision

    using ..CommonFunctions
    using ..Config
    using ..Equations


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

    export StochasticDataSeries, SStochasticDataSeries

    include("solutions/stochasticdataseries.jl")

    export SemiMartingale
    export WienerProcess, generate_wienerprocess!

    include("solutions/wienerprocess.jl")


    export AtomisticSolution, AtomisticSolutionODE, AtomisticSolutionPODE,
           AtomisticSolutionDAE, AtomisticSolutionPDAE

    include("solutions/atomistic_solution.jl")
    include("solutions/atomistic_solution_ode.jl")
    include("solutions/atomistic_solution_pode.jl")
    include("solutions/atomistic_solution_dae.jl")
    include("solutions/atomistic_solution_pdae.jl")


    export Solution, StochasticSolution, SolutionODE, SolutionPODE, SolutionDAE, SolutionPDAE, SolutionSDE, SolutionPSDE,
           get_initial_conditions, get_initial_conditions!, set_initial_conditions!,
           create_hdf5
    export PSolutionPDAE, SSolutionPDAE

    include("solutions/solution.jl")

    export AtomisticSolution, AtomisticSolutionODE, AtomisticSolutionPODE,
           AtomisticSolutionDAE, AtomisticSolutionPDAE
    export update!, cut_periodic_solution!

    include("solutions/atomistic_solution.jl")
    include("solutions/atomistic_solution_ode.jl")
    include("solutions/atomistic_solution_pode.jl")
    include("solutions/atomistic_solution_dae.jl")
    include("solutions/atomistic_solution_pdae.jl")

    export SolutionODE, SolutionPODE, SolutionDAE, SolutionPDAE, SolutionSDE, SolutionPSDE
    export PSolutionPDAE, SSolutionPDAE
    export get_initial_conditions, get_initial_conditions!, set_initial_conditions!,
           create_hdf5

    include("solutions/solution_ode.jl")
    include("solutions/solution_pode.jl")
    include("solutions/solution_dae.jl")
    include("solutions/solution_pdae.jl")
    include("solutions/solution_sde.jl")
    include("solutions/solution_psde.jl")

end
