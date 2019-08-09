__precompile__()

module Solutions

    using HDF5
    using SharedArrays

    using Base: TwicePrecision

    using ..CommonFunctions
    using ..Equations


    export DEFAULT_NSAVE, DEFAULT_NWRITE

    const DEFAULT_NSAVE = 1
    const DEFAULT_NWRITE = 0


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

    export Solution, StochasticSolution, SolutionODE, SolutionPODE, SolutionDAE, SolutionPDAE, SolutionSDE, SolutionPSDE,
           get_initial_conditions, get_initial_conditions!, set_initial_conditions!,
           create_hdf5
    export PSolutionPDAE, SSolutionPDAE

    include("solutions/solutions.jl")
    include("solutions/solutions_ode.jl")
    include("solutions/solutions_pode.jl")
    include("solutions/solutions_dae.jl")
    include("solutions/solutions_pdae.jl")
    include("solutions/solutions_sde.jl")
    include("solutions/solutions_psde.jl")

end
