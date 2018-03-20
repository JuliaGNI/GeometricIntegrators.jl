__precompile__()

module Solutions

    using DoubleDouble
    using HDF5

    using ..CommonFunctions
    using ..Equations

    export DataSeries, get_data!, set_data!, reset!
    export PDataSeries, SDataSeries

    include("solutions/dataseries.jl")

    export TimeSeries, compute_timeseries!

    include("solutions/timeseries.jl")

    export StochasticDataSeries, get_data!, set_data!, reset!
    export SStochasticDataSeries

    include("solutions/stochasticdataseries.jl")

    export Solution, SolutionODE, SolutionPODE, SolutionDAE, SolutionPDAE, SolutionSDE,
           copy_solution!, reset!,
           get_initial_conditions!, set_initial_conditions!,
           create_hdf5, write_to_hdf5
    export PSolutionPDAE, SSolutionPDAE

    include("solutions/solutions.jl")
    include("solutions/solutions_ode.jl")
    include("solutions/solutions_pode.jl")
    include("solutions/solutions_dae.jl")
    include("solutions/solutions_pdae.jl")
    include("solutions/solutions_sde.jl")

end
