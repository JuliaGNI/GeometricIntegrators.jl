__precompile__()

module Solutions

    using DoubleDouble

    using ..Equations

    export DataSeries, get_data!, set_data!, reset!
    export PDataSeries, SDataSeries

    include("solutions/dataseries.jl")

    export TimeSeries, compute_timeseries!

    include("solutions/timeseries.jl")

    export Solution, SolutionODE, SolutionPODE, SolutionDAE, SolutionPDAE,
           copy_solution!, reset!,
           get_initial_conditions!, set_initial_conditions!,
           createHDF5, writeSolutionToHDF5
    export PSolutionPDAE, SSolutionPDAE

    include("solutions/solutions.jl")
    include("solutions/solutions_ode.jl")
    include("solutions/solutions_pode.jl")
    include("solutions/solutions_dae.jl")
    include("solutions/solutions_pdae.jl")

end
