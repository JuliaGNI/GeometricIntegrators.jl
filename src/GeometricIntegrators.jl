module GeometricIntegrators

using Reexport

@reexport using GeometricBase
@reexport using GeometricEquations
@reexport using GeometricIntegratorsBase
@reexport using GeometricSolutions

# import GeometricBase.Config

include("Discontinuities.jl")
@reexport using .Discontinuities
include("Integrators.jl")
@reexport using .Integrators
include("Simulations.jl")
@reexport using .Simulations

include("SPARK.jl")

# function __init__()
#     add_config(:verbosity, 1)
# end

end
