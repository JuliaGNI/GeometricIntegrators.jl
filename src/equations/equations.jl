
abstract Equation{dType <: Number, tType <: Real}

include("ode.jl")
include("iode.jl")
include("pode.jl")
include("dae.jl")
include("pdae.jl")

# TODO Add functions and vectors for invariants (momentum maps, energy, ...).
