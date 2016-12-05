
abstract Equation{dType <: Number, tType <: Real}

function do_nothing(x...) end

include("ode.jl")
include("pode.jl")
include("sode.jl")
include("dae.jl")
include("pdae.jl")

# TODO Add functions and vectors for invariants (momentum maps, energy, ...).
