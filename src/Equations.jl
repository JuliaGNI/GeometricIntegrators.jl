__precompile__()

module Equations

    export Equation, ODE, IODE, PODE, DAE, IDAE, PDAE

    include("equations/equations.jl")

    include("equations/ode.jl")
    include("equations/iode.jl")
    include("equations/pode.jl")
    include("equations/vode.jl")
    include("equations/dae.jl")
    include("equations/idae.jl")
    include("equations/pdae.jl")

end
