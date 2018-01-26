__precompile__()

module Equations

    export Equation, ODE, IODE, PODE, SODE, VODE, DAE, HDAE, IDAE, PDAE, SDE

    include("equations/equations.jl")

    include("equations/ode.jl")
    include("equations/iode.jl")
    include("equations/pode.jl")
    include("equations/sode.jl")
    include("equations/vode.jl")
    include("equations/dae.jl")
    include("equations/hdae.jl")
    include("equations/idae.jl")
    include("equations/pdae.jl")
    include("equations/sde.jl")

end
