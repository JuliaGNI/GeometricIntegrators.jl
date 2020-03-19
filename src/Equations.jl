module Equations

    using ..CommonFunctions

    export Equation
    export AbstractEquationODE, AbstractEquationPODE,
           AbstractEquationDAE, AbstractEquationPDAE,
           AbstractEquationSDE, AbstractEquationPSDE

    export ODE, IODE, PODE, HODE, VODE, SODE
    export DAE, IDAE, PDAE, HDAE, VDAE, SPDAE
    export SDE, PSDE, SPSDE

    export get_function_tuple

    include("equations/equations.jl")

    include("equations/ode.jl")
    include("equations/hode.jl")
    include("equations/iode.jl")
    include("equations/pode.jl")
    include("equations/vode.jl")

    include("equations/dae.jl")
    include("equations/hdae.jl")
    include("equations/idae.jl")
    include("equations/pdae.jl")
    include("equations/vdae.jl")

    include("equations/sode.jl")
    include("equations/spdae.jl")

    include("equations/sde.jl")
    include("equations/psde.jl")
    include("equations/spsde.jl")

end
