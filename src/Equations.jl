module Equations

    using ..Common

    export Equation
    export AbstractEquationODE, AbstractEquationPODE,
           AbstractEquationDAE, AbstractEquationPDAE,
           AbstractEquationSDE, AbstractEquationPSDE

    export ODE, IODE, PODE, HODE, LODE, SODE
    export DAE, IDAE, PDAE, HDAE, LDAE, SPDAE
    export SDE, PSDE, SPSDE

    export initial_conditions
    export get_function_tuple, get_solution_tuple, hassolution, hasvectorfield

    include("equations/equations.jl")

    include("equations/ode.jl")
    include("equations/hode.jl")
    include("equations/iode.jl")
    include("equations/lode.jl")
    include("equations/pode.jl")

    include("equations/dae.jl")
    include("equations/hdae.jl")
    include("equations/idae.jl")
    include("equations/ldae.jl")
    include("equations/pdae.jl")

    include("equations/sode.jl")
    include("equations/spdae.jl")

    include("equations/sde.jl")
    include("equations/psde.jl")
    include("equations/spsde.jl")

end
