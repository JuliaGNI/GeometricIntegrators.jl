__precompile__()

module Problems

    using ..Equations

    export exponential_growth_ode

    include("problems/exponential_growth.jl")

    export lotka_volterra_2d_ode, lotka_volterra_2d_iode, lotka_volterra_2d_idae

    include("problems/lotka_volterra_2d.jl")

    export pendulum_ode, pendulum_pode, pendulum_iode, pendulum_idae

    include("problems/pendulum.jl")

end
