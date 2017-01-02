__precompile__()

module Problems

    using ..Equations

    export pendulum_ode, pendulum_pode, pendulum_iode, pendulum_idae

    include("problems/pendulum.jl")

end
