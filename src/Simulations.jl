module Simulations

    using ..Equations
    using ..Integrators
    using ..Solutions

    export Simulation, run!

    include("simulations/simulation.jl")

end
