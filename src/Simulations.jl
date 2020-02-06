module Simulations

    using ProgressMeter

    using ..CommonFunctions
    using ..Equations
    using ..Integrators
    using ..Solutions

    export Simulation, ParallelSimulation
    export run!

    include("simulations/simulation.jl")
    include("simulations/parallel_simulation.jl")

end
