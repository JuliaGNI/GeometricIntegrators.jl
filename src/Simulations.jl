module Simulations

    using ProgressMeter
    using RungeKutta
    using RungeKutta: AbstractTableau

    using GeometricBase
    using GeometricEquations

    using ..Integrators
    using ..Solutions

    import ..Solutions: DEFAULT_NSAVE, DEFAULT_NWRITE

    export Simulation, ParallelSimulation
    export run!

    include("simulations/simulation.jl")
    include("simulations/parallel_simulation.jl")

end
