module TestProblems

    using Reexport

    include("problems/oscillator.jl")
    include("problems/kubo_oscillator.jl")
    include("problems/lotka_volterra.jl")

    @reexport using .Oscillator
    @reexport using .KuboOscillator
    @reexport using .LotkaVolterra

end
