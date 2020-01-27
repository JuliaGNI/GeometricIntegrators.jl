module TestProblems

    using Reexport

    include("problems/harmonic_oscillator.jl")
    include("problems/kubo_oscillator.jl")
    include("problems/lotka_volterra_2d.jl")

    @reexport using .HarmonicOscillator
    @reexport using .KuboOscillator
    @reexport using .LotkaVolterra2d

end
