
using Plots

@userplot PlotOscillator
@recipe function f(p::PlotOscillator)
    if length(p.args) != 2 || !(typeof(p.args[1]) <: AbstractVector) || !(typeof(p.args[2]) <: AbstractMatrix)
        error("Oscillator plots should be called with two arrays. Got: $(typeof(p.args))")
    end
    local t = p.args[1]
    local q = p.args[2]
    local h = zero(t)

    for i in eachindex(t)
        h[i] = hamiltonian(t[i], q[:,i], k)
    end

    size   := (800,600)
    layout := @layout [xPlot
                       yPlot
                       hPlot]
    legend := :none

    @series begin
        subplot := 1
        xlabel := "t"
        ylabel := "x"
        t, q[1,:]
    end

    @series begin
        subplot := 2
        xlabel := "t"
        ylabel := "y"
        t, q[2,:]
    end

    @series begin
        subplot := 3
        xlabel := "t"
        ylabel := "h"
        t, h
    end
end


@userplot PlotOscillatorOverview
@recipe function f(p::PlotOscillatorOverview)
    if length(p.args) != 1 || !(typeof(p.args[1]) <: AbstractArray)
        error("Oscillator plots should be called with an array. Got: $(typeof(p.args))")
    end

    local q = p.args[1]

    size   := (800, div(size(q,2), 2)*300)

    layout := (div(size(q,2), 2), 2)
    legend := :none

    for k in axes(q,2)
        @series begin
            subplot := k
            xlabel  := "x"
            ylabel  := "y"
            q[1,k,:], q[2,k,:]
        end
    end
end
