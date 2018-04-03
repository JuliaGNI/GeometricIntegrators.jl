
abstract type SemiMartingale{dType, tType, N} end


#Type for holding the increments of a Wiener process:
# nd   - dimension of the Wiener process
# nt   - number of increments in the DataSeries
# ns   - number of sample paths of the Wiener process
# Δt   - time increment of the TimeSeries
# ΔW   - variable storing the increments of the Wiener process over Δt
# ΔZ   - variable holding the time integral of the Wiener process \int_{tk}^{tk+1} (W(t)-W(tk))dt
# dType- type of the elements of the increments of the Wiener process
# tType- type of the time steps
# N    - the number of dimensions of the arrays holding data in ΔW and ΔZ
struct WienerProcess{dType, tType, N} <: SemiMartingale{dType, tType, N}
    nd::Int
    nt::Int
    ns::Int
    Δt::tType
    ΔW::SStochasticDataSeries{dType,N}
    ΔZ::SStochasticDataSeries{dType,N}

    function WienerProcess{dType, tType, N}(nd, nt, ns, Δt) where {dType, tType, N}

        @assert dType <: Number
        @assert tType <: Number
        @assert nd > 0
        @assert nt ≥ 1
        @assert ns > 0
        @assert N ∈ (1,2,3)

        if N == 1
            chi = randn(nt)
            eta = randn(nt)
        elseif N == 2
            chi = randn(nd, nt)
            eta = randn(nd, nt)
        elseif N == 3
            chi = randn(nd, nt, ns)
            eta = randn(nd, nt, ns)
        end

        dW = chi*√Δt
        dZ = Δt^(3/2)/2 * (chi+eta/√3)

        # Creating ΔW, ΔZ
        # nt-1, because SStochasticDataSeries adds one time step
        # ni=1, because we cosider a single IC for the Wiener process
        ΔW = SStochasticDataSeries{dType,N}(nd, nt-1, ns, 1, dW)
        ΔZ = SStochasticDataSeries{dType,N}(nd, nt-1, ns, 1, dZ)

        new(nd, nt, ns, Δt, ΔW, ΔZ)
    end


    function WienerProcess{dType, tType, N}(nd, nt, ns, Δt, dW, dZ) where {dType, tType, N}
        new(nd, nt, ns, Δt, SStochasticDataSeries(dW), SStochasticDataSeries(dZ))
    end

end


function WienerProcess(dType, nd::Int, nt::Int, ns::Int, Δt::tType) where {tType}

    if nd==ns==1
        N=1
    elseif ns==1
        N=2
    else
        N=3
    end

    WienerProcess{dType, tType, N}(nd, nt, ns, Δt)
end


function WienerProcess(Δt::tType, dW::Array{dType,1}, dZ::Array{dType,1}) where {dType, tType}

    @assert size(dW,1)==size(dZ,1)

    nd = 1
    nt = size(dW,1)
    ns = 1

    return WienerProcess{dType, tType, 1}(nd, nt, ns, Δt, dW, dZ)
end


function WienerProcess(Δt::tType, dW::Array{dType,2}, dZ::Array{dType,2}) where {dType, tType}

    @assert size(dW,1)==size(dZ,1)
    @assert size(dW,2)==size(dZ,2)

    nd = size(dW,1)
    nt = size(dW,2)
    ns = 1

    return WienerProcess{dType, tType, 2}(nd, nt, ns, Δt, dW, dZ)
end


function WienerProcess(Δt::tType, dW::Array{dType,3}, dZ::Array{dType,3}) where {dType, tType}

    @assert size(dW,1)==size(dZ,1)
    @assert size(dW,2)==size(dZ,2)
    @assert size(dW,3)==size(dZ,3)

    nd = size(dW,1)
    nt = size(dW,2)
    ns = size(dW,3)

    return WienerProcess{dType, tType, 3}(nd, nt, ns, Δt, dW, dZ)
end
