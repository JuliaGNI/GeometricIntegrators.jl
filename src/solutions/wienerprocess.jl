
using Random: GLOBAL_RNG

abstract type SemiMartingale{dType, tType, N} end


"""
Type for holding the increments of a Wiener process

### Fields:
* `nd`: dimension of the Wiener process
* `nt`: number of increments in the DataSeries
* `ns`: number of sample paths of the Wiener process
* `Δt`: time increment of the TimeSeries
* `ΔW`: variable storing the increments of the Wiener process over Δt, or the discrete random variable \\hat{I}
* `ΔZ`: variable holding the time integral of the Wiener process \\int_{tk}^{tk+1} (W(t)-W(tk))dt, or the discrete random variable \\tilde{I}

### Parameters:
* `dType`: type of the elements of the increments of the Wiener process
* `tType`: type of the time steps
* `N`:     the number of dimensions of the arrays holding data in ΔW and ΔZ
* `CONV`:  mode of convergence: `:strong` or `:weak`
"""
struct WienerProcess{dType, tType, N, CONV} <: SemiMartingale{dType, tType, N}
    nd::Int
    nt::Int
    ns::Int
    Δt::tType
    ΔW::Array{dType,N}
    ΔZ::Array{dType,N}

    function WienerProcess{dType, tType, N, CONV}(nd, nt, ns, Δt,
                    ΔW = (ns == 1 ? zeros(dType, nd, nt) : zeros(dType, nd, nt, ns)),
                    ΔZ = (ns == 1 ? zeros(dType, nd, nt) : zeros(dType, nd, nt, ns))) where {dType, tType, N, CONV}

        @assert CONV==:strong || CONV==:weak || CONV==:null
        @assert nd > 0
        @assert nt ≥ 1
        @assert ns > 0
        @assert N ∈ (2,3)

        if CONV != :null
            @assert nd == size(ΔW, 1) == size(ΔZ, 1)
            @assert nt == size(ΔW, 2) == size(ΔZ, 2)
            @assert ns == size(ΔW, 3) == size(ΔZ, 3)
        end

        new(nd, nt, ns, Δt, ΔW, ΔZ)
    end

end


function WienerProcess(dType, nd, nt, ns, Δt::tType, conv=DEFAULT_SCONV; rng=GLOBAL_RNG) where {tType <: Number}
    wp = WienerProcess{dType, tType, ns==1 ? 2 : 3, conv}(nd, nt, ns, Δt)
    generate_wienerprocess!(wp, rng)
    return wp
end


function WienerProcess(Δt::tType, dW::Array{dType,1}, dZ::Array{dType,1}, conv=DEFAULT_SCONV) where {dType, tType}
    @assert size(dW,1) == size(dZ,1)

    nd = 1
    nt = size(dW,1)
    ns = 1

    return WienerProcess{dType, tType, 2, conv}(nd, nt, ns, Δt, reshape(dW, (1,nt)), reshape(dZ, (1,nt)))
end


function WienerProcess(Δt::tType, dW::Array{dType,N}, dZ::Array{dType,N}, conv=DEFAULT_SCONV) where {dType, tType, N}
    @assert size(dW) == size(dZ)

    nd = size(dW,1)
    nt = size(dW,2)
    ns = size(dW,3)

    return WienerProcess{dType, tType, N, conv}(nd, nt, ns, Δt, dW, dZ)
end


Base.:(==)(wp1::WienerProcess{DT1, TT1, N1, C1}, wp2::WienerProcess{DT2, TT2, N2, C2}) where {DT1, TT1, N1, C1, DT2, TT2, N2, C2} = (
                                DT1 == DT2
                             && TT1 == TT2
                             && N1  == N2
                             && C1  == C2
                             && wp1.nd == wp2.nd
                             && wp1.nt == wp2.nt
                             && wp1.ns == wp2.ns
                             && wp1.Δt == wp2.Δt
                             && wp1.ΔW == wp2.ΔW
                             && wp1.ΔZ == wp2.ΔZ)


function set_chi_and_eta!(W::WienerProcess{dType, tType, N, :strong}, chi::AbstractArray{dType}, eta::AbstractArray{dType}) where {dType, tType, N}
    @. W.ΔW = chi * √W.Δt
    @. W.ΔZ = W.Δt^(3/2)/2 * (chi+eta/√3)
    return nothing
end

function set_chi_and_eta!(W::WienerProcess{dType, tType, N, :weak}, chi::AbstractArray{dType}, eta::AbstractArray{dType}) where {dType, tType, N}
    W.ΔW .= 0
    W.ΔZ .= sqrt(W.Δt)

    # Generating the random variable \\hat I
    # P(\\hat I=+-sqrt(3*Δt)) = 1/6
    # P(\\hat I=0) = 2/3
    indx        = findall(x->x<1. / 6., chi)
    W.ΔW[indx] .= -sqrt(3*W.Δt)
    indx        = findall(x->x>5. / 6., chi)
    W.ΔW[indx] .= sqrt(3*W.Δt)

    # Generating the random variable \\tilde I
    # P(\\tilde I=+-sqrt(Δt)) = 1/2
    indx        = findall(x->x<0.5, eta)
    W.ΔZ[indx] .= -sqrt(W.Δt)

    return nothing
end


# Generates a new series of increments for the strong Wiener process W
function generate_wienerprocess!(W::WienerProcess{DT, TT, 2, :strong}, rng=GLOBAL_RNG) where {DT,TT}
    chi = randn(rng, DT, W.nd, W.nt)
    eta = randn(rng, DT, W.nd, W.nt)
    set_chi_and_eta!(W, chi, eta)
end

# Generates a new series of increments for the strong Wiener process W
function generate_wienerprocess!(W::WienerProcess{DT, TT, 3, :strong}, rng=GLOBAL_RNG) where {DT,TT}
    chi = randn(rng, DT, W.nd, W.nt, W.ns)
    eta = randn(rng, DT, W.nd, W.nt, W.ns)
    set_chi_and_eta!(W, chi, eta)
end

# Generates a new series of increments for the weak Wiener process W
function generate_wienerprocess!(W::WienerProcess{DT, TT, 2, :weak}, rng=GLOBAL_RNG) where {DT,TT}
    chi = rand(rng, DT, W.nd, W.nt)
    eta = rand(rng, DT, W.nd, W.nt)
    set_chi_and_eta!(W, chi, eta)
end

# Generates a new series of increments for the weak Wiener process W
function generate_wienerprocess!(W::WienerProcess{DT, TT, 3, :weak}, rng=GLOBAL_RNG) where {DT,TT}
    chi = rand(rng, DT, W.nd, W.nt, W.ns)
    eta = rand(rng, DT, W.nd, W.nt, W.ns)
    set_chi_and_eta!(W, chi, eta)
end


Base.ndims(ds::WienerProcess{DT,TT,N}) where {DT,TT,N} = N
