
abstract type SemiMartingale{dType, tType, N} end


#Type for holding the increments of a Wiener process:
# conv - mode of convergence: :strong or :weak
# nd   - dimension of the Wiener process
# nt   - number of increments in the DataSeries
# ns   - number of sample paths of the Wiener process
# Δt   - time increment of the TimeSeries
# ΔW   - variable storing the increments of the Wiener process over Δt, or the discrete random variable \\hat{I}
# ΔZ   - variable holding the time integral of the Wiener process \\int_{tk}^{tk+1} (W(t)-W(tk))dt, or the discrete random variable \\tilde{I}
# dType- type of the elements of the increments of the Wiener process
# tType- type of the time steps
# N    - the number of dimensions of the arrays holding data in ΔW and ΔZ
struct WienerProcess{dType, tType, N, conv} <: SemiMartingale{dType, tType, N}
    nd::Int
    nt::Int
    ns::Int
    Δt::tType
    ΔW::Array{dType,N}
    ΔZ::Array{dType,N}

    function WienerProcess{dType, tType, N, conv}(nd, nt, ns, Δt, rng=MersenneTwister()) where {dType <: Number, tType <: Number, N, conv}

        @assert conv==:strong || conv==:weak
        @assert nd > 0
        @assert nt ≥ 1
        @assert ns > 0
        @assert N ∈ (2,3)

        if conv==:strong
            if N == 2
                chi = randn(rng, dType, nd, nt)
                eta = randn(rng, dType, nd, nt)
            elseif N == 3
                chi = randn(rng, dType, nd, nt, ns)
                eta = randn(rng, dType, nd, nt, ns)
            end

            ΔW = chi*√Δt
            ΔZ = Δt^(3/2)/2 * (chi+eta/√3)
        else
            if N == 2
                chi = rand(rng, dType, nd, nt)
                eta = rand(rng, dType, nd, nt)
                ΔW  = zeros(dType, nd, nt)
                ΔZ  = sqrt(Δt)*ones(dType, nd, nt)
            elseif N == 3
                chi = rand(rng, dType, nd, nt, ns)
                eta = rand(rng, dType, nd, nt, ns)
                ΔW  = zeros(dType, nd, nt, ns)
                ΔZ  = sqrt(Δt)*ones(dType, nd, nt, ns)
            end

            # Generating the random variable \hat I
            # P(\\hat I=+-sqrt(3*Δt)) = 1/6
            # P(\\hat I=0) = 2/3
            indx      = findall(x->x<1. / 6., chi)
            ΔW[indx] .= -sqrt(3*Δt)
            indx      = findall(x->x>5. / 6., chi)
            ΔW[indx] .= sqrt(3*Δt)

            # Generating the random variable \tilde I
            # P(\tilde I=+-sqrt(Δt)) = 1/2
            indx      = findall(x->x<0.5, eta)
            ΔZ[indx] .= -sqrt(Δt)
        end

        new(nd, nt, ns, Δt, ΔW, ΔZ)
    end


    function WienerProcess{dType, tType, N, conv}(nd, nt, ns, Δt, dW, dZ) where {dType, tType, N, conv}

        @assert conv==:strong || conv==:weak

        new(conv, nd, nt, ns, Δt, SStochasticDataSeries(dW), SStochasticDataSeries(dZ))
    end

end


function WienerProcess(dType, nd::Int, nt::Int, ns::Int, Δt::tType, conv=:strong) where {tType}
    ns==1 ? N=2 : N=3
    WienerProcess{dType, tType, N, conv}(nd, nt, ns, Δt)
end


function WienerProcess(Δt::tType, dW::Array{dType,1}, dZ::Array{dType,1}, conv=:strong) where {dType, tType}

    @assert size(dW,1)==size(dZ,1)

    nd = 1
    nt = size(dW,1)
    ns = 1

    return WienerProcess{dType, tType, 2, conv}(nd, nt, ns, Δt, reshape(dW, (1,nt)), reshape(dZ, (1,nt)))
end


function WienerProcess(Δt::tType, dW::Array{dType,2}, dZ::Array{dType,2}, conv=:strong) where {dType, tType}

    @assert size(dW,1)==size(dZ,1)
    @assert size(dW,2)==size(dZ,2)

    nd = size(dW,1)
    nt = size(dW,2)
    ns = 1

    return WienerProcess{dType, tType, 2, conv}(nd, nt, ns, Δt, dW, dZ)
end


function WienerProcess(Δt::tType, dW::Array{dType,3}, dZ::Array{dType,3}, conv=:strong) where {dType, tType}

    @assert size(dW,1)==size(dZ,1)
    @assert size(dW,2)==size(dZ,2)
    @assert size(dW,3)==size(dZ,3)

    nd = size(dW,1)
    nt = size(dW,2)
    ns = size(dW,3)

    return WienerProcess{dType, tType, 3, conv}(nd, nt, ns, Δt, dW, dZ)
end


# Generates a new series of increments for the strong Wiener process W
function generate_wienerprocess!(W::WienerProcess{dType, tType, N, :strong}, rng=MersenneTwister()) where {dType, tType, N}
    @assert N ∈ (2,3)

    if N == 2
        chi = randn(rng, dType, W.nd, W.nt)
        eta = randn(rng, dType, W.nd, W.nt)
    elseif N == 3
        chi = randn(rng, dType, W.nd, W.nt, W.ns)
        eta = randn(rng, dType, W.nd, W.nt, W.ns)
    end

    dW = chi*√W.Δt
    dZ = W.Δt^(3/2)/2 * (chi+eta/√3)

    W.ΔW.d .= dW
    W.ΔZ.d .= dZ

    return nothing
end


# Generates a new series of increments for the weak Wiener process W
function generate_wienerprocess!(W::WienerProcess{dType, tType, N, :weak}, rng=MersenneTwister()) where {dType, tType, N}
    @assert N ∈ (2,3)

    if N == 2
        chi = rand(rng, dType, W.nd, W.nt)
        eta = rand(rng, dType, W.nd, W.nt)
        dW  = zeros(dType, W.nd, W.nt)
        dZ  = sqrt(W.Δt)*ones(dType, W.nd, W.nt)
    elseif N == 3
        chi = rand(rng, dType, W.nd, W.nt, W.ns)
        eta = rand(rng, dType, W.nd, W.nt, W.ns)
        dW  = zeros(dType, W.nd, W.nt, W.ns)
        dZ  = sqrt(W.Δt)*ones(dType, W.nd, W.nt, W.ns)
    end

    # Generating the random variable \\hat I
    # P(\\hat I=+-sqrt(3*Δt)) = 1/6
    # P(\\hat I=0) = 2/3
    indx      = findall(x->x<1. / 6., chi)
    dW[indx] .= -sqrt(3*W.Δt)
    indx      = findall(x->x>5. / 6., chi)
    dW[indx] .= sqrt(3*W.Δt)

    # Generating the random variable \\tilde I
    # P(\\tilde I=+-sqrt(Δt)) = 1/2
    indx      = findall(x->x<0.5, eta)
    dZ[indx] .= -sqrt(W.Δt)

    W.ΔW.d .= dW
    W.ΔZ.d .= dZ

    return nothing
end


Base.ndims(ds::WienerProcess{DT,TT,N}) where {DT,TT,N} = N
