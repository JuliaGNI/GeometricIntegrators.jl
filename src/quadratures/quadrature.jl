
import FastGaussQuadrature
import FastTransforms

using DecFP: Dec128

using ..CommonFunctions


struct Quadrature{T,N}
    order::Int
    nodes::Vector{T}
    weights::Vector{T}
end

function Quadrature(order, nodes::Vector{T}, weights::Vector{T}) where {T}
    @assert length(nodes) == length(weights)
    Quadrature{T, length(nodes)}(order, nodes, weights)
end

CommonFunctions.nnodes(Q::Quadrature{T,N}) where {T,N} = N
CommonFunctions.nodes(Q::Quadrature) = Q.nodes
CommonFunctions.order(Q::Quadrature) = Q.order
weights(Q::Quadrature) = Q.weights


"Scale nodes and weights from the interval [-1,+1] to the interval [0,1]"
function shift!(b,c)
    b .= b ./ 2
    c .= (c .+ 1) ./ 2
end


function GaussLegendreQuadrature(s, T=Float64)
    c, b = FastGaussQuadrature.gausslegendre(s)
    shift!(b,c)
    Quadrature(2s, c, b)
end


function LobattoLegendreQuadrature(s, T=Float64)
    c, b = FastGaussQuadrature.gausslobatto(s)
    shift!(b,c)
    Quadrature(2s-2, c, b)
end


function get_nodes_chebychev_gauss(s, k=1, T=Float64)
    x = zeros(T, s)
    if k == 1
        for i in eachindex(x)
            x[i] = @dec128 ( sin( π*(s-2i+1) / (2s) ) + 1 ) / 2
        end
        return x[end:-1:1]
    elseif k == 2
        for i in eachindex(x)
            x[i] = @dec128 ( cos( π*i / (s-1) ) + 1 ) / 2
        end
        return x
    end
end

function GaussChebyshevQuadrature(s, k=1, T=Float64)
    c = get_nodes_chebychev_gauss(s,k,T)
    b = zeros(c)

    tj::Dec128 = 0
    th::Dec128 = 0

    if k == 1
        for i in eachindex(b)
            tj = 0
            th = @dec128 π * (2i-1) / (2s)
            for j in 1:div(s,2)
                tj += @dec128 cos(2j*th) / (4j^2 - 1)
            end
            b[i] = @dec128 (1 - 2tj) / s
        end
    elseif k == 2
        for i in eachindex(b)
            tj = 0
            th = @dec128 π * i / (s+1)
            for j in 1:div(s+1,2)
                tj += @dec128 sin((2j-1) * th) / (2j - 1)
            end
            b[i] = @dec128 (2tj * sin(th)) / (s+1)
        end
    else
        error("GaussChebyshevQuadrature: k=", k, " not supported.")
    end

    Quadrature(2s-2,c,b)
end


function get_nodes_chebychev_lobatto(s, T=Float64)
    x = zeros(T, s)
    for i in eachindex(x)
        x[i] = @dec128 ( cos( (i-1)*π / (s-1) ) + 1 ) / 2
    end
    return x[end:-1:1]
end

function LobattoChebyshevQuadrature(s, T=Float64)
    c = get_nodes_chebychev_lobatto(s)
    b = zeros(c)

    # TODO compute weights

    error("LobattoChebyshevQuadrature not implemented, yet!")
end


function ClenshawCurtisQuadrature(s, T=Float64)
    c, b = FastTransforms.clenshawcurtis(s, zero(T), zero(T))
    c = flipdim(c,1)
    shift!(b,c)
    Quadrature(s^2, c, b)
end
