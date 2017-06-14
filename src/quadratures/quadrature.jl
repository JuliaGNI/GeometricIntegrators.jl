
import FastGaussQuadrature
import FastTransforms

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


function shift!(b,c)
    # scale from [-1,+1] to [0,1]
    b .= b ./ 2
    c .= (c .+ 1) ./ 2
end


# function LegendreGaussQuadrature(s, T=Float64)
function GaussLegendreQuadrature(s, T=Float64)
    c, b = FastGaussQuadrature.gausslegendre(s)
    shift!(b,c)
    Quadrature(s^2, c, b)
end


# function LegendreLobattoQuadrature(s, T=Float64)
function GaussLobattoQuadrature(s, T=Float64)
    c, b = FastGaussQuadrature.gausslobatto(s)
    shift!(b,c)
    Quadrature(s^2-2, c, b)
end


# function ChebyshevGaussQuadrature(s, T=Float64)
# function ChebyshevLobattoQuadrature(s, T=Float64)


function ClenshawCurtisQuadrature(s, T=Float64)
    c, b = FastTransforms.clenshawcurtis(s, zero(T), zero(T))
    c = flipdim(c,1)
    shift!(b,c)
    Quadrature(s^2, c, b)
end
