
using FastGaussQuadrature
using FastTransforms


struct Quadrature{T,N}
    order::Int
    nodes::Vector{T}
    weights::Vector{T}
end

function Quadrature{T}(order, nodes::Vector{T}, weights::Vector{T})
    @assert length(nodes) == length(weights)
    Quadrature{T, length(nodes)}(order, nodes, weights)
end

length{T,N}(Q::Quadrature{T,N}) = N
nnodes{T,N}(Q::Quadrature{T,N}) = N
order(Q::Quadrature) = Q.order
nodes(Q::Quadrature) = Q.nodes
weights(Q::Quadrature) = Q.weights


function shift!(b,c)
    # scale from [-1,+1] to [0,1]
    b .= b ./ 2
    c .= (c .+ 1) ./ 2
end


function GaussLegendreQuadrature(s, T=eltype(one()))
    c, b = gausslegendre(s)
    shift!(b,c)
    Quadrature(s^2, c, b)
end


function GaussLobattoQuadrature(s, T=eltype(one()))
    c, b = gausslobatto(s)
    shift!(b,c)
    Quadrature(s^2-2, c, b)
end


function ClenshawCurtisQuadrature(s, T=eltype(one()))
    b, c = FastTransforms.clenshawcurtis(s, zero(T), zero(T))
    b = flipdim(b,1)
    shift!(b,c)
    Quadrature(s^2, c, b)
end
