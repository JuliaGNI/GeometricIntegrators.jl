
"""
Compute p(x) where p is the unique polynomial of degree length(xi),
such that p(x[i]) = y[i]) for all i.

    ti: interpolation nodes
    xi: interpolation values
    t:  evaluation point
    x:  evaluation value
"""
function aitken_neville(ti::Vector{TT}, xi::Matrix, t::TT, x::Vector) where {TT}
    @assert length(ti) == size(xi,2)
    @assert length(x)  == size(xi,1)

    for j in axes(t)
        for i in 1:(length(t)-j)
            for k in axes(x,1)
                xi[k,i] = xi[k,i+1] + (xi[k,i] - xi[k,i+1]) * (ti[i+j] - t) / (ti[i+j] - ti[i])
            end
        end
    end
    for k in eachindex(x)
        x[k] = xi[k,1]
    end
end
