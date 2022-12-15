"""
Compute p(x) where p is the unique polynomial of degree length(xi),
such that p(x[i]) = y[i]) for all i.
Call with
```julia
aitken_neville!(x::AbstractVector, t::Real, ti::AbstractVector, xi::AbstractMatrix)
```
where
* `x`:  evaluation value
* `t`:  evaluation point
* `ti`: interpolation nodes
* `xi`: interpolation values
"""
function aitken_neville!(x::AbstractVector, t::TT, ti::AbstractVector{TT}, xi::AbstractMatrix) where {TT}
    @assert length(ti) == size(xi,2)
    @assert length(x)  == size(xi,1)

    for j in eachindex(ti)
        for i in 1:(length(ti)-j)
            for k in axes(x,1)
                xi[k,i] = xi[k,i+1] + (xi[k,i] - xi[k,i+1]) * (ti[i+j] - t) / (ti[i+j] - ti[i])
            end
        end
    end
    for k in eachindex(x)
        x[k] = xi[k,1]
    end
end
