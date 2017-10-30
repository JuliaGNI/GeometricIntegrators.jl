
struct Flux{T<:AbstractFloat, PT<:FluxPath, QN}
    path::PT
    quadrature::Quadrature{T,QN}
end

function Flux(path::PT, quadrature::Quadrature{T,QN}) where {T,PT,QN}
    Flux{T,PT,QN}(path, quadrature)
end
