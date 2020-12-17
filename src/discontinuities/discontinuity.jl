
struct Discontinuity{T<:AbstractFloat, PT<:PathIntegral, QN}
    path::PT
    quadrature::QuadratureRule{T,QN}
end

function Discontinuity(path::PT, quadrature::QuadratureRule{T,QN}) where {T,PT,QN}
    Discontinuity{T,PT,QN}(path, quadrature)
end
