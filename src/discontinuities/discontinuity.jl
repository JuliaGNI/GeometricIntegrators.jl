
struct Discontinuity{T<:AbstractFloat, PT<:PathIntegral, QN}
    path::PT
    quadrature::Quadrature{T,QN}
end

function Discontinuity(path::PT, quadrature::Quadrature{T,QN}) where {T,PT,QN}
    Discontinuity{T,PT,QN}(path, quadrature)
end
