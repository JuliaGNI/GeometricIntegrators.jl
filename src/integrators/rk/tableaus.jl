
@define HeaderCoefficientsRK begin
    name::Symbol
    o::Int
    s::Int
end

@define CoefficientsRK begin
    a::Matrix{T}
    b::Vector{T}
    c::Vector{T}

    â::Matrix{T}
    b̂::Vector{T}
    ĉ::Vector{T}
end
