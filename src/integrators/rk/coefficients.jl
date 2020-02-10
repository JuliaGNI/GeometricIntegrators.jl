
@define HeaderCoefficientsRK begin
    name::Symbol
    o::Int
    s::Int
end

@define CoefficientsRK begin
    a::Matrix{T}
    b::Vector{T}
    c::Vector{T}
end

"Holds the coefficients of a Runge-Kutta method."
struct CoefficientsRK{T} <: AbstractCoefficients{T}
    @HeaderCoefficientsRK
    @CoefficientsRK
    â::Matrix{T}
    b̂::Vector{T}
    ĉ::Vector{T}

    function CoefficientsRK{T}(name,o,s,a,b,c) where {T}
        new(name,o,s,a,b,c,zero(a),zero(b),zero(c))
    end

    function CoefficientsRK{T}(name::Symbol, o::Int, s::Int, a, b, c, â, b̂, ĉ) where {T <: Real}
        @assert s > 0 "Number of stages must be > 0"
        @assert s==size(a,1)==size(a,2)==length(b)==length(c)
        @assert s==size(â,1)==size(â,2)==length(b̂)==length(ĉ)

        if !get_config(:tab_compensated_summation)
            â = zero(â)
            b̂ = zero(b̂)
            ĉ = zero(ĉ)
        end

        new(name,o,s,a,b,c,â,b̂,ĉ)
    end
end

function CoefficientsRK(T::Type, name::Symbol, order::Int, a::AbstractArray{CT,2}, b::AbstractArray{CT,1}, c::AbstractArray{CT,1}) where {CT}
    a̅ = Matrix{T}(a)
    b̅ = Vector{T}(b)
    c̅ = Vector{T}(c)

    if get_config(:tab_compensated_summation)
        â = Matrix{T}(a-Matrix{eltype(a)}(a̅))
        b̂ = Vector{T}(b-Vector{eltype(b)}(b̅))
        ĉ = Vector{T}(c-Vector{eltype(c)}(c̅))
    else
        â = zero(a̅)
        b̂ = zero(b̅)
        ĉ = zero(c̅)
    end

    CoefficientsRK{T}(name, order, length(c), a̅, b̅, c̅, â, b̂, ĉ)
end

function CoefficientsRK(name::Symbol, order::Int, a::AbstractArray{T,2}, b::AbstractArray{T,1}, c::AbstractArray{T,1}) where {T}
    CoefficientsRK{T}(name, order, length(c), a, b, c)
end

function CoefficientsRK(name::Symbol, order::Int,
        a::AbstractArray{T,2}, b::AbstractArray{T,1}, c::AbstractArray{T,1},
        â::AbstractArray{T,2}, b̂::AbstractArray{T,1}, ĉ::AbstractArray{T,1}) where {T}
    CoefficientsRK{T}(name, order, length(c), a, b, c, â, b̂, ĉ)
end

Base.hash(tab::CoefficientsRK, h::UInt) = hash(tab.o, hash(tab.s, hash(tab.a, hash(tab.b, hash(tab.c, hash(:CoefficientsRK, h))))))

Base.:(==)(tab1::CoefficientsRK, tab2::CoefficientsRK) = (tab1.o == tab2.o
                                                       && tab1.s == tab2.s
                                                       && tab1.a == tab2.a
                                                       && tab1.b == tab2.b
                                                       && tab1.c == tab2.c)

Base.isequal(tab1::CoefficientsRK{T1}, tab2::CoefficientsRK{T2}) where {T1,T2} = (tab1 == tab2 && T1 == T2 && typeof(tab1) == typeof(tab2))

"Print Runge-Kutta coefficients."
function Base.show(io::IO, tab::CoefficientsRK)
    print(io, "Runge-Kutta Coefficients ", tab.name, " with ", tab.s, " stages and order ", tab.o)
    print(io, "  a = ", tab.a)
    print(io, "  b = ", tab.b)
    print(io, "  c = ", tab.c)
end





function get_symplectic_conjugate_coefficients(a::Matrix{T}, b::Vector{T}, a̅::Matrix{T}) where {T}
    @assert size(a) == size(a̅)
    @assert length(b) == size(a,1) == size(a,2)

    for i in 1:size(a̅, 1)
        for j in 1:size(a̅, 2)
            a̅[i,j] = b[j] * ( one(T) - a[j,i] / b[i] )
        end
    end
    nothing
end


function get_symplectic_conjugate_coefficients(coeff::CoefficientsRK)
    a̅ = zero(coeff.a)
    ã = zero(coeff.a)

    get_symplectic_conjugate_coefficients(coeff.a, coeff.b, a̅)

    if all(x -> x≠0, coeff.b̂)
        get_symplectic_conjugate_coefficients(coeff.â, coeff.b̂, ã)
    end

    CoefficientsRK(coeff.name, coeff.s, a̅, coeff.b, coeff.c, ã, coeff.b̂, coeff.ĉ)
end


function symplecticize(coeff::CoefficientsRK; name=nothing, T=Float64)
    name == nothing ? Symbol(string(coeff.name)*"S") : nothing
    a̅ = zero(coeff.a)
    get_symplectic_conjugate_coefficients(coeff.a, coeff.b, a̅)
    CoefficientsRK(T, name, coeff.o, 0.5*(coeff.a + a̅), coeff.b, coeff.c)
end


function compute_symplecticity_error(coeff::CoefficientsRK)
    [coeff.b[i] * coeff.a[i,j] + coeff.b[j] * coeff.a[j,i] - coeff.b[i] * coeff.b[j] for i in 1:size(coeff.a,1), j in 1:size(coeff.a,2)]
end


function check_symplecticity(coeff::CoefficientsRK{T}) where {T}
    symplectic = falses(coeff.s, coeff.s)

    for i in 1:size(coeff.a, 1)
        for j in 1:size(coeff.a, 2)
            symplectic[i,j] = isapprox(coeff.b[i] * coeff.a[i,j] + coeff.b[j] * coeff.a[j,i], coeff.b[i] * coeff.b[j], atol=16*eps(T), rtol=16*eps(T))
        end
    end

    symplectic
end


function check_symmetry(coeff::CoefficientsRK{T}) where {T}
    symmetric = falses(coeff.s, coeff.s)

    for i in 1:size(coeff.a, 1)
        for j in 1:size(coeff.a, 2)
            symmetric[i,j] = isapprox(coeff.a[coeff.s+1-i, coeff.s+1-j] + coeff.a[i,j], coeff.b[j], atol=16*eps(T), rtol=16*eps(T))
        end
    end

    symmetric
end


function check_order_conditions_B(coeff::CoefficientsRK{T}, k) where {T}
    local res::T = 0

    for i in 1:coeff.s
        res += coeff.b[i] * coeff.c[i]^(k-1)
    end

    return isapprox(res, 1/k, atol=16*eps(T), rtol=16*eps(T))
end


function check_order_conditions_C(coeff::CoefficientsRK{T}, k) where {T}
    local order  = falses(coeff.s)
    local res::T

    for i in 1:size(coeff.a, 1)
        res = 0
        for j in 1:size(coeff.a, 2)
            res += coeff.a[i,j] * coeff.c[j]^(k-1)
        end
        order[i] = isapprox(res, coeff.c[i]^k/k, atol=16*eps(T), rtol=16*eps(T))
        # println(res - coeff.c[i]^k/k)
    end

    order
end


function check_order_conditions_D(coeff::CoefficientsRK{T}, k) where {T}
    local order  = falses(coeff.s)
    local res::T

    for i in 1:size(coeff.a, 1)
        res = 0
        for j in 1:size(coeff.a, 2)
            res += coeff.b[j] * coeff.c[j]^(k-1) * coeff.a[j,i]
        end
        order[i] = isapprox(res, coeff.b[i] * (1 - coeff.c[i]^k)/k, atol=16*eps(T), rtol=16*eps(T))
        # println(res - coeff.b[i] * (1 - coeff.c[i]^k)/k)
    end

    order
end
