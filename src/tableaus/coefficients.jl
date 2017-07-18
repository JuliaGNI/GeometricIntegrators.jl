
abstract type AbstractCoefficients{T} end

Base.isequal(coeff1::AbstractCoefficients{T1}, coeff2::AbstractCoefficients{T2}) where {T1,T2} = (coeff1 == coeff2 && T1 == T2 && typeof(coeff1) == typeof(coeff2))


@define HeaderCoefficientsRK begin
    name::Symbol
    o::Int
    s::Int
end

@define HeaderCoefficientsARK begin
    @HeaderCoefficientsRK
    r::Int
end

@define CoefficientsRK begin
    a::Matrix{T}
    b::Vector{T}
    c::Vector{T}
end

@define CoefficientsARK begin
    @CoefficientsRK
    α::Matrix{T}
    β::Vector{T}
end

@define CoefficientsPRK begin
    a::Matrix{T}
    c::Vector{T}
    α::Matrix{T}
end


"Holds the coefficients of a Runge-Kutta method."
struct CoefficientsRK{T} <: AbstractCoefficients{T}
    @HeaderCoefficientsRK
    @CoefficientsRK
    â::Matrix{T}
    b̂::Vector{T}
    ĉ::Vector{T}

    function CoefficientsRK{T}(name,o,s,a,b,c) where {T}
        new(name,o,s,a,b,c,zeros(a),zeros(b),zeros(c))
    end

    function CoefficientsRK{T}(name,o,s,a,b,c,â,b̂,ĉ) where {T}
        @assert T <: Real
        @assert isa(name, Symbol)
        @assert isa(o, Integer)
        @assert isa(s, Integer)
        @assert s > 0 "Number of stages must be > 0"
        @assert s==size(a,1)==size(a,2)==length(b)==length(c)
        @assert s==size(â,1)==size(â,2)==length(b̂)==length(ĉ)

        if !get_config(:tab_compensated_summation)
            â = zeros(â)
            b̂ = zeros(b̂)
            ĉ = zeros(ĉ)
        end
        
        new(name,o,s,a,b,c,â,b̂,ĉ)
    end
end

function CoefficientsRK(T::Type, name::Symbol, order::Int, a::Matrix, b::Vector, c::Vector)
    a̅ = Matrix{T}(a)
    b̅ = Vector{T}(b)
    c̅ = Vector{T}(c)

    â = Matrix{T}(a-Matrix{eltype(a)}(a̅))
    b̂ = Vector{T}(b-Vector{eltype(b)}(b̅))
    ĉ = Vector{T}(c-Vector{eltype(c)}(c̅))

    CoefficientsRK{T}(name, order, length(c), a̅, b̅, c̅, â, b̂, ĉ)
end

function CoefficientsRK(name::Symbol, order::Int, a::Matrix{T}, b::Vector{T}, c::Vector{T}) where {T}
    CoefficientsRK{T}(name, order, length(c), a, b, c)
end

function CoefficientsRK(name::Symbol, order::Int, a::Matrix{T}, b::Vector{T}, c::Vector{T}, â::Matrix{T}, b̂::Vector{T}, ĉ::Vector{T}) where {T}
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


"Holds the coefficients of an additive Runge-Kutta method."
struct CoefficientsARK{T} <: AbstractCoefficients{T}
    @HeaderCoefficientsARK
    @CoefficientsARK

    function CoefficientsARK{T}(name,o,s,r,a,b,c,α,β) where {T}
        @assert T <: Real
        @assert isa(name, Symbol)
        @assert isa(o, Integer)
        @assert isa(s, Integer)
        @assert isa(r, Integer)
        @assert s > 0 "Number of stages s must be > 0"
        @assert r > 0 "Number of stages r must be > 0"
        @assert s==size(a,1)==size(a,2)==size(α,1)==length(b)==length(c)
        @assert r==size(α,2)==length(β)
        new(name,o,s,r,a,b,c,α,β)
    end
end

function CoefficientsARK(name::Symbol, order::Int, a::Matrix{T}, b::Vector{T}, c::Vector{T}, α::Matrix{T}, β::Vector{T}) where {T}
    CoefficientsARK{T}(name, order, length(b), length(β), a, b, c, α, β)
end

Base.hash(tab::CoefficientsARK, h::UInt) = hash(tab.o, hash(tab.s, hash(tab.r, hash(tab.a, hash(tab.b, hash(tab.c, hash(tab.α, hash(tab.β, hash(:CoefficientsARK, h)))))))))

Base.:(==)(tab1::CoefficientsARK, tab2::CoefficientsARK) = (tab1.o == tab2.o
                                                         && tab1.s == tab2.s
                                                         && tab1.r == tab2.r
                                                         && tab1.a == tab2.a
                                                         && tab1.b == tab2.b
                                                         && tab1.c == tab2.c
                                                         && tab1.α == tab2.α
                                                         && tab1.β == tab2.β)

"Print additive Runge-Kutta coefficients."
function Base.show(io::IO, tab::CoefficientsARK)
    print(io, "Additive Runge-Kutta Coefficients ", tab.name, "with ", tab.s, " internal stages, ", tab.r, " projective stages and order ", tab.o)
    print(io, "  a = ", tab.a)
    print(io, "  b = ", tab.b)
    print(io, "  c = ", tab.c)
    print(io, "  α = ", tab.α)
    print(io, "  β = ", tab.β)
end


"Holds the coefficients of a projective Runge-Kutta method."
struct CoefficientsPRK{T} <: AbstractCoefficients{T}
    @HeaderCoefficientsARK
    @CoefficientsPRK

    function CoefficientsPRK{T}(name,o,s,r,a,c,α) where {T}
        @assert T <: Real
        @assert isa(name, Symbol)
        @assert isa(o, Integer)
        @assert isa(s, Integer)
        @assert isa(r, Integer)
        @assert s > 0 "Number of stages s must be > 0"
        @assert r > 0 "Number of stages r must be > 0"
        @assert r==size(a,1)==size(α,1)==size(α,2)==length(c)
        @assert s==size(a,2)
        new(name,o,s,r,a,c,α)
    end
end

function CoefficientsPRK(name::Symbol, order::Int, a::Matrix{T}, c::Vector{T}, α::Matrix{T}) where {T}
    CoefficientsPRK{T}(name, order, size(a,2), length(c), a, c, α)
end

Base.hash(tab::CoefficientsPRK, h::UInt) = hash(tab.o, hash(tab.s, hash(tab.r, hash(tab.a, hash(tab.c, hash(tab.α, hash(:CoefficientsPRK, h)))))))

Base.:(==)(tab1::CoefficientsPRK, tab2::CoefficientsPRK) = (tab1.o == tab2.o
                                                         && tab1.s == tab2.s
                                                         && tab1.r == tab2.r
                                                         && tab1.a == tab2.a
                                                         && tab1.c == tab2.c
                                                         && tab1.α == tab2.α)

"Print projective Runge-Kutta coefficients."
function Base.show(io::IO, tab::CoefficientsPRK)
    print(io, "Projective Runge-Kutta Coefficients ", tab.name, "with ", tab.s, " internal stages, ", tab.r, " projective stages and order ", tab.o)
    print(io, "  a = ", tab.a)
    print(io, "  c = ", tab.c)
    print(io, "  α = ", tab.α)
end


"Holds the multiplier Runge-Kutta coefficients."
struct CoefficientsMRK{T}
    name::Symbol
    r::Int
    b::Vector{T}
    c::Vector{T}

    function CoefficientsMRK{T}(name,r,b,c) where {T}
        @assert T <: Real
        @assert isa(name, Symbol)
        @assert isa(r, Integer)
        @assert r > 0 "Number of stages r must be > 0"
        @assert r==length(b)==length(c)
        new(name,r,b,c)
    end
end

function CoefficientsMRK(name::Symbol, b::Vector{T}, c::Vector{T}) where {T}
    CoefficientsMRK{T}(name, length(c), b, c)
end

Base.hash(tab::CoefficientsMRK, h::UInt) = hash(tab.r, hash(tab.b, hash(tab.c, hash(:CoefficientsMRK, h))))

Base.:(==)(tab1::CoefficientsMRK, tab2::CoefficientsMRK) = (tab1.r == tab2.r
                                                         && tab1.b == tab2.b
                                                         && tab1.c == tab2.c)

"Print multiplier Runge-Kutta coefficients."
function Base.show(io::IO, tab::CoefficientsMRK)
    print(io, "Multiplier Runge-Kutta coefficients ", tab.name, "with ", tab.r, " projective stages")
    print(io, "  b = ", tab.b)
    print(io, "  c = ", tab.c)
end


"Holds the coefficients of a projected Gauss-Legendre Runge-Kutta method."
struct CoefficientsPGLRK{T} <: AbstractCoefficients{T}
    @HeaderCoefficientsRK
    @CoefficientsRK

    P::Matrix{T}
    Q::Matrix{T}
    X::Matrix{T}
    W::Matrix{T}
    A::Matrix{T}

    function CoefficientsPGLRK{T}(name,o,s,a,b,c,P,X,W) where {T}
        @assert T <: Real
        @assert isa(name, Symbol)
        @assert isa(o, Integer)
        @assert isa(s, Integer)
        @assert s ≥ 2 "Number of stages must be ≥ 2"
        @assert s==size(a,1)==size(a,2)==length(b)==length(c)
        @assert s==size(P,1)==size(P,2)
        @assert s==size(X,1)==size(X,2)
        @assert s==size(W,1)==size(W,2)

        Q = inv(P)
        A = zeros(a)
        B = zeros(a)

        simd_mult!(B, W, Q)
        simd_mult!(A, P, B)

        new(name,o,s,a,b,c,P,Q,X,W,A)
    end
end

function CoefficientsPGLRK(name::Symbol, order::Int, a::Matrix{T}, b::Vector{T}, c::Vector{T}, P::Matrix{T}, X::Matrix{T}, W::Matrix{T}) where {T}
    CoefficientsPGLRK{T}(name, order, length(c), a, b, c, P, X, W)
end

Base.hash(tab::CoefficientsPGLRK, h::UInt) = hash(tab.o, hash(tab.s, hash(tab.a, hash(tab.b, hash(tab.c, hash(tab.P, hash(tab.Q, hash(tab.X, hash(tab.W, hash(tab.A, h))))))))))

Base.:(==)(tab1::CoefficientsPGLRK, tab2::CoefficientsPGLRK) = (tab1.o == tab2.o
                                                             && tab1.s == tab2.s
                                                             && tab1.a == tab2.a
                                                             && tab1.b == tab2.b
                                                             && tab1.c == tab2.c
                                                             && tab1.P == tab2.P
                                                             && tab1.Q == tab2.Q
                                                             && tab1.X == tab2.X
                                                             && tab1.W == tab2.W
                                                             && tab1.A == tab2.A)

Base.isequal(tab1::CoefficientsPGLRK{T1}, tab2::CoefficientsPGLRK{T2}) where {T1,T2} = (tab1 == tab2 && T1 == T2 && typeof(tab1) == typeof(tab2))

"Print Runge-Kutta coefficients."
function Base.show(io::IO, tab::CoefficientsPGLRK)
    print(io, "Projected Gauss-Legendre Runge-Kutta Coefficients ", tab.name, " with ", tab.s, " stages and order ", tab.o)
    print(io, "  a = ", tab.a)
    print(io, "  b = ", tab.b)
    print(io, "  c = ", tab.c)
    print(io, "  P = ", tab.P)
    print(io, "  Q = ", tab.Q)
    print(io, "  X = ", tab.X)
    print(io, "  W = ", tab.W)
    print(io, "  A = ", tab.A)
end
