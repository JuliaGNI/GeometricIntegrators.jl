
@define HeaderCoefficientsARK begin
    @HeaderCoefficientsRK
    r::Int
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


"Holds the coefficients of a SPARK method."
struct CoefficientsSPARK{T,N} <: AbstractCoefficients{T}
    @HeaderCoefficientsRK

    σ::Int

    a::Tuple{Vararg{Matrix{T},N}}
    b::Tuple{Vararg{Vector{T},N}}
    c::Vector{T}

    function CoefficientsSPARK(name::Symbol, o::Int, s::Int, σ::Int, a::Tuple{Vararg{Matrix{T},N}}, b::Tuple{Vararg{Vector{T},N}}, c::Vector{T}) where {T,N}
        @assert T <: Real
        @assert isa(name, Symbol)
        @assert isa(o, Integer)
        @assert isa(s, Integer)
        @assert isa(σ, Integer)
        @assert s > 0 "Number of stages s must be > 0"
        @assert σ > 0 "Number of stages σ must be > 0"
        @assert s==length(c)

        for α in a
            @assert s==size(α,1)
            @assert s==size(α,2) || r==size(α,2)
        end

        for β in b
            @assert s==length(β) || r==length(β)
        end

        new{T,N}(name,o,s,σ,a,b,c)
    end
end

Base.hash(tab::CoefficientsSPARK, h::UInt) = hash(tab.o, hash(tab.s, hash(tab.r, hash(tab.a, hash(tab.b, hash(tab.c, hash(:CoefficientsSPARK, h)))))))

Base.:(==)(tab1::CoefficientsSPARK, tab2::CoefficientsSPARK) = (tab1.o == tab2.o
                                                             && tab1.s == tab2.s
                                                             && tab1.σ == tab2.σ
                                                             && tab1.a == tab2.a
                                                             && tab1.b == tab2.b
                                                             && tab1.c == tab2.c)

"Print SPARK coefficients."
function Base.show(io::IO, tab::CoefficientsSPARK)
    print(io, "SPARK Coefficients ", tab.name, "with ", tab.s, " internal stages, ", tab.σ, " projective stages and order ", tab.o)
    print(io, "  a = ", tab.a)
    print(io, "  b = ", tab.b)
    print(io, "  c = ", tab.c)
end
