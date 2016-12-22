
@define HeaderCoefficients begin
    name::Symbol
    o::Int
    s::Int
end

@define CoefficientsRK begin
    a::Matrix{T}
    b::Vector{T}
    c::Vector{T}
end

@define CoefficientsARK begin
    α::Matrix{T}
    β::Vector{T}
end

"Holds the coefficients of a Runge-Kutta method."
immutable CoefficientsRK{T}
    @HeaderCoefficients
    @CoefficientsRK

    function CoefficientsRK(name,o,s,a,b,c)
        @assert T <: Real
        @assert isa(name, Symbol)
        @assert isa(s, Integer)
        @assert isa(o, Integer)
        @assert s > 0 "Number of stages must be > 0"
        @assert s==size(a,1)==size(a,2)==length(b)==length(c)
        new(name,o,s,a,b,c)
    end
end

function CoefficientsRK{T}(name::Symbol, order::Int, a::Matrix{T}, b::Vector{T}, c::Vector{T})
    CoefficientsRK{T}(name, order, length(c), a, b, c)
end

Base.hash(tab::CoefficientsRK, h::UInt) = hash(tab.o, hash(tab.a, hash(tab.b, hash(tab.c, hash(:TableauRK, h)))))

Base.:(==){T1, T2}(tab1::CoefficientsRK{T1}, tab2::CoefficientsRK{T2}) = (tab1.o == tab2.o
                                                             && tab1.s == tab2.s
                                                             && tab1.a == tab2.a
                                                             && tab1.b == tab2.b
                                                             && tab1.c == tab2.c)

Base.isequal{T1, T2}(tab1::CoefficientsRK{T1}, tab2::CoefficientsRK{T2}) = (tab1 == tab2 && T1 == T2 && typeof(tab1) == typeof(tab2))

"Print Runge-Kutta tableau to standard output."
function show_coefficients{T}(tab::CoefficientsRK{T})
    println("Runge-Kutta Method ", tab.name, "with ", tab.s, " stages and order ", tab.o)
    println("  a = ", tab.a)
    println("  b = ", tab.b)
    println("  c = ", tab.c)
end
