
using GeomDAE.Fields: Field, Q, P, V, F

"Special Ordinary Differential Equation"
immutable SODE{T} <: Equation{T}
    d::Int

    q₀::Array{T, 1}
    p₀::Array{T, 1}
    t₀::T

    q::Function
    p::Function
    v::Function
    f::Function

    fields::Tuple{Vararg{Field}}
    functs::Tuple{Vararg{Field}}

    function SODE(d, q₀, p₀, t₀, q, p, v, f)
        @assert d == length(q₀) == length(p₀)
        @assert T == eltype(q₀) == eltype(p₀)

        fields = []
        functs = []

        for (field, symb) in ((q, Q()),
                              (p, P()),
                              (v, V()),
                              (f, F()))
            if field == do_nothing
                push!(fields, symb)
            else
                push!(functs, symb)
            end
        end

        @assert length(fields) == length(functs) == 2

        new(d, q₀, p₀, t₀, q, p, v, f, tuple(fields...), tuple(functs...))
    end
end

function SODE{T}(q₀::Vector{T}, p₀::Vector{T}, t₀::Real=0;
                 q::Function=do_nothing, p::Function=do_nothing,
                 v::Function=do_nothing, f::Function=do_nothing)
    @assert length(q₀) == length(p₀)
    SODE{T}(length(q₀), q₀, p₀, t₀, q, p, v, f)
end
