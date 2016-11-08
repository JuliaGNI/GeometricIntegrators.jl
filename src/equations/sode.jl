
using GeomDAE.Fields: Field, Q, P, V, F

"Special Ordinary Differential Equation"
immutable SODE{T} <: Equation{T}
    d::Int
    n::Int

    q₀::Array{T, 2}
    p₀::Array{T, 2}
    t₀::T

    q::Function
    p::Function
    v::Function
    f::Function

    fields::Tuple{Vararg{Field}}
    functs::Tuple{Vararg{Field}}

    function SODE(d, n, q₀, p₀, t₀, q, p, v, f)
        @assert d == size(q₀,1) == size(p₀,1)
        @assert n == size(q₀,2) == size(p₀,2)
        @assert T == eltype(q₀) == eltype(p₀)

        if ndims(q₀) == 1
            q₀ = reshape(q₀, d, n)
        end

        if ndims(p₀) == 1
            p₀ = reshape(p₀, d, n)
        end

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

        new(d, n, q₀, p₀, t₀, q, p, v, f, tuple(fields...), tuple(functs...))
    end
end

function SODE{T}(q₀::DenseArray{T}, p₀::DenseArray{T}, t₀::Real=0;
                 q::Function=do_nothing, p::Function=do_nothing,
                 v::Function=do_nothing, f::Function=do_nothing)
    @assert size(q₀) == size(p₀)
    SODE{T}(size(q₀, 1), size(q₀, 2), q₀, p₀, t₀, q, p, v, f)
end

Base.hash(ode::SODE, h::UInt) = hash(ode.d, hash(ode.n,
                                hash(ode.t₀, hash(ode.q₀, hash(ode.p₀,
                                hash(ode.q, hash(ode.p, hash(ode.v, hash(ode.f,
                                hash(ode.fields, hash(ode.functs, h)))))))))))

Base.:(==){T1, T2}(ode1::SODE{T1}, ode2::SODE{T2}) = (ode1.d == ode2.d
                                                   && ode1.n == ode2.n
                                                   && ode1.t₀ == ode2.t₀
                                                   && ode1.q₀ == ode2.q₀
                                                   && ode1.p₀ == ode2.p₀
                                                   && ode1.q == ode2.q
                                                   && ode1.p == ode2.p
                                                   && ode1.v == ode2.v
                                                   && ode1.f == ode2.f
                                                   && ode1.fields == ode2.fields
                                                   && ode1.functs == ode2.functs)
