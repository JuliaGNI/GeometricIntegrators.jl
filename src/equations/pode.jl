
"Partitioned Ordinary Differential Equation"
immutable PODE{T} <: Equation{T}
    d::Int
    n::Int
    f::Function
    g::Function
    t₀::T
    q₀::Array{T, 2}
    p₀::Array{T, 2}

    function PODE(d, n, f, g, t₀, q₀, p₀)
        @assert d == size(q₀,1) == size(p₀,1)
        @assert n == size(q₀,2) == size(p₀,2)
        @assert T == eltype(q₀) == eltype(p₀)

        if ndims(q₀) == 1
            q₀ = reshape(q₀, d, n)
        end

        if ndims(p₀) == 1
            p₀ = reshape(p₀, d, n)
        end

        new(d, n, f, g, t₀, q₀, p₀)
    end
end


function PODE{T}(f::Function, g::Function, t₀::Real, q₀::DenseArray{T}, p₀::DenseArray{T})
    @assert size(q₀) == size(p₀)
    PODE{T}(size(q₀, 1), size(q₀, 2), f, g, t₀, q₀, p₀)
end

function PODE{T}(f::Function, g::Function, q₀::DenseArray{T}, p₀::DenseArray{T})
    PODE(f, g, 0, q₀, p₀)
end

Base.hash(ode::PODE, h::UInt) = hash(ode.d, hash(ode.n, hash(ode.f, hash(ode.g, hash(ode.t₀, hash(ode.q₀, hash(ode.p₀, h)))))))
Base.:(==){T1, T2}(ode1::PODE{T1}, ode2::PODE{T2}) = (ode1.d == ode2.d
                                                   && ode1.n == ode2.n
                                                   && ode1.f == ode2.f
                                                   && ode1.g == ode2.g
                                                   && ode1.t₀ == ode2.t₀
                                                   && ode1.q₀ == ode2.q₀
                                                   && ode1.p₀ == ode2.p₀)
