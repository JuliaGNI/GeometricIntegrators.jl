

immutable ODE{T} <: Equation{T}
    d::Int
    n::Int
    f::Function
    t₀::T
    q₀::Array{T, 2}

    function ODE(d, n, f, t₀, q₀)
        @assert d == size(q₀,1)
        @assert n == size(q₀,2)
        @assert T == eltype(q₀)
        @assert ndims(q₀) ∈ (1,2)

        if ndims(q₀) == 1
            q₀ = reshape(q₀, d, n)
        end

        new(d, n, f, t₀, q₀)
    end
end

function ODE{T}(f::Function, t₀::Real, q₀::DenseArray{T})
    ODE{T}(size(q₀, 1), size(q₀, 2), f, t₀, q₀)
end

function ODE{T}(f::Function, q₀::DenseArray{T})
    ODE(f, 0, q₀)
end

Base.hash(ode::ODE, h::UInt) = hash(ode.d, hash(ode.n, hash(ode.f, hash(ode.t₀, hash(ode.q₀, h)))))
Base.:(==){T1, T2}(ode1::ODE{T1}, ode2::ODE{T2}) = (ode1.d == ode2.d
                                                 && ode1.n == ode2.n
                                                 && ode1.f == ode2.f
                                                 && ode1.t₀ == ode2.t₀
                                                 && ode1.q₀ == ode2.q₀)
