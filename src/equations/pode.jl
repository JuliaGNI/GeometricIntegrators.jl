
"Partitioned Ordinary Differential Equation"
immutable PODE{dType, tType, fqType, fpType} <: Equation{dType, tType}
    d::Int
    n::Int
    f::fqType
    g::fpType
    t₀::tType
    q₀::Array{dType, 2}
    p₀::Array{dType, 2}

    function PODE(d, n, f, g, t₀, q₀, p₀)
        @assert d == size(q₀,1) == size(p₀,1)
        @assert n == size(q₀,2) == size(p₀,2)
        @assert dType == eltype(q₀) == eltype(p₀)

        if ndims(q₀) == 1
            q₀ = reshape(q₀, d, n)
        end

        if ndims(p₀) == 1
            p₀ = reshape(p₀, d, n)
        end

        new(d, n, f, g, t₀, q₀, p₀)
    end
end


function PODE{DT,TT,FQ,FP}(f::FQ, g::FP, t₀::TT, q₀::DenseArray{DT}, p₀::DenseArray{DT})
    @assert size(q₀) == size(p₀)
    PODE{DT,TT,FQ,FP}(size(q₀, 1), size(q₀, 2), f, g, t₀, q₀, p₀)
end

function PODE(f, g, q₀, p₀)
    PODE(f, g, zero(eltype(q₀)), q₀, p₀)
end

Base.hash(ode::PODE, h::UInt) = hash(ode.d, hash(ode.n, hash(ode.f, hash(ode.g, hash(ode.t₀, hash(ode.q₀, hash(ode.p₀, h)))))))
Base.:(==){DT1, DT2, TT1, TT2, FQ1, FQ2, FP1, FP2}(ode1::PODE{DT1,TT1,FQ1,FP1}, ode2::PODE{DT2,TT2,FQ2,FP2}) = (
                                ode1.d == ode2.d
                             && ode1.n == ode2.n
                             && ode1.f == ode2.f
                             && ode1.g == ode2.g
                             && ode1.t₀ == ode2.t₀
                             && ode1.q₀ == ode2.q₀
                             && ode1.p₀ == ode2.p₀)
