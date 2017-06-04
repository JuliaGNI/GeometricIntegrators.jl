
"Holds the tableau of an partitioned additive Runge-Kutta method."
struct TableauPARK{T} <: AbstractTableau{T}
    name::Symbol
    o::Int
    s::Int
    r::Int

    q::CoefficientsARK{T}
    p::CoefficientsARK{T}

    q̃::CoefficientsPRK{T}
    p̃::CoefficientsPRK{T}

    λ::CoefficientsMRK{T}

    function TableauPARK{T}(name, o, s, r, q, p, q̃, p̃, λ) where {T}
        @assert isa(name, Symbol)
        @assert isa(s, Integer)
        @assert isa(r, Integer)
        @assert isa(o, Integer)

        @assert s > 0 "Number of stages s must be > 0"
        @assert r > 0 "Number of stages r must be > 0"

        @assert s==q.s==p.s==q̃.s==p̃.s
        @assert r==q.r==p.r==q̃.r==p̃.r==λ.r

        new(name, o, s, r, q, p, q̃, p̃, λ)
    end
end

function TableauPARK(name::Symbol, order::Int,
                     a_q::Matrix{T}, a_p::Matrix{T},
                     α_q::Matrix{T}, α_p::Matrix{T},
                     a_q̃::Matrix{T}, a_p̃::Matrix{T},
                     α_q̃::Matrix{T}, α_p̃::Matrix{T},
                     b_q::Vector{T}, b_p::Vector{T},
                     β_q::Vector{T}, β_p::Vector{T},
                     c_q::Vector{T}, c_p::Vector{T},
                     c_λ::Vector{T}, d_λ::Vector{T}) where {T <: Real}

    s = length(c_q)
    r = length(c_λ)

    @assert s > 0 "Number of stages s must be > 0"
    @assert r > 0 "Number of stages r must be > 0"

    @assert s==size(a_q,1)==size(a_q,2)==length(b_q)==length(c_q)
    @assert s==size(a_p,1)==size(a_p,2)==length(b_p)==length(c_p)
    @assert s==size(α_q,1)==size(α_p,1)==length(β_q)==length(β_p)
    @assert r==size(α_q,2)==size(α_p,2)
    @assert r==length(c_λ)==length(d_λ)
    @assert r==size(a_q̃,1)==size(a_p̃,1)
    @assert r==size(α_q̃,1)==size(α_q̃,2)
    @assert r==size(α_p̃,1)==size(α_p̃,2)
    @assert s==size(a_q̃,2)==size(a_p̃,2)

    q = CoefficientsARK{T}(name, order, s, r, a_q, b_q, c_q, α_q, β_q)
    p = CoefficientsARK{T}(name, order, s, r, a_p, b_p, c_p, α_p, β_p)
    q̃ = CoefficientsPRK{T}(name, order, s, r, a_q̃, c_λ, α_q̃)
    p̃ = CoefficientsPRK{T}(name, order, s, r, a_p̃, c_λ, α_p̃)
    λ = CoefficientsMRK{T}(name, r, d_λ, c_λ)

    TableauPARK{T}(name, order, s, r, q, p, q̃, p̃, λ)
end

# TODO function readTableauPARKFromFile(dir::AbstractString, name::AbstractString)
