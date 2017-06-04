
"Holds the tableau of an variational special partitioned additive Runge-Kutta method."
struct TableauVSPARK{T} <: AbstractTableau{T}
    name::Symbol
    o::Int
    s::Int
    r::Int
    ρ::Int

    q::CoefficientsARK{T}
    p::CoefficientsARK{T}

    q̃::CoefficientsPRK{T}
    p̃::CoefficientsPRK{T}

    λ::CoefficientsMRK{T}

    ω::Matrix{T}
    d::Vector{T}

    function TableauVSPARK{T}(name, o, s, r, ρ, q, p, q̃, p̃, λ, ω, d) where {T}
        @assert isa(name, Symbol)
        @assert isa(s, Integer)
        @assert isa(r, Integer)
        @assert isa(ρ, Integer)
        @assert isa(o, Integer)

        @assert s > 0 "Number of stages s must be > 0"
        @assert r > 0 "Number of stages r must be > 0"
        @assert ρ > 0 && ρ ≤ r

        @assert s==q.s==p.s==q̃.s==p̃.s==length(d)
        @assert r==q.r==p.r==q̃.r==p̃.r==λ.r
        @assert ρ==size(ω, 1)
        @assert r==size(ω, 2)

        new(name, o, s, r, ρ, q, p, q̃, p̃, λ, ω, d)
    end

    function TableauVSPARK{T}(name, o, s, r, ρ, q, p, q̃, p̃, λ, ω) where {T}
        @assert isa(name, Symbol)
        @assert isa(s, Integer)
        @assert isa(r, Integer)
        @assert isa(ρ, Integer)
        @assert isa(o, Integer)

        @assert s > 0 "Number of stages s must be > 0"
        @assert r > 0 "Number of stages r must be > 0"
        @assert ρ > 0 && ρ ≤ r

        @assert s==q.s==p.s==q̃.s==p̃.s
        @assert r==q.r==p.r==q̃.r==p̃.r==λ.r
        @assert ρ==size(ω, 1)
        @assert r==size(ω, 2)

        new(name, o, s, r, ρ, q, p, q̃, p̃, λ, ω)
    end
end

function TableauVSPARK(name::Symbol, order::Int,
                         a_q::Matrix{T}, a_p::Matrix{T},
                         α_q::Matrix{T}, α_p::Matrix{T},
                         a_q̃::Matrix{T}, a_p̃::Matrix{T},
                         α_q̃::Matrix{T}, α_p̃::Matrix{T},
                         b_q::Vector{T}, b_p::Vector{T},
                         β_q::Vector{T}, β_p::Vector{T},
                         c_q::Vector{T}, c_p::Vector{T},
                         c_λ::Vector{T}, d_λ::Vector{T},
                         ω_λ::Matrix{T}, d::Vector{T}) where {T <: Real}

    s = length(c_q)
    r = length(c_λ)
    ρ = size(ω_λ, 1)

    @assert s > 0 "Number of stages s must be > 0"
    @assert r > 0 "Number of stages r must be > 0"

    @assert s==size(a_q,1)==size(a_q,2)==length(b_q)==length(c_q)
    @assert s==size(a_p,1)==size(a_p,2)==length(b_p)==length(c_p)
    @assert s==size(α_q,1)==size(α_p,1)
    @assert r==size(α_q,2)==size(α_p,2)
    @assert s==length(d)
    @assert r==length(c_λ)==length(d_λ)
    @assert r==size(a_q̃,1)==size(a_p̃,1)
    @assert s==size(a_q̃,2)==size(a_p̃,2)
    @assert r==size(α_q̃,1)==size(α_q̃,2)==length(β_q)
    @assert r==size(α_p̃,1)==size(α_p̃,2)==length(β_p)

    q = CoefficientsARK{T}(name, order, s, r, a_q, b_q, c_q, α_q, β_q)
    p = CoefficientsARK{T}(name, order, s, r, a_p, b_p, c_p, α_p, β_p)
    q̃ = CoefficientsPRK{T}(name, order, s, r, a_q̃, c_λ, α_q̃)
    p̃ = CoefficientsPRK{T}(name, order, s, r, a_p̃, c_λ, α_p̃)
    λ = CoefficientsMRK{T}(name, r, d_λ, c_λ)

    TableauVSPARK{T}(name, order, s, r, ρ, q, p, q̃, p̃, λ, ω_λ, d)
end


function TableauVSPARK(name::Symbol, order::Int,
                         a_q::Matrix{T}, a_p::Matrix{T},
                         α_q::Matrix{T}, α_p::Matrix{T},
                         a_q̃::Matrix{T}, a_p̃::Matrix{T},
                         α_q̃::Matrix{T}, α_p̃::Matrix{T},
                         b_q::Vector{T}, b_p::Vector{T},
                         β_q::Vector{T}, β_p::Vector{T},
                         c_q::Vector{T}, c_p::Vector{T},
                         c_λ::Vector{T}, d_λ::Vector{T},
                         ω_λ::Matrix{T}) where {T <: Real}

    s = length(c_q)
    r = length(c_λ)
    ρ = size(ω_λ, 1)

    @assert s > 0 "Number of stages s must be > 0"
    @assert r > 0 "Number of stages r must be > 0"

    @assert s==size(a_q,1)==size(a_q,2)==length(b_q)==length(c_q)
    @assert s==size(a_p,1)==size(a_p,2)==length(b_p)==length(c_p)
    @assert s==size(α_q,1)==size(α_p,1)
    @assert r==size(α_q,2)==size(α_p,2)
    @assert r==length(c_λ)==length(d_λ)
    @assert r==size(a_q̃,1)==size(a_p̃,1)
    @assert s==size(a_q̃,2)==size(a_p̃,2)
    @assert r==size(α_q̃,1)==size(α_q̃,2)==length(β_q)
    @assert r==size(α_p̃,1)==size(α_p̃,2)==length(β_p)

    q = CoefficientsARK{T}(name, order, s, r, a_q, b_q, c_q, α_q, β_q)
    p = CoefficientsARK{T}(name, order, s, r, a_p, b_p, c_p, α_p, β_p)
    q̃ = CoefficientsPRK{T}(name, order, s, r, a_q̃, c_λ, α_q̃)
    p̃ = CoefficientsPRK{T}(name, order, s, r, a_p̃, c_λ, α_p̃)
    λ = CoefficientsMRK{T}(name, r, d_λ, c_λ)

    TableauVSPARK{T}(name, order, s, r, ρ, q, p, q̃, p̃, λ, ω_λ)
end

# TODO function readTableauVSPARKFromFile(dir::AbstractString, name::AbstractString)
