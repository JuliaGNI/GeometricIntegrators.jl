"Holds the tableau of an Specialised Partitioned Additive Runge-Kutta method for Variational systems."
struct AbstractTableauSPARK{IT,DT} <: AbstractTableau{DT}
    name::Symbol
    o::Int
    s::Int
    r::Int
    ρ::Int

    q::CoefficientsARK{DT}
    p::CoefficientsARK{DT}

    q̃::CoefficientsPRK{DT}
    p̃::CoefficientsPRK{DT}

    λ::CoefficientsMRK{DT}

    ω::Matrix{DT}
    δ::Matrix{DT}
    d::Vector{DT}

    function AbstractTableauSPARK{IT,DT}(name, o, s, r, ρ, q, p, q̃, p̃, λ, ω, δ, d) where {IT,DT}
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
        @assert size(ω, 1)==r-ρ
        @assert size(δ, 1)==ρ

        new(name, o, s, r, ρ, q, p, q̃, p̃, λ, ω, δ, d)
    end

    function AbstractTableauSPARK{IT,DT}(name, o, s, r, ρ, q, p, q̃, p̃, λ, ω, δ) where {IT,DT}
        @assert isa(name, Symbol)
        @assert isa(s, Integer)
        @assert isa(r, Integer)
        @assert isa(ρ, Integer)
        @assert isa(o, Integer)

        @assert s > 0 "Number of stages s must be > 0"
        @assert r > 0 "Number of stages r must be > 0"
        @assert ρ ≥ 0 && ρ ≤ r

        @assert s==q.s==p.s==q̃.s==p̃.s
        @assert r==q.r==p.r==q̃.r==p̃.r==λ.r
        @assert size(ω, 1)==r-ρ
        @assert size(δ, 1)==ρ

        new(name, o, s, r, ρ, q, p, q̃, p̃, λ, ω, δ)
    end
end

function AbstractTableauSPARK{IT}(name::Symbol, order::Int,
                         a_q::Matrix{DT}, a_p::Matrix{DT},
                         α_q::Matrix{DT}, α_p::Matrix{DT},
                         a_q̃::Matrix{DT}, a_p̃::Matrix{DT},
                         α_q̃::Matrix{DT}, α_p̃::Matrix{DT},
                         b_q::Vector{DT}, b_p::Vector{DT},
                         β_q::Vector{DT}, β_p::Vector{DT},
                         c_q::Vector{DT}, c_p::Vector{DT},
                         c_λ::Vector{DT}, d_λ::Union{Vector{DT},Matrix{DT}},
                         ω_λ::Matrix{DT}, δ_λ::Matrix{DT}, d::Vector{DT}) where {IT, DT <: Real}

    s = length(c_q)
    r = length(c_λ)
    ρ = size(δ_λ, 1)

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

    q = CoefficientsARK{DT}(name, order, s, r, a_q, b_q, c_q, α_q, β_q)
    p = CoefficientsARK{DT}(name, order, s, r, a_p, b_p, c_p, α_p, β_p)
    q̃ = CoefficientsPRK{DT}(name, order, s, r, a_q̃, c_λ, α_q̃)
    p̃ = CoefficientsPRK{DT}(name, order, s, r, a_p̃, c_λ, α_p̃)
    λ = CoefficientsMRK{DT}(name, r, d_λ, c_λ)

    AbstractTableauSPARK{IT,DT}(name, order, s, r, ρ, q, p, q̃, p̃, λ, ω_λ, δ_λ, d)
end


function AbstractTableauSPARK{IT}(name::Symbol, order::Int,
                         a_q::Matrix{DT}, a_p::Matrix{DT},
                         α_q::Matrix{DT}, α_p::Matrix{DT},
                         a_q̃::Matrix{DT}, a_p̃::Matrix{DT},
                         α_q̃::Matrix{DT}, α_p̃::Matrix{DT},
                         b_q::Vector{DT}, b_p::Vector{DT},
                         β_q::Vector{DT}, β_p::Vector{DT},
                         c_q::Vector{DT}, c_p::Vector{DT},
                         c_λ::Vector{DT}, d_λ::Union{Vector{DT},Matrix{DT}},
                         ω_λ::Matrix{DT}, δ_λ::Matrix{DT}) where {IT, DT <: Real}

    s = length(c_q)
    r = length(c_λ)
    ρ = size(δ_λ, 1)

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

    q = CoefficientsARK{DT}(name, order, s, r, a_q, b_q, c_q, α_q, β_q)
    p = CoefficientsARK{DT}(name, order, s, r, a_p, b_p, c_p, α_p, β_p)
    q̃ = CoefficientsPRK{DT}(name, order, s, r, a_q̃, c_λ, α_q̃)
    p̃ = CoefficientsPRK{DT}(name, order, s, r, a_p̃, c_λ, α_p̃)
    λ = CoefficientsMRK{DT}(name, r, d_λ, c_λ)

    AbstractTableauSPARK{IT,DT}(name, order, s, r, ρ, q, p, q̃, p̃, λ, ω_λ, δ_λ)
end

# TODO function readAbstractTableauSPARKFromFile(dir::AbstractString, name::AbstractString)
