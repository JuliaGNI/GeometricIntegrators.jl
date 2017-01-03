
"Holds the tableau of a additive Runge-Kutta method."
immutable TableauARK{T} <: AbstractTableau{T}
    name::Symbol
    o::Int
    s::Int
    r::Int

    a_q::Matrix{T}
    α_q::Matrix{T}

    a_q̃::Matrix{T}
    α_q̃::Matrix{T}

    b_q::Vector{T}
    β_q::Vector{T}

    c_q::Vector{T}
    c_λ::Vector{T}

    function TableauARK(name, o, s, r,
                         a_q, α_q, a_q̃, α_q̃,
                         b_q, β_q, c_q, c_λ)
        @assert T <: Real
        @assert isa(name, Symbol)
        @assert isa(s, Integer)
        @assert isa(r, Integer)
        @assert isa(o, Integer)
        @assert s > 0 "Number of stages s must be > 0"
        @assert r > 0 "Number of stages r must be > 0"
        @assert s==size(a_q,1)==size(a_q,2)==length(b_q)==length(c_q)
        @assert s==size(α_q,1)==length(β_q)
        @assert r==size(α_q,2)==size(α_p,2)
        @assert r==length(c_λ)
        @assert r==size(a_q̃,1)==size(α_q̃,1)==size(α_q̃,2)
        @assert s==size(a_q̃,2)
        new(name, o, s, r, a_q, α_q, a_q̃, α_q̃, b_q, β_q, c_q, c_λ)
    end
end

# TODO Add external constructor for TableauARK.

# TODO function readTableauARKFromFile(dir::AbstractString, name::AbstractString)
