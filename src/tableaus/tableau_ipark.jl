
"Holds the tableau of an implicit partitioned additive Runge-Kutta method."
immutable TableauIPARK{T} <: AbstractTableau{T}
    name::Symbol
    o::Int
    s::Int
    r::Int

    a_q::Matrix{T}
    a_p::Matrix{T}
    α_q::Matrix{T}
    α_p::Matrix{T}

    a_q̃::Matrix{T}
    a_p̃::Matrix{T}
    α_q̃::Matrix{T}
    α_p̃::Matrix{T}

    b_q::Vector{T}
    b_p::Vector{T}
    β_q::Vector{T}
    β_p::Vector{T}

    c_q::Vector{T}
    c_p::Vector{T}
    c_λ::Vector{T}

    d_λ::Vector{T}

    function TableauIPARK(name, o, s, r,
                          a_q, a_p, α_q, α_p, a_q̃, a_p̃, α_q̃, α_p̃,
                          b_q, b_p, β_q, β_p,
                          c_q, c_p, c_λ, d_λ)
        @assert T <: Real
        @assert isa(name, Symbol)
        @assert isa(s, Integer)
        @assert isa(r, Integer)
        @assert isa(o, Integer)
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
        new(name, o, s, r,
            a_q, a_p, α_q, α_p, a_q̃, a_p̃, α_q̃, α_p̃,
            b_q, b_p, β_q, β_p,
            c_q, c_p, c_λ, d_λ)
    end
end

function TableauIPARK{T}(name::Symbol, order::Int,
                        a_q::Matrix{T}, a_p::Matrix{T},
                        α_q::Matrix{T}, α_p::Matrix{T},
                        a_q̃::Matrix{T}, a_p̃::Matrix{T},
                        α_q̃::Matrix{T}, α_p̃::Matrix{T},
                        b_q::Vector{T}, b_p::Vector{T},
                        β_q::Vector{T}, β_p::Vector{T},
                        c_q::Vector{T}, c_p::Vector{T},
                        c_λ::Vector{T}, d_λ::Vector{T})
    @assert length(c_q) == length(c_p)
    TableauIPARK{T}(name, order, length(c_q), length(c_λ),
                    a_q, a_p, α_q, α_p, a_q̃, a_p̃, α_q̃, α_p̃,
                    b_q, b_p, β_q, β_p,
                    c_q, c_p, c_λ, d_λ)
end

# TODO function readTableauIPARKFromFile(dir::AbstractString, name::AbstractString)
