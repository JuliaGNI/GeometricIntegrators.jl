
"Holds the tableau of a General Partitioned Additive Runge-Kutta method."
struct TableauGPARK{T} <: AbstractTableau{T}
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

    ω_q::Matrix{T}
    ω_p::Matrix{T}
    ω_λ::Matrix{T}

    function TableauGPARK{T}(name, o, s, r,
                             a_q, a_p, α_q, α_p, a_q̃, a_p̃, α_q̃, α_p̃,
                             b_q, b_p, β_q, β_p,
                             c_q, c_p, c_λ,
                             ω_q, ω_p, ω_λ) where {T}
        # TODO Make ω_q, ω_p, ω_λ optional arguments.
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
        @assert r==length(c_λ)
        @assert r==size(a_q̃,1)==size(a_p̃,1)
        @assert r==size(α_q̃,1)==size(α_q̃,2)
        @assert r==size(α_p̃,1)==size(α_p̃,2)
        @assert s==size(a_q̃,2)==size(a_p̃,2)
        # TODO Add assertions on ω_q, ω_p, ω_λ to be (S-1)x(S) or (R-1)x(R) if set.
        new(name, o, s, r,
            a_q, a_p, α_q, α_p, a_q̃, a_p̃, α_q̃, α_p̃,
            b_q, b_p, β_q, β_p,
            c_q, c_p, c_λ,
            ω_q, ω_p, ω_λ)
    end
end

# TODO Add external constructor for TableauGPARK.

# TODO function readTableauGPARKFromFile(dir::AbstractString, name::AbstractString)

# TODO function writeTableauToFile{Name, T}(dir::AbstractString, tab::TableauGPARK{Name, T})


"General Partitioned Additive Runge Kutta integrator."
struct IntegratorGPARK{DT,TT,VT,FT,UT,GT,ΦT} <: AbstractIntegratorSPARK{DT,TT}
    equation::PDAE{DT,TT,VT,FT,UT,GT,ΦT}
    tableau::TableauGPARK{TT}
    Δt::TT

    solver::NonlinearSolver{DT}

    x::Array{DT,1}
    y::Array{DT,1}
    X::Array{DT,2}
    Y::Array{DT,2}
    F::Array{DT,2}
end

"Integrate partitioned DAE with General Additive Runge Kutta integrator."
function integrate!(int::IntegratorGPARK, s::SolutionPDAE)
    # TODO
end
