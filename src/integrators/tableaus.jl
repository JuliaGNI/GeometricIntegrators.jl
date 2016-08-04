
"Tableau: Holds the information for the various methods' tableaus."
abstract Tableau{Name, T}

"TableauRK: Holds the tableau of a Runge-Kutta method."
abstract TableauRK{Name, S, T} <: Tableau{Name, T}

"TableauERK: Holds the tableau of an explicit Runge-Kutta method."
type TableauERK{Name, S, T} <: TableauRK{Name, S, T}
    order::Integer
    a::Matrix{T}
    b::Vector{T}
    c::Vector{T}

    function TableauERK(order,a,b,c)
        @assert T <: Real
        @assert isa(Name, Symbol)
        @assert isa(S, Integer)
        @assert isa(order, Integer)
        @assert S==size(a,1)==size(a,2)==length(b)==length(c)
        @assert c[1] == 0
        @assert istrilstrict(a)
        new(order,a,b,c)
    end
end

function TableauERK{T}(name::Symbol, order::Int,
                       a::Matrix{T}, b::Vector{T}, c::Vector{T})
    TableauERK{name, length(c), T}(order, a, b, c)
end


"TableauIRK: Holds the tableau of a linearly implicit Runge-Kutta method."
type TableauIRK{Name, S, T} <: TableauRK{Name, S, T}
    order::Integer
    a::Matrix{T}
    b::Vector{T}
    c::Vector{T}

    function TableauIRK(order,a,b,c)
        @assert T <: Real
        @assert isa(Name, Symbol)
        @assert isa(S, Integer)
        @assert isa(order, Integer)
        @assert S==size(a,1)==size(a,2)==length(b)==length(c)
        @assert istril(a)
        new(order,a,b,c)
    end
end

function TableauIRK{T}(name::Symbol, order::Int,
                       a::Matrix{T}, b::Vector{T}, c::Vector{T})
    TableauIRK{name, length(c), T}(order, a, b, c)
end


"TableauIRK: Holds the tableau of a nonlinearly implicit Runge-Kutta method."
type TableauNLIRK{Name, S, T} <: TableauRK{Name, S, T}
    order::Integer
    a::Matrix{T}
    b::Vector{T}
    c::Vector{T}

    function TableauNLIRK(order,a,b,c)
        @assert T <: Real
        @assert isa(Name, Symbol)
        @assert isa(S, Integer)
        @assert isa(order, Integer)
        @assert S==size(a,1)==size(a,2)==length(b)==length(c)
        new(order,a,b,c)
    end
end

function TableauNLIRK{T}(name::Symbol, order::Int,
                         a::Matrix{T}, b::Vector{T}, c::Vector{T})
    TableauNLIRK{name, length(c), T}(order, a, b, c)
end


"TableauPRK: Holds the tableau of a partitioned Runge-Kutta method."
# TODO Need explicit and implicit version?
type TableauPRK{Name, S, T} <: TableauRK{Name, S, T}
    order::Integer
    a_q::Matrix{T}
    a_p::Matrix{T}
    b_q::Vector{T}
    b_p::Vector{T}
    c_q::Vector{T}
    c_p::Vector{T}

    function TableauPRK(order, a_q, a_p, b_q, b_p, c_q, c_p)
        @assert T <: Real
        @assert isa(Name, Symbol)
        @assert isa(S, Integer)
        @assert isa(order, Integer)
        @assert S==size(a_q,1)==size(a_q,2)==length(b_q)==length(c_q)
        @assert S==size(a_p,1)==size(a_p,2)==length(b_p)==length(c_p)
        new(order, a_q, a_p, b_q, b_p, c_q, c_p)
    end
end

function TableauPRK{T}(name::Symbol, order::Int,
                       a_q::Matrix{T}, a_p::Matrix{T},
                       b_q::Vector{T}, b_p::Vector{T},
                       c_q::Vector{T}, c_p::Vector{T})
    @assert length(c_q)==length(c_p)
    TableauPRK{name, length(c_q), T}(order, a_q, a_p, b_q, b_p, c_q, c_p)
end


"TableauSPARK: Holds the tableau of a spezialized partitioned additive
 Runge-Kutta method."
type TableauSPARK{Name, S, R, T} <: TableauRK{Name, S, T}
    order::Integer
    a_q::Matrix{T}
    a_p::Matrix{T}
    α_q::Matrix{T}
    α_p::Matrix{T}

    a_qᵠ::Matrix{T}
    a_pᵠ::Matrix{T}
    α_qᵠ::Matrix{T}
    α_pᵠ::Matrix{T}

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

    function TableauSPARK(order, a_q, a_p, α_q, α_p, a_qᵠ, a_pᵠ, α_qᵠ, α_pᵠ,
                                 b_q, b_p, β_q, β_p,
                                 c_q, c_p, c_λ,
                                 ω_q, ω_p, ω_λ)
        # TODO Make ω_q, ω_p, ω_λ optional arguments.
        @assert T <: Real
        @assert isa(Name, Symbol)
        @assert isa(S, Integer)
        @assert isa(R, Integer)
        @assert isa(order, Integer)
        @assert S==size(a_q,1)==size(a_q,2)==length(b_q)==length(c_q)
        @assert S==size(a_p,1)==size(a_p,2)==length(b_p)==length(c_p)
        @assert S==size(α_q,1)==size(α_p,1)==length(β_q)==length(β_p)
        @assert R==size(α_q,2)==size(α_p,2)
        @assert R==length(c_λ)
        @assert R==size(a_qᵠ,1)==size(a_pᵠ,1)
        @assert R==size(α_qᵠ,1)==size(α_qᵠ,2)
        @assert R==size(α_pᵠ,1)==size(α_pᵠ,2)
        @assert S==size(a_qᵠ,2)==size(a_pᵠ,2)
        # TODO Add assertions on ω_q, ω_p, ω_λ to be (S-1)x(S) or (R-1)x(R) if set.
        new(order, a_q, a_p, α_q, α_p, a_qᵠ, a_pᵠ, α_qᵠ, α_pᵠ,
                   b_q, b_p, β_q, β_p,
                   c_q, c_p, c_λ,
                   ω_q, ω_p, ω_λ)
    end
end


"TableauGLM: Holds the tableau of a general linear method."
type TableauGLM{Name, S, R, T} <: Tableau{Name, T}
    order::Integer
    a::Matrix{T}
    b::Matrix{T}
    u::Matrix{T}
    v::Matrix{T}
    c::Vector{T}

    function TableauGLM(order, a, b, u, v, c)
        @assert T <: Real
        @assert isa(Name, Symbol)
        @assert isa(S, Integer)
        @assert isa(R, Integer)
        @assert isa(order, Integer)
        @assert S==length(c)
        @assert S==size(a,1)==size(a,2)==length(u,1)==length(b,2)
        @assert R==size(v,1)==size(v,2)==length(u,2)==length(b,1)
        new(order, a, b, u, v, c)
    end
end


# TODO Add TableauAGLM.
# TODO Add TableauTSRK.
# TODO Add TableauSplitting.
