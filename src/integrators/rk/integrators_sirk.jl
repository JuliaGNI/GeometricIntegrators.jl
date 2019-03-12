
"Holds the tableau of a singly implicit Runge-Kutta method."
struct TableauSIRK{T} <: AbstractTableauIRK{T}
    # TODO
end

function TableauSIRK(name::Symbol, order::Int, a::Matrix{T}, b::Vector{T}, c::Vector{T}) where {T}
    TableauSIRK{T}(name, order, length(c), a, b, c)
end


"Singly implicit Runge-Kutta integrator."
struct IntegratorSIRK{DT,TT,FT} <: Integrator{DT,TT}
    equation::ODE{DT,TT,FT}
    tableau::TableauSIRK{TT}
    Δt::TT

    x::Vector{Vector{TwicePrecision{DT}}}
    X::Vector{Vector{DT}}
    Y::Vector{Vector{DT}}
    F::Vector{Vector{DT}}

    function IntegratorSIRK{DT,TT,FT}(equation, tableau, Δt) where {DT,TT,FT}
        D = equation.d
        M = equation.n
        S = tableau.s
        X = create_internal_stage_vector(DT, D, S)
        Y = create_internal_stage_vector(DT, D, S)
        F = create_internal_stage_vector(DT, D, S)
        new(equation, tableau, Δt, create_solution_vector(DT,D,M), X, Y, F)
    end
end

function IntegratorSIRK(equation::Equation, tableau::TableauSIRK, Δt)
    T = eltype(equation.q₀)
    IntegratorSIRK{T}(equation, tableau, Δt)
end

"Integrate ODE with singly implicit Runge-Kutta integrator."
function integrate!(int::IntegratorSIRK, s::SolutionODE)
    # TODO
end

"Integrate partitioned ODE with singly implicit Runge-Kutta integrator."
function integrate!(int::IntegratorSIRK, s::SolutionPODE)
    # TODO
end
