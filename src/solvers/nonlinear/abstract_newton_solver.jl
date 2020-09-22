
abstract type AbstractNewtonSolver{T} <: NonlinearSolver{T} end

@define newton_solver_variables begin
    x::Vector{T}
    J::Matrix{T}

    x₀::Vector{T}
    x₁::Vector{T}
    y₀::Vector{T}
    y₁::Vector{T}
    δx::Vector{T}
    δy::Vector{T}

    F!::FT
    Jparams::TJ

    linear::TL

    params::NonlinearSolverParameters{T}
    status::NonlinearSolverStatus{T}
end

function setInitialConditions!(s::AbstractNewtonSolver{T}, x₀::Vector{T}) where {T}
    s.x[:] = x₀
end

status(solver::AbstractNewtonSolver) = solver.status
params(solver::AbstractNewtonSolver) = solver.params

computeJacobian(s::AbstractNewtonSolver) = computeJacobian(s.x, s.J, s.Jparams)

function check_jacobian(s::AbstractNewtonSolver)
    println("Condition Number of Jacobian: ", cond(s.J))
    println("Determinant of Jacobian:      ", det(s.J))
    println("minimum(|Jacobian|):          ", minimum(abs.(s.J)))
    println("maximum(|Jacobian|):          ", maximum(abs.(s.J)))
    println()
end

function print_jacobian(s::AbstractNewtonSolver)
    display(s.J)
    println()
end
