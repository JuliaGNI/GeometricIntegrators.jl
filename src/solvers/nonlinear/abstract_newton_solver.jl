
abstract AbstractNewtonSolver{T} <: NonlinearSolver{T}

@define newton_solver_variables begin
    x::Vector{T}
    J::Matrix{T}

    x₀::Vector{T}
    x₁::Vector{T}
    y₀::Vector{T}
    y₁::Vector{T}
    δx::Vector{T}
    δy::Vector{T}

    Fparams::TF
    Jparams::TJ

    linear::TL

    params::NonlinearSolverParameters{T}
    status::NonlinearSolverStatus{T}
end


function setInitialConditions!{T}(s::AbstractNewtonSolver{T}, x₀::Vector{T})
    s.x[:] = x₀
end
