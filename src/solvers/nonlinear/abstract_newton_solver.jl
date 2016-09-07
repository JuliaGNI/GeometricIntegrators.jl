
abstract AbstractNewtonSolver{T} <: NonlinearSolver{T}


@define newton_solver_variables begin
    x::Vector{T}

    F::Function
    J::Function

    linear::LinearSolver{T}

    params::NonlinearSolverParameters{T}
    status::NonlinearSolverStatus{T}
end


function setInitialConditions!{T}(s::AbstractNewtonSolver{T}, x₀::Vector{T})
    s.x[:] = x₀
end
