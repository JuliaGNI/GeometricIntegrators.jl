
abstract AbstractNewtonSolver{T} <: NonlinearSolver{T}


@define newton_solver_variables begin
    z::Vector{T}

    F::Function
    J::Function

    linear::LinearSolver{T}

    params::NonlinearSolverParameters{T}
    status::NonlinearSolverStatus{T}
end


function setInitialConditions!{T}(s::AbstractNewtonSolver{T}, z₀::Vector{T})
    s.z[:] = z₀
end
