
abstract AbstractNewtonSolver{T} <: NonlinearSolver{T}

@define newton_solver_variables begin
    x::Vector{T}

    Fparams::TF
    Jparams::TJ

    linear::TL

    params::NonlinearSolverParameters{T}
    status::NonlinearSolverStatus{T}
end


function setInitialConditions!{T}(s::AbstractNewtonSolver{T}, x₀::Vector{T})
    s.x[:] = x₀
end
