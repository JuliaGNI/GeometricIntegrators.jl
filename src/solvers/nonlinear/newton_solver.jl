
struct NewtonSolver{T, FT, TJ, TL} <: AbstractNewtonSolver{T}
    @newton_solver_variables
    function NewtonSolver{T,FT,TJ,TL}(x, F!, Jparams, linear_solver) where {T,FT,TJ,TL}
        J  = zeros(linear_solver.A)
        x₀ = zeros(x)
        x₁ = zeros(x)
        y₁ = zeros(x)
        y₀ = zeros(x)
        δx = zeros(x)
        δy = zeros(x)

        nls_params = NonlinearSolverParameters(T)
        nls_status = NonlinearSolverStatus{T}(length(x))

        new(x, J, x₀, x₁, y₀, y₁, δx, δy, F!, Jparams, linear_solver,
            nls_params, nls_status)
    end
end


function NewtonSolver(x::Vector, F!::Function; J=nothing)
    T = eltype(x)
    n = length(x)
    Jparams = getJacobianParameters(J, F!, T, n)
    linear_solver = getLinearSolver(T, n)
    NewtonSolver{T, typeof(F!), typeof(Jparams), typeof(linear_solver)}(x, F!, Jparams, linear_solver)
end


function solve!(s::NewtonSolver{T}; n::Int=0) where {T}
    local nmax::Int = n > 0 ? nmax = n : s.params.nmax

    s.F!(s.x, s.linear.b)
    residual_initial!(s.status, s.linear.b)
    s.status.i  = 0

    if s.status.rₐ ≥ s.params.atol²
        for s.status.i = 1:nmax
            computeJacobian(s.x, s.linear.A, s.Jparams)
            factorize!(s.linear)
            scale!(s.linear.b, -one(T))
            solve!(s.linear)
            simd_xpy!(s.linear.b, s.x)
            s.F!(s.x, s.linear.b)
            residual!(s.status, s.linear.b)

            if solverConverged(s.status, s.params) && s.status.i ≥ s.params.nmin && !(n > 0)
                break
            end
        end
    end
end
