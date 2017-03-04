
immutable NewtonSolver{T, TF, TJ, TL} <: AbstractNewtonSolver{T}
    @newton_solver_variables
    function NewtonSolver(x, Fparams, Jparams, linear_solver, nmax, atol, rtol, stol)
        J  = zeros(linear_solver.A)
        x₀ = zeros(x)
        x₁ = zeros(x)
        y₁ = zeros(x)
        y₀ = zeros(x)
        δx = zeros(x)
        δy = zeros(x)
        new(x, J, x₀, x₁, y₀, y₁, δx, δy, Fparams, Jparams, linear_solver, NonlinearSolverParameters{T}(nmax, atol, rtol, stol), NonlinearSolverStatus{T}(length(x)))
    end
end


function NewtonSolver(x::Vector, Fparams::NonlinearFunctionParameters; J=nothing, linear_solver=nothing, nmax=DEFAULT_nmax, atol=DEFAULT_atol, rtol=DEFAULT_rtol, stol=DEFAULT_stol, ϵ=DEFAULT_ϵ, autodiff=false)
    T = eltype(x)
    n = length(x)
    Jparams = getJacobianParameters(J, Fparams, ϵ, T, n, autodiff)
    linear_solver = getLinearSolver(T, n, linear_solver)
    NewtonSolver{T, typeof(Fparams), typeof(Jparams), typeof(linear_solver)}(x, Fparams, Jparams, linear_solver, nmax, atol, rtol, stol)
end


function solve!{T}(s::NewtonSolver{T})
    function_stages!(s.x, s.linear.b, s.Fparams)
    residual_initial!(s.status, s.linear.b)
    s.status.i  = 0

    if s.status.rₐ ≥ s.params.atol²
        for s.status.i = 1:s.params.nmax
            computeJacobian(s.x, s.linear.A, s.Jparams)
            factorize!(s.linear)
            scale!(s.linear.b, -one(T))
            solve!(s.linear)
            simd_xpy!(s.linear.b, s.x)
            function_stages!(s.x, s.linear.b, s.Fparams)
            residual!(s.status, s.linear.b)

            if solverConverged(s.status, s.params)
                break
            end
        end
    end
    nothing
end
