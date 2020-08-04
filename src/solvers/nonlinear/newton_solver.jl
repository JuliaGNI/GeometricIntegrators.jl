
struct NewtonSolver{T, FT, TJ, TL} <: AbstractNewtonSolver{T}
    @newton_solver_variables
    function NewtonSolver{T,FT,TJ,TL}(x, F!, Jparams, linear_solver) where {T,FT,TJ,TL}
        J  = zero(linear_solver.A)
        x₀ = zero(x)
        x₁ = zero(x)
        y₁ = zero(x)
        y₀ = zero(x)
        δx = zero(x)
        δy = zero(x)

        nls_params = NonlinearSolverParameters(T)
        nls_status = NonlinearSolverStatus{T}(length(x))

        new(x, J, x₀, x₁, y₀, y₁, δx, δy, F!, Jparams, linear_solver,
            nls_params, nls_status)
    end
end


function NewtonSolver(x::AbstractVector{T}, F!::Function; J!::Union{Function,Nothing}=nothing) where {T}
    n = length(x)
    Jparams = getJacobianParameters(J!, F!, T, n)
    linear_solver = getLinearSolver(x)
    NewtonSolver{T, typeof(F!), typeof(Jparams), typeof(linear_solver)}(x, F!, Jparams, linear_solver)
end


function solve!(s::NewtonSolver{T}; n::Int=0) where {T}
    local nmax::Int = n > 0 ? nmax = n : s.params.nmax

    s.F!(s.x, s.y₀)
    residual_initial!(s.status, s.x, s.y₀)
    s.status.i  = 0

    if s.status.rₐ ≥ s.params.atol²
        for s.status.i = 1:nmax
            computeJacobian(s)
            s.linear.A .= s.J
            s.linear.b .= -s.y₀
            factorize!(s.linear)
            solve!(s.linear)
            s.δx .= s.linear.b
            s.x .+= s.δx
            s.F!(s.x, s.y₀)
            residual!(s.status, s.δx, s.x, s.y₀)

            if check_solver_converged(s.status, s.params) && s.status.i ≥ s.params.nmin && !(n > 0)
                break
            end
        end
    end
end
