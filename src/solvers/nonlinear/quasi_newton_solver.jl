
immutable QuasiNewtonSolver{T} <: AbstractNewtonSolver{T}
    @newton_solver_variables
    function QuasiNewtonSolver(x, F, J, linear_solver, nmax, atol, rtol, stol)
        new(x, F, J, linear_solver, NonlinearSolverParameters{T}(nmax, atol, rtol, stol), NonlinearSolverStatus{T}())
    end
end


function QuasiNewtonSolver(x::AbstractVector, F::Function; J=nothing, linear_solver=nothing, nmax=DEFAULT_nmax, atol=DEFAULT_atol, rtol=DEFAULT_rtol, stol=DEFAULT_stol, ϵ=DEFAULT_ϵ, autodiff=false)
    J = getComputeJacobianFunction(J, F, ϵ, autodiff)
    linear_solver = getLinearSolver(eltype(x), length(x), linear_solver)
    QuasiNewtonSolver{eltype(x)}(x, F, J, linear_solver, nmax, atol, rtol, stol)
end


function solve!{T}(s::QuasiNewtonSolver{T})
    s.F(s.x, s.linear.b)
    s.linear.b[:] *= -1
    s.status.rₐ = residual_absolute(s.linear.b)
    s.status.r₀ = s.status.rₐ

    if s.status.r₀ ≥ s.params.atol²
        s.J(s.x, s.linear.A)
        factorize!(s.linear)
        for s.status.i = 1:s.params.nmax
            solve!(s.linear)
            s.x[:] += s.linear.b[:]
            s.status.rᵣ = residual_relative(s.linear.b, s.x)
            s.F(s.x, s.linear.b)
            s.linear.b[:] *= -1
            s.status.rₛ = s.status.rₐ
            s.status.rₐ = residual_absolute(s.linear.b)
            s.status.rₛ = abs(s.status.rₛ - s.status.rₐ)/s.status.r₀

            if s.status.rₐ < s.params.atol² || s.status.rᵣ < s.params.rtol² || s.status.rₛ < s.params.stol²
                break
            end
        end
    end
end
