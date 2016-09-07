
immutable NewtonSolver{T} <: NonlinearSolver{T}
    z::Vector{T}

    F::Function
    J::Function

    linear::LinearSolver{T}

    atol::T
    rtol::T
    stol::T

    atol²::T
    rtol²::T
    stol²::T

    nmax::Int

    status::NonlinearSolverStatus{T}

    function NewtonSolver(z, F, J, linear_solver, atol, rtol, stol, nmax)
        # @assert F <: Function
        # @assert J <: Function
        #
        # @assert linear_solver <: LinearSolver{T}

        @assert atol > 0
        @assert rtol > 0
        @assert stol > 0
        @assert nmax > 0

        new(z, F, J, linear_solver, atol, rtol, stol, atol^2, rtol^2, stol^2, nmax, NonlinearSolverStatus{T}())
    end
end


function NewtonSolver(z::AbstractVector, F::Function; J=nothing, linear_solver=nothing, atol=DEFAULT_atol, rtol=DEFAULT_rtol, stol=DEFAULT_stol, nmax=DEFAULT_nmax, ϵ=DEFAULT_ϵ, autodiff=false)
    T = eltype(z)
    n = length(z)

    if J == nothing
        if autodiff
            function computeJacobianADF{T}(x::Vector{T}, A::Matrix{T})
                computeJacobianAD(x, A, F, ϵ)
            end
            J = computeJacobianADF
        else
            function computeJacobianFDF{T}(x::Vector{T}, A::Matrix{T})
                computeJacobianFD(x, A, F, ϵ)
            end
            J = computeJacobianFDF
        end
    end

    if linear_solver == nothing
        linear_solver = LUSolver(zeros(T, n, n), zeros(T, n))
    end

    NewtonSolver{T}(z, F, J, linear_solver, atol, rtol, stol, nmax)
end


function solve!{T}(s::NewtonSolver{T})
    s.F(s.z, s.linear.b)
    s.linear.b[:] *= -1
    s.status.rₐ = residual_absolute(s.linear.b)
    s.status.r₀ = s.status.rₐ

    if s.status.r₀ ≥ s.atol²
        for s.status.i = 1:s.nmax
            s.J(s.z, s.linear.A)
            factorize!(s.linear)
            solve!(s.linear)
            s.z[:] += s.linear.b[:]
            s.status.rᵣ = residual_relative(s.linear.b, s.z)
            s.F(s.z, s.linear.b)
            s.linear.b[:] *= -1
            s.status.rₛ = s.status.rₐ
            s.status.rₐ = residual_absolute(s.linear.b)
            s.status.rₛ = abs(s.status.rₛ - s.status.rₐ)/s.status.r₀

            if s.status.rₐ < s.atol² || s.status.rᵣ < s.rtol² || s.status.rₛ < s.stol²
                break
            end
        end
    end
end


function solve!{T}(s::NewtonSolver{T}, z₀::Vector{T})
    s.z[:] = z₀
    solve!(s)
end
