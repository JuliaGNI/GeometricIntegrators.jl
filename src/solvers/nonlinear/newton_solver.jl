
type NewtonSolver{T} <: NonlinearSolver{T}
    x::Vector{T}

    F::Function
    J::Function

    linear::LinearSolver{T}

    atol::T
    rtol::T
    stol::T

    atol²::T
    rtol²::T
    stol²::T

    n_max::Int

    i::Int
    r₀::T
    rₐ::T
    rᵣ::T
    rₛ::T

    function NewtonSolver(x, F, J, linear_solver, atol, rtol, stol, n_max)
        # @assert F <: Function
        # @assert J <: Function
        #
        # @assert linear_solver <: LinearSolver{T}

        @assert atol > 0
        @assert rtol > 0
        @assert stol > 0
        @assert n_max > 0

        new(x, F, J, linear_solver, atol, rtol, stol, atol^2, rtol^2, stol^2, n_max, 0, 0, 0, 0, 0)
    end
end


function NewtonSolver{T}(x::Vector{T}, F::Function; J=nothing, linear_solver=nothing, atol=1E-15, rtol=1E-20, stol=1E-20, n_max=100, ϵ=1E-6, autodiff=false)
    if J == nothing
        function computeJacobian(tx::Vector, A::Matrix)
            computeJacobianFD(tx, A, F, ϵ)
        end

        # TODO
        # if autodiff
        #     function computeJacobian(x::Vector, J::Matrix)
        #         computeJacobianAD(x, J, F, ϵ)
        #     end
        # else
        #     function computeJacobian(x::Vector, J::Matrix)
        #         computeJacobianFD(x, J, F, ϵ)
        #     end
        # end

        J = computeJacobian
    end

    if linear_solver == nothing
        n = length(x)
        linear_solver = LUSolver(zeros(T, n, n), zeros(T, n))
    end

    NewtonSolver{T}(x, F, J, linear_solver, atol, rtol, stol, n_max)
end


function solve!{T}(s::NewtonSolver{T})
    s.F(s.x, s.linear.b)
    s.linear.b[:] *= -1
    s.rₐ = residual_absolute(s.linear.b)
    s.r₀ = s.rₐ

    # δx::Vector{T} = zeros(T, length(s.x))
    #
    # s.i = 0
    # println(s.x, ", ", δx, ", ", s.i, ", ", s.rₐ,", ",  s.rᵣ,", ",  s.rₛ)

    if s.r₀ ≥ s.atol²
        for s.i = 1:s.n_max
            s.J(s.x, s.linear.A)
            factorize!(s.linear)
            solve!(s.linear)
            s.x[:] += s.linear.b[:]
            # δx[:] = s.linear.b
            s.rᵣ = residual_relative(s.linear.b, s.x)
            s.F(s.x, s.linear.b)
            s.linear.b[:] *= -1
            s.rₛ = s.rₐ
            s.rₐ = residual_absolute(s.linear.b)
            s.rₛ = abs(s.rₛ - s.rₐ)/s.r₀

            # println(s.x, ", ", s.i, ", ", s.rₐ,", ",  s.rᵣ,", ",  s.rₛ)
            # println(s.x, ", ", δx, ", ", s.i, ", ", s.rₐ,", ",  s.rᵣ,", ",  s.rₛ)

            if s.rₐ < s.atol² || s.rᵣ < s.rtol² || s.rₛ < s.stol²
            # if s.rₐ < s.atol²
                break
            end
        end
    end
end


function solve!{T}(s::NewtonSolver{T}, x₀::Vector{T})
    s.x[:] = x₀
    solve!(s)
end
