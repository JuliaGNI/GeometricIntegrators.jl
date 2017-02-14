
immutable QuasiNewtonSolver{T, TF, TJ, TL} <: AbstractNewtonSolver{T}
    @newton_solver_variables
    function QuasiNewtonSolver(x, Fparams, Jparams, linear_solver, nmax, atol, rtol, stol)
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

const DEFAULT_NonlinearSolver=QuasiNewtonSolver
const DEFAULT_QUASINEWTON_REFACTORIZE=5
const DEFAULT_LINESEARCH_nmax=10
const DEFAULT_ARMIJO_λ₀ = 1.0
const DEFAULT_ARMIJO_σ₀ = 0.1
const DEFAULT_ARMIJO_σ₁ = 0.5
const DEFAULT_ARMIJO_ϵ  = 0.5


function QuasiNewtonSolver(x::Vector, Fparams::NonlinearFunctionParameters; J=nothing, linear_solver=nothing, nmax=DEFAULT_nmax, atol=DEFAULT_atol, rtol=DEFAULT_rtol, stol=DEFAULT_stol, ϵ=DEFAULT_ϵ, autodiff=false)
    T = eltype(x)
    n = length(x)
    Jparams = getJacobianParameters(J, Fparams, ϵ, T, n, autodiff)
    linear_solver = getLinearSolver(T, n, linear_solver)
    QuasiNewtonSolver{T, typeof(Fparams), typeof(Jparams), typeof(linear_solver)}(x, Fparams, Jparams, linear_solver, nmax, atol, rtol, stol)
end


function solve!{T}(s::QuasiNewtonSolver{T})
    local λ::T
    local λₜ::T
    local y₀norm::T
    local y₁norm::T
    local lsiter::Int
    local p₀::T
    local p₁::T
    local p₂::T

    function_stages!(s.x, s.y₀, s.Fparams)
    residual_initial!(s.status, s.y₀)
    s.status.i  = 0

    if s.status.rₐ ≥ s.params.atol²
        computeJacobian(s.x, s.J, s.Jparams)
        simd_copy!(s.J, s.linear.A)
        factorize!(s.linear)

        for s.status.i = 1:s.params.nmax
            simd_copy!(s.x, s.x₀)
            y₀norm = vecnorm(s.y₀, 2)

            # b = - y₀
            simd_copy_scale!(-one(T), s.y₀, s.linear.b)

            # solve J δx = -f(x)
            solve!(s.linear)

            # TODO Separate line search algorithms into independent modules.

            # δx = b
            simd_copy!(s.linear.b, s.δx)

            # δy = Jδx
            simd_mult!(s.δy, s.J, s.δx)

            # set λ to default initial value
            λ = DEFAULT_ARMIJO_λ₀

            # x₁ = x₀ + λ δx
            simd_waxpy!(s.x₁, λ, s.δx, s.x₀)

            # y₁ = f(x₁)
            function_stages!(s.x₁, s.y₁, s.Fparams)

            # compute norms of solutions
            y₀norm = vecnorm(s.y₀, 2)
            y₁norm = vecnorm(s.y₁, 2)

            if y₁norm < (one(T)-DEFAULT_ARMIJO_σ₀*λ)*y₀norm
                nothing
            else
                # determine coefficients of polynomial p(λ) = p₀ + p₁λ + p₂λ²
                p₀ = y₀norm^2
                p₁ = 2(⋅(s.y₀, s.δy))
                p₂ = (y₁norm^2 - y₀norm^2 - p₁*λ)/(λ^2)

                # compute minimum λₜ of p(λ)
                λₜ = - p₁/(2p₂)

                if λₜ < DEFAULT_ARMIJO_σ₀ * λ
                    λ = DEFAULT_ARMIJO_σ₀ * λ
                elseif λₜ > DEFAULT_ARMIJO_σ₁ * λ
                    λ = DEFAULT_ARMIJO_σ₁ * λ
                else
                    λ = λₜ
                end
            end

            simd_scale!(s.δx, λ)

            # simple Armijo line search
            #
            # λ = DEFAULT_ARMIJO_λ₀
            #
            # # δx = λ b
            # simd_copy_scale!(λ, s.linear.b, s.δx)
            #
            # for lsiter in 1:DEFAULT_LINESEARCH_nmax
            #     # x₁ = x₀ + λ δx
            #     simd_wxpy!(s.x₁, s.δx, s.x₀)
            #
            #     try
            #         # y₁ = f(x₁)
            #         function_stages!(s.x₁, s.y₁, s.Fparams)
            #
            #         # y₁norm = ||y₁||₂
            #         y₁norm = vecnorm(s.y₁, 2)
            #
            #         # exit loop if ||f(x₁)||₂ < (1-λϵ) ||f(x₀)||₂
            #         if y₁norm < (one(T)-λ*DEFAULT_ARMIJO_ϵ)*y₀norm
            #             break
            #         end
            #     catch DomainError
            #         # in case the new function value results in some DomainError
            #         # (e.g., for functions f(x) containing sqrt's or log's),
            #         # decrease λ and retry
            #         println("WARNING: Quasi-Newton Solver encountered Domain Error.")
            #         # if lsiter == DEFAULT_LINESEARCH_nmax
            #         #     λ *= 0
            #         #     simd_scale!(s.δx, 0)
            #         # end
            #     end
            #
            #     λ *= DEFAULT_ARMIJO_σ₁
            #     simd_scale!(s.δx, DEFAULT_ARMIJO_σ₁)
            # end

            simd_xpy!(s.δx, s.x)
            function_stages!(s.x, s.y₀, s.Fparams)
            residual_absolute!(s.status, s.y₀)
            residual_relative!(s.status, s.y₀)

            if solverConverged(s.status, s.params)
                if s.status.i > DEFAULT_nwarn
                    println("WARNING: Quasi-Newton Solver took ", s.status.i, " iterations.")
                end
                break
            end

            if s.status.i % DEFAULT_QUASINEWTON_REFACTORIZE == 0
                computeJacobian(s.x, s.J, s.Jparams)
                simd_copy!(s.J, s.linear.A)
                factorize!(s.linear)
            end
        end
    end
    nothing
end
