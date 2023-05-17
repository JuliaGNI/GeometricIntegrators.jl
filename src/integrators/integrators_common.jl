
# general helper functions for integrators.

_default_solver() = SimpleSolvers.NewtonSolver

"""
Create nonlinear solver object for a system of `N` equations with data type `DT`.
The function ``f(x)=0`` to be solved for is determined by a julia function
`function_stages!(x, b, params)`, where `x` is the current solution and `b` is
the output vector, s.th. ``b = f(x)``. `params` are a set of parameters depending
on the equation and integrator that is used.
The solver type is obtained from the config dictionary (`:nls_solver`).
"""
function create_nonlinear_solver(DT, N, params, caches, solver::Type{<:NonlinearSolver}=_default_solver(), F::Function=function_stages!)
    # create solution vector for nonlinear solver
    x = zeros(DT, N)

    # create wrapper function f!(x,b) that calls `F(x, b, params)`
    # with the appropriate `params`
    f! = (b,x) -> F(x, b, params, caches)

    # create nonlinear solver with solver type obtained from config dictionary
    solver(x, zero(x), f!; linesearch = Backtracking(), config = Options(min_iterations = 1, x_abstol = 8eps(), f_abstol = 8eps()))
end

function create_nonlinear_solver(DT, N, params, caches, i::Int, solver::Type{<:NonlinearSolver}=_default_solver(), F::Function=function_stages!)
    # create solution vector for nonlinear solver
    x = zeros(DT, N)

    # create wrapper function f(x,b) that calls `F(x, b, params)`
    # with the appropriate `params`
    f! = (b,x) -> F(x, b, params, caches, i)

    # create nonlinear solver with solver type obtained from config dictionary
    solver(x, zero(x), f!; linesearch = Backtracking(), config = Options(min_iterations = 1, x_abstol = 8eps(), f_abstol = 8eps()))
end

function create_nonlinear_solver_with_jacobian(DT, N, params, caches, solver::Type{<:NonlinearSolver}=_default_solver(), F=function_stages!, J=jacobian!)
    # create solution vector for nonlinear solver
    x = zeros(DT, N)

    # create wrapper function f!(x,b) that calls `F(x, b, params)`
    # with the appropriate `params`
    f! = (b,x) -> F(x, b, params, caches)

    # create wrapper function j!(x,df) that calls `J(x, df, params)`
    # with the appropriate `params`
    j! = (df,x) -> J(x, df, caches[DT], params)

    # create nonlinear solver with solver type obtained from config dictionary
    solver(x, zero(x), f!; J! = j!, linesearch = Backtracking(), config = Options(min_iterations = 1, x_abstol = 8eps(), f_abstol = 8eps()))
end


function update_solution!(x::Vector{T}, xₑᵣᵣ::Vector{T}, ẋ::Matrix{T}, b::AbstractVector, Δt) where {T}
    @assert length(x) == length(xₑᵣᵣ)
    @assert length(x) == size(ẋ, 1)
    @assert length(b) == size(ẋ, 2)

    for k in axes(ẋ, 1)
        for i in axes(ẋ, 2)
            x[k], xₑᵣᵣ[k] = compensated_summation(Δt * b[i] * ẋ[k,i], x[k], xₑᵣᵣ[k])
        end
    end
end

function update_solution!(x::Vector{T}, xₑᵣᵣ::Vector{T}, ẋ::Vector{Vector{T}}, b::AbstractVector, Δt) where {T}
    @assert length(b) == length(ẋ)
    @assert length(x) == length(ẋ[1])
    @assert length(x) == length(xₑᵣᵣ)

    for i in eachindex(ẋ)
        for k in eachindex(ẋ[i])
            x[k], xₑᵣᵣ[k] = compensated_summation(Δt * b[i] * ẋ[i][k], x[k], xₑᵣᵣ[k])
        end
    end
end

function update_solution!(x::SolutionVector{T}, ẋ::Matrix{T}, b::AbstractVector, Δt) where {T}
    @assert length(x) == size(ẋ, 1)
    @assert length(b) == size(ẋ, 2)

    local Δx::eltype(x)

    for k in axes(ẋ, 1)
        Δx = 0
        for i in axes(ẋ, 2)
            Δx += b[i] * ẋ[k,i]
        end
        x[k] += Δt * Δx
    end
end


function update_solution!(x::SolutionVector{T}, ẋ::Vector{Vector{T}}, b::AbstractVector, Δt) where {T}
    @assert length(b) == length(ẋ)
    @assert length(x) == length(ẋ[1])

    local Δx::T

    for k in eachindex(x)
        Δx = 0
        for i in eachindex(ẋ)
            Δx += b[i] * ẋ[i][k]
        end
        x[k] += Δt * Δx
    end
end

function update_solution!(x::Vector{T}, xₑᵣᵣ::Vector{T}, ẋ::Union{Matrix{T},Vector{Vector{T}}}, b::AbstractVector, b̂::AbstractVector, Δt) where {T}
    update_solution!(x, xₑᵣᵣ, ẋ, b, Δt)
    update_solution!(x, xₑᵣᵣ, ẋ, b̂, Δt)
end

function update_solution!(x::SolutionVector{T}, ẋ::Union{Matrix{T},Vector{Vector{T}}}, b::AbstractVector, b̂::AbstractVector, Δt) where {T}
    update_solution!(x, ẋ, b, Δt)
    update_solution!(x, ẋ, b̂, Δt)
end

function update_multiplier!(λ::SolutionVector{T}, Λ::Vector{Vector{T}}, b::AbstractVector) where {T}
    for Λᵢ in Λ
        @assert length(λ) == length(Λᵢ)
    end
    local t::T
    @inbounds for i in eachindex(λ)
        t = zero(T)
        for j in eachindex(b,Λ)
            t += b[j] * Λ[j][i]
        end
        λ[i] = t
    end
    nothing
end
