
# general helper functions for integrators.

_default_solver() = SimpleSolvers.get_config(:nls_solver)

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
    solver(x, zero(x), f!)
end

function create_nonlinear_solver(DT, N, params, caches, i::Int, solver::Type{<:NonlinearSolver}=_default_solver(), F::Function=function_stages!)
    # create solution vector for nonlinear solver
    x = zeros(DT, N)

    # create wrapper function f(x,b) that calls `F(x, b, params)`
    # with the appropriate `params`
    f! = (b,x) -> F(x, b, params, caches, i)

    # create nonlinear solver with solver type obtained from config dictionary
    solver(x, zero(x), f!)
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
    solver(x, zero(x), f!; J! = j!)
end


"""
Create a vector of S solution vectors of type DT to store the solution of S
internal stages for a problem with `D` dimensions.
"""
function create_internal_stage_vector(DT, D, S)
    [zeros(DT,D) for i in 1:S]
end


"""
Create a vector of S solution matrices of type DT to store the solution of S
internal stages for a problem with `DxM` dimensions.
"""
function create_internal_stage_matrix(DT, D, M, S)
    [zeros(DT,D,M) for i in 1:S]
end


"""
Create a vector of S+1 solution vectors of type DT to store the solution of S
internal stages and the solution of the previous timestep for a problem with `D`
    dimensions.
"""
function create_internal_stage_vector_with_zero(DT, D, S)
    a = OffsetArray{Vector{DT}}(undef, 0:S)

    for i in 0:S
        a[i] = zeros(DT,D)
    end

    return a
end


"""
Create a vector of (S,M+1) solution vectors of type DT to store the solution of S
internal stages and M random processes for a problem with `D` dimensions.
"""
function create_internal_stage_vector_with_zero(DT, D, M, S)
    a = OffsetArray{Vector{DT}}(undef, 0:M, 1:S)

    for i in 0:M
        for j in 1:S
            a[i,j] = zeros(DT,D)
        end
    end

    return a
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


function GeometricBase.cut_periodic_solution!(x::SolutionVector{T}, periodicity::AbstractVector{T}) where {T}
    @assert length(x) == length(periodicity)

    for k in eachindex(x, periodicity)
        if periodicity[k] ≠ 0
            while x[k] < 0
                x[k] += periodicity[k]
            end
            while x[k] ≥ periodicity[k]
                x[k] -= periodicity[k]
            end
        end
    end
end

function GeometricBase.cut_periodic_solution!(x::SolutionVector{T}, periodicity::AbstractVector{T}, shift::AbstractVector{T}) where {T}
    @assert length(x) == length(periodicity)
    shift .= 0
    for k in eachindex(x, periodicity, shift)
        if periodicity[k] ≠ 0
            while x[k] + shift[k] < 0
                shift[k] += periodicity[k]
            end
            while x[k] + shift[k] ≥ periodicity[k]
                shift[k] -= periodicity[k]
            end
        end
    end
end
