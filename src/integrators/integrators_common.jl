
# general helper functions for integrators.

"""
Create nonlinear solver object for a system of `N` equations with data type `DT`.
The function ``f(x)=0`` to be solved for is determined by a julia function
`function_stages!(x, b, params)`, where `x` is the current solution and `b` is
the output vector, s.th. ``b = f(x)``. `params` are a set of parameters depending
on the equation and integrator that is used.
The solver type is obtained from the config dictionary (`:nls_solver`).
"""
function create_nonlinear_solver(DT, N, params; F=function_stages!)
    # create solution vector for nonlinear solver
    x = zeros(DT, N)

    # create wrapper function f!(x,b) that calls `function_stages!(x, b, params)`
    # with the appropriate `params`
    f! = (x,b) -> F(x, b, params)

    # create nonlinear solver with solver type obtained from config dictionary
    s = get_config(:nls_solver)(x, f!)
end


function create_nonlinear_solver_with_jacobian(DT, N, params)
    # create solution vector for nonlinear solver
    x = zeros(DT, N)

    # create wrapper function f!(x,b) that calls `function_stages!(x, b, params)`
    # with the appropriate `params`
    f! = (x,b) -> function_stages!(x, b, params)

    # create wrapper function j!(x,df) that calls `jacobian!(x, df, params)`
    # with the appropriate `params`
    cache = IntegratorCache(params)
    j! = (x,df) -> jacobian!(x, df, cache, params)

    # create nonlinear solver with solver type obtained from config dictionary
    s = get_config(:nls_solver)(x, f!; J! = j!)
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
function create_internal_stage_vector(DT, D, M, S)
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


function update_solution!(x::Vector{T}, xₑᵣᵣ::Vector{T}, ẋ::Matrix{T}, b::Vector{T}, Δt::T) where {T}
    @assert length(x) == length(xₑᵣᵣ)
    @assert length(x) == size(ẋ, 1)
    @assert length(b) == size(ẋ, 2)

    for k in axes(ẋ, 1)
        for i in axes(ẋ, 2)
            x[k], xₑᵣᵣ[k] = compensated_summation(Δt * b[i] * ẋ[k,i], x[k], xₑᵣᵣ[k])
        end
    end
end

function update_solution!(x::Vector{T}, xₑᵣᵣ::Vector{T}, ẋ::Vector{Vector{T}}, b::Vector{T}, Δt::T) where {T}
    @assert length(b) == length(ẋ)
    @assert length(x) == length(ẋ[1])
    @assert length(x) == length(xₑᵣᵣ)

    for i in eachindex(ẋ)
        for k in eachindex(ẋ[i])
            x[k], xₑᵣᵣ[k] = compensated_summation(Δt * b[i] * ẋ[i][k], x[k], xₑᵣᵣ[k])
        end
    end
end

function update_solution!(x::SolutionVector{T}, ẋ::Matrix{T}, b::Vector{T}, Δt::T) where {T}
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


function update_solution!(x::SolutionVector{T}, ẋ::Vector{Vector{T}}, b::Vector{T}, Δt::T) where {T}
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

function update_solution!(x::Vector{T}, xₑᵣᵣ::Vector{T}, ẋ::Union{Matrix{T},Vector{Vector{T}}}, b::Vector{T}, b̂::Vector, Δt::T) where {T}
    update_solution!(x, xₑᵣᵣ, ẋ, b, Δt)
    update_solution!(x, xₑᵣᵣ, ẋ, b̂, Δt)
end

function update_solution!(x::SolutionVector{T}, ẋ::Union{Matrix{T},Vector{Vector{T}}}, b::Vector{T}, b̂::Vector, Δt::T) where {T}
    update_solution!(x, ẋ, b, Δt)
    update_solution!(x, ẋ, b̂, Δt)
end

function update_multiplier!(λ::SolutionVector{T}, Λ::Vector{Vector{T}}, b::Vector{T}) where {T}
    @assert length(λ) == size(Λ, 1)
    local t::T
    @inbounds for i in eachindex(λ)
        t = zero(T)
        for j=eachindex(b)
            t += b[j] * Λ[j][i]
        end
        λ[i] = t
    end
    nothing
end


function CommonFunctions.cut_periodic_solution!(x::SolutionVector{T}, periodicity::Vector{T}) where {T}
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

function CommonFunctions.cut_periodic_solution!(x::SolutionVector{T}, periodicity::Vector{T}, shift::Vector{T}) where {T}
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
