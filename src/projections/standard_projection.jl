struct StandardProjection{DT} <: ProjectionMethod 
    RU::Vector{DT}
    RG::Vector{DT}
    
    function StandardProjection(RU, RG, R∞ = 1)
        DT, RU, RG = _projection_weights(RU, RG, R∞)
        new{DT}(RU, RG)
    end
end

PostProjection(method::GeometricMethod) = ProjectedMethod(StandardProjection([0,1], [0,1]), method)
SymplecticProjection(method::Union{PRKMethod,VPRKMethod}) = ProjectedMethod(StandardProjection([0,1], [0,1], tableau(method).R∞), method)
VariationalProjectionOnP(method::GeometricMethod) = ProjectedMethod(StandardProjection([0,1], [1,0]), method)
VariationalProjectionOnQ(method::GeometricMethod) = ProjectedMethod(StandardProjection([1,0], [0,1]), method)

# description(::PostProjection) = "Post projection"
# description(::SymplecticProjection) = "Symplectic Projection"
# description(::VariationalProjectionOnP) = @doc raw"Variational projection on $(q_{n}, p_{n+1})$"
# description(::VariationalProjectionOnQ) = @doc raw"Variational projection on $(p_{n}, q_{n+1})$"

const StandardProjectionIntegrator{PT} = ProjectionIntegrator{<:ProjectedMethod{<:StandardProjection, <:GeometricMethod}, PT} where {PT <: AbstractProblem}

function Cache{ST}(problem::EquationProblem, method::ProjectedMethod{<:StandardProjection}; kwargs...) where {ST}
    ProjectionCache{ST}(problem, method; kwargs...)
end

@inline CacheType(ST, problem::EquationProblem, method::ProjectedMethod{<:StandardProjection}) =
    ProjectionCache{ST, timetype(problem), typeof(problem), ndims(problem), nconstraints(problem), solversize(problem, parent(method))}


default_solver(::ProjectedMethod{<:StandardProjection}) = Newton()


function split_nlsolution(x::AbstractVector, int::StandardProjectionIntegrator)
    D = ndims(int)
    M = nconstraints(int)
    N = solversize(problem(int), parent(method(int)))

    x̄ = @view x[1:N]
    x̃ = @view x[N+1:N+D+M]

    return (x̄, x̃)
end


function initsolver(::NewtonMethod, config::Options, ::ProjectedMethod{<:StandardProjection}, caches::CacheDict; kwargs...)
    x̄, x̃ = split_nlsolution(cache(caches))
    NewtonSolver(zero(x̃), zero(x̃); config = config, kwargs...)
end


# function Base.show(io::IO, int::ProjectedMethod{<:StandardProjection})
#     print(io, "\nProjection method with:\n")
#     print(io, "   Timestep: $(timestep(int))\n")
#     print(io, "   Tableau:  $(description(tableau(int)))\n")
#     print(io, "   $(string(tableau(int)))")
#     # print(io, reference(int.params.tab))
# end


function initial_guess!(sol, history, params, int::StandardProjectionIntegrator)
    # compute initial guess for parent method
    initial_guess!(sol, history, params, subint(int))

    # set initial guess for Lagrange multiplier to zero
    cache(int).x̃[ndims(int)+1:end] .= 0
end


function components!(x::AbstractVector{ST}, sol, params, int::StandardProjectionIntegrator{<:DAEProblem}) where {ST}
    # get cache for internal stages
    local C = cache(int, ST)

    # copy x to q
    for k in eachindex(C.q)
        C.q[k] = x[k]
    end

    # copy x to λ
    for k in eachindex(C.λ)
        C.λ[k] = x[ndims(int)+k]
    end

    # compute u = u(q,λ)
    equations(int).u(C.u, sol.t, C.q, C.λ, params)
    C.U[1] .= projection(method(int)).RU[1] .* C.u
    C.U[2] .= projection(method(int)).RU[2] .* C.u

    # compute ϕ = ϕ(q)
    equations(int).ϕ(C.ϕ, sol.t, C.q, params)
end


function components!(x::AbstractVector{ST}, sol, params, int::StandardProjectionIntegrator{<:Union{IODEProblem,LODEProblem}}) where {ST}
    # get cache for internal stages
    local C = cache(int, ST)

    # copy x to q
    for k in eachindex(C.q)
        C.q[k] = x[k]
    end

    # copy x to λ
    for k in eachindex(C.λ)
        C.λ[k] = x[ndims(int)+k]
    end

    # compute u = λ
    C.U[1] .= projection(method(int)).RU[1] .* C.λ
    C.U[2] .= projection(method(int)).RU[2] .* C.λ

    # compute g = ∇ϑ(q)⋅λ
    equations(int).g(C.g, sol.t, C.q, C.v, C.λ, params)
    # TODO: This cache(int, ST).v does not necessary have a proper value
    # (important for the non-degenerate case)
    C.G[1] .= projection(method(int)).RG[1] .* C.g
    C.G[2] .= projection(method(int)).RG[2] .* C.g

    # project p
    C.p .= sol.p .+ timestep(int) .* C.G[2]

    # compute ϕ = ϑ(q) - p
    equations(int).ϑ(C.ϕ, sol.t, C.q, C.v, params)
    C.ϕ .-= C.p
end

function residual!(b::AbstractVector{ST}, sol, params, int::StandardProjectionIntegrator) where {ST}
    # get cache for internal stages
    local C = cache(int, ST)

    # compute b = q̄ - q - Δt * U
    for k in 1:ndims(int)
        b[k] = C.q[k] - sol.q[k] - timestep(int) * C.U[2][k]
    end

    # compute b = ϕ(q) or b = ϕ(q,p) or b = ϕ(...)
    for k in 1:nconstraints(int)
        b[ndims(int)+k] = C.ϕ[k]
    end
end

function residual!(b::AbstractVector{ST}, x::AbstractVector{ST}, sol, params, int::StandardProjectionIntegrator) where {ST}
    # check that x and b are compatible
    @assert axes(x) == axes(b)

    # compute stages
    components!(x, sol, params, int)

    # compute residual of projection method
    residual!(b, sol, params, int)
end


# function update!(sol, params, x::AbstractVector{ST}, int::StandardProjectionIntegrator) where {ST}
#     # add perturbation for next time step to solution
#     # (same vector field as previous time step)
#     project!(sol, cache(int).U[1], cache(int).G[1], int)

#     # compute update of parent integrator
#     update!(sol, params, x̄, subint(int))

#     # add projection to solution
#     project!(sol, cache(int).U[2], cache(int).G[2], int)
# end


function integrate_step!(sol, history, params, int::StandardProjectionIntegrator)
    # add perturbation for next time step to solution
    # (same vector field as previous time step)
    project!(sol, cache(int).U[1], cache(int).G[1], int)

    # integrate one step with parent method
    integrate_step!(sol, history, params, subint(int))

    # copy initial guess for projected solution to common solution vector
    cache(int).x̃[1:ndims(int)] .= sol.q

    # call nonlinear solver for projection
    x̄, x̃ = split_nlsolution(nlsolution(int), int)
    solve!(x̃, (b,x) -> residual!(b, x, sol, params, int), solver(int))

    # check_jacobian(solver(int))
    # print_jacobian(solver(int))

    # print solver status
    # println(status(solver(int)))

    # check if solution contains NaNs or error bounds are violated
    # println(meets_stopping_criteria(status(solver(int))))

    # TODO: copy λ/U/G to internal variables

    # add projection to solution
    components!(x̃, sol, params, int)
    project!(sol, cache(int).U[2], cache(int).G[2], int)

    # println()
    # println("λ    = $(cache(int).λ)")
    # println("U[1] = $(cache(int).U[1])")
    # println("U[2] = $(cache(int).U[2])")
    # println("G[1] = $(cache(int).G[1])")
    # println("G[2] = $(cache(int).G[2])")
    # println()

    # # update solution step
    # update!(sol, params, nlsolution(int), int)
end
