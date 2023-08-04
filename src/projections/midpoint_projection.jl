
function Cache{ST}(problem::EquationProblem, method::ProjectedMethod{<:MidpointProjection}; kwargs...) where {ST}
    ProjectionCache{ST, ndims(problem), nconstraints(problem), solversize(problem, parent(method))}(; kwargs...)
end

@inline CacheType(ST, problem::EquationProblem, method::ProjectedMethod{<:MidpointProjection}) = 
    ProjectionCache{ST, ndims(problem), nconstraints(problem), solversize(problem, parent(method))}


default_solver(::ProjectedMethod{<:MidpointProjection}) = Newton()
default_iguess(::ProjectedMethod{<:MidpointProjection}) = HermiteExtrapolation()


const IntegratorMidpointProjection{DT,TT} = Integrator{<:EquationProblem{DT,TT}, <:ProjectedMethod{<:MidpointProjection}}

# TODO: Try to disable this, once everything works!
function initsolver(::NewtonMethod, ::ProjectedMethod{<:MidpointProjection}, caches::CacheDict; kwargs...)
    x = zero(nlsolution(caches))
    y = zero(nlsolution(caches))
    NewtonSolver(x, y; linesearch = Backtracking(), config = Options(min_iterations = 1, f_abstol = 2eps(eltype(nlsolution(caches)))), kwargs...)
end


# function Base.show(io::IO, int::ProjectedMethod{<:MidpointProjection})
#     print(io, "\nProjection method with:\n")
#     print(io, "   Timestep: $(timestep(int))\n")
#     print(io, "   Tableau:  $(description(tableau(int)))\n")
#     print(io, "   $(string(tableau(int)))")
#     # print(io, reference(int.params.tab))
# end


function split_nlsolution(x::AbstractVector, int::IntegratorMidpointProjection)
    D = ndims(int)
    M = nconstraints(int)
    N = solversize(problem(int), parent(method(int)))

    x̄ = @view x[1:N]
    x̃ = @view x[N+1:N+D+M]

    return (x̄, x̃)
end


function initial_guess!(int::IntegratorMidpointProjection)
    # compute initial guess for parent method
    initial_guess!(subint(int))

    # copy initial guess for parent method to common solution vector
    cache(int).x̄ .= nlsolution(subint(int))

    # compute initial guess for projected solution
    initialguess!(solstep(int).t̄ + timestep(int)/2, cache(int).q̃, cache(int).v, solstep(int), problem(int), iguess(int))

    # copy initial guess for projected solution to common solution vector
    cache(int).x̃[1:ndims(int)] .= cache(int).q̃

    # set initial guess for Lagrange multiplier to zero
    cache(int).x̃[ndims(int)+1:end] .= 0
end


function components!(
    x::AbstractVector{DT},
    solstep::SolutionStep, 
    problem::DAEProblem,
    method::ProjectedMethod{<:MidpointProjection}, 
    caches::CacheDict) where {DT}

    q̃ = caches[DT].q̃
    λ = caches[DT].λ
    u = caches[DT].u
    U = caches[DT].U

    # copy x to q
    for k in eachindex(q̃)
        q̃[k] = x[k]
    end

    # copy x to λ
    for k in eachindex(λ)
        λ[k] = x[ndims(problem)+k]
    end

    # compute u=λ and g=∇ϑ(q)⋅λ
    functions(problem).u(u, solstep.t, q̃, λ)
    U[1] .= projection(method).RU[1] .* u
    U[2] .= projection(method).RU[2] .* u
end


function components!(
    x::AbstractVector{ST},
    solstep::SolutionStep, 
    problem::Union{IODEProblem,LODEProblem},
    method::ProjectedMethod{<:MidpointProjection}, 
    caches::CacheDict) where {ST}

    q̃ = caches[ST].q̃
    λ = caches[ST].λ
    g = caches[ST].g
    U = caches[ST].U
    G = caches[ST].G

    # copy x to q
    for k in eachindex(q̃)
        q̃[k] = x[k]
    end

    # copy x to λ
    for k in eachindex(λ)
        λ[k] = x[ndims(problem)+k]
    end

    # compute u = λ
    U[1] .= projection(method).RU[1] .* λ
    U[2] .= projection(method).RU[2] .* λ

    # compute g = ∇ϑ(q)⋅λ
    functions(problem).g(g, solstep.t, q̃, solstep.v, λ)
    G[1] .= projection(method).RG[1] .* g
    G[2] .= projection(method).RG[2] .* g
end


function constraint!(solstep::SolutionStep, problem::DAEProblem, ::IntegratorMidpointProjection, cache::ProjectionCache)
    # compute ϕ = ϕ(q)
    functions(problem).ϕ(cache.ϕ, solstep.t, cache.q)
end


function constraint!(solstep::SolutionStep, problem::Union{IODEProblem,LODEProblem}, ::IntegratorMidpointProjection, cache::ProjectionCache)
    # compute ϕ = ϑ(q) - p
    functions(problem).ϑ(cache.ϕ, solstep.t, cache.q, solstep.v)
    cache.ϕ .-= cache.p
end


function components!(x::AbstractVector{ST}, int::IntegratorMidpointProjection) where {ST}
    # TODO: Further generalise for non-RK methods
    # Need to implement update_vector! for integrators

    # split x and b
    x̄, x̃ = split_nlsolution(x, int)

    # compute stages
    components!(x̃, solstep(int), problem(int), method(int), caches(int))

    # compute initial projection (perturbation)
    project!(cache(int, ST), cache(int, ST).U[1], cache(int, ST).G[1], timestep(int))

    # copy projected solution to cache of subint
    reset!(cache(subint(int), ST), current(cache(int, ST), solstep(int))...)

    # call components method of parent integrator
    components!(x̄, subint(int))

    # update solution with vectorfield of parent integrator
    update!(cache(int, ST), cache(subint(int), ST), tableau(subint(int)), timestep(subint(int)))

    # compute final projection (perturbation)
    project!(cache(int, ST), cache(int, ST).U[2], cache(int, ST).G[2], timestep(int))
end


function residual!(b::AbstractVector{ST}, int::IntegratorMidpointProjection) where {ST}
    # compute b = q̃ - (q + q̄) / 2
    for k in 1:ndims(int)
        b[k] = cache(int, ST).q̃[k] - ( cache(int, ST).q[k] + solstep(int).q̄[k] ) / 2
    end

    # compute b = ϕ(q) or b = ϕ(q,p) or b = ϕ(...)
    for k in 1:nconstraints(int)
        b[ndims(int)+k] = cache(int, ST).ϕ[k]
    end
end


# Compute stages of variational partitioned Runge-Kutta methods.
function residual!(
    b::AbstractVector{ST},
    x::AbstractVector{ST},
    int::IntegratorMidpointProjection) where {ST}

    @assert axes(x) == axes(b)

    # copy previous solution from solstep to cache
    reset!(cache(int, ST), current(solstep(int))...)

    # update solstep from nonlinear solution vector
    components!(x, int)

    # update constraint
    constraint!(solstep(int), problem(int), int, cache(int, ST))

    # split b
    b̄, b̃ = split_nlsolution(b, int)

    # compute residual of parent method
    residual!(b̄, subint(int))

    # compute residual of projection method
    residual!(b̃, int)
end


function update!(x::AbstractVector{ST}, int::IntegratorMidpointProjection) where {ST}
    # split x and b
    x̄, x̃ = split_nlsolution(x, int)

    # copy previous solution from solstep to cache
    reset!(cache(int, ST), current(solstep(int))...)

    # compute stages
    components!(x̃, solstep(int), problem(int), method(int), caches(int))

    # compute initial projection (perturbation)
    project!(solstep(int), problem(int), method(int), cache(int, ST).U[1], cache(int, ST).G[1], cache(int, ST).λ)

    # compute update of parent integrator
    update!(x̄, subint(int))

    # compute final projection (perturbation)
    project!(solstep(int), problem(int), method(int), cache(int, ST).U[2], cache(int, ST).G[2], cache(int, ST).λ)
end


function integrate_step!(int::Integrator{<:EquationProblem, <:ProjectedMethod{<:MidpointProjection}})
    # call nonlinear solver for projection
    solve!(nlsolution(int), (b,x) -> residual!(b, x, int), solver(int))

    # check_jacobian(solver(int))
    # print_jacobian(solver(int))

    # print solver status
    # println(status(solver(int)))

    # check if solution contains NaNs or error bounds are violated
    # println(meets_stopping_criteria(status(solver(int))))

    # copy solver status
    # get_solver_status!(solver(int), solstep(int).internal[:solver])

    # update solution step
    update!(nlsolution(int), int)

    # update vector field for initial guess
    update_vector_fields!(solstep(int), problem(int))
end
