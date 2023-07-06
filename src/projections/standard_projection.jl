
struct StandardProjectionCache{DT,D,M} <: IODEIntegratorCache{DT,D}
    x::Vector{DT}
    q::Vector{DT}
    p::Vector{DT}
    λ::Vector{DT}
    ϕ::Vector{DT}
    u::Vector{DT}
    g::Vector{DT}
    U::Vector{Vector{DT}}
    G::Vector{Vector{DT}}

    function StandardProjectionCache{DT,D,M}() where {DT,D,M}
        x = zeros(DT, D+M)
        q = zeros(DT, D)
        p = zeros(DT, D)
        λ = zeros(DT, M)
        ϕ = zeros(DT, M)
        g = zeros(DT, D)
        u = zeros(DT, D)
        U = [zeros(DT, D), zeros(DT, D)]
        G = [zeros(DT, D), zeros(DT, D)]
        new(x, q, p, λ, ϕ, u, g, U, G)
    end
end

nlsolution(cache::StandardProjectionCache) = cache.x


function Cache{ST}(problem::GeometricProblem, ::ProjectedMethod{<:StandardProjection}; kwargs...) where {ST}
    StandardProjectionCache{ST, ndims(problem), nconstraints(problem)}(; kwargs...)
end

@inline CacheType(ST, problem::GeometricProblem, ::ProjectedMethod{<:StandardProjection}) = StandardProjectionCache{ST, ndims(problem), nconstraints(problem)}


default_solver(::ProjectedMethod{<:StandardProjection}) = Newton()


const IntegratorStandardProjection{DT,TT} = Integrator{<:GeometricProblem{DT,TT}, <:ProjectedMethod{<:StandardProjection}}


# function Base.show(io::IO, int::ProjectedMethod{<:StandardProjection})
#     print(io, "\nProjection method with:\n")
#     print(io, "   Timestep: $(timestep(int))\n")
#     print(io, "   Tableau:  $(description(tableau(int)))\n")
#     print(io, "   $(string(tableau(int)))")
#     # print(io, reference(int.params.tab))
# end


function initial_guess!(int::IntegratorStandardProjection)
    # initial_guess!(subint(int))
    cache(int).x[1:ndims(int)] .= solstep(int).q
    cache(int).x[ndims(int)+1:end] .= 0
end


function initial_guess!(
    solstep::SolutionStep,
    ::GeometricProblem,
    ::ProjectedMethod{<:StandardProjection},
    caches::CacheDict,
    ::NonlinearSolver,
    ::NoInitialGuess)
    
    caches[datatype(solstep)].x .= 0
end


function components!(
    x::AbstractVector{DT},
    solstep::SolutionStep, 
    problem::DAEProblem,
    method::ProjectedMethod{<:StandardProjection}, 
    caches::CacheDict) where {DT}

    q = caches[DT].q
    λ = caches[DT].λ
    ϕ = caches[DT].ϕ
    u = caches[DT].u
    U = caches[DT].U

    # copy x to q
    for k in eachindex(q)
        q[k] = x[k]
    end

    # copy x to λ
    for k in eachindex(λ)
        λ[k] = x[ndims(problem)+k]
    end

    # compute u = u(q,λ)
    functions(problem).u(u, solstep.t, q, λ)
    U[1] .= projection(method).RU[1] .* u
    U[2] .= projection(method).RU[2] .* u

    # compute ϕ = ϕ(q)
    functions(problem).ϕ(ϕ, solstep.t, q)
end


function components!(
    x::AbstractVector{ST},
    solstep::SolutionStep, 
    problem::Union{IODEProblem,LODEProblem},
    method::ProjectedMethod{<:StandardProjection}, 
    caches::CacheDict) where {ST}

    q = caches[ST].q
    p = caches[ST].p
    λ = caches[ST].λ
    ϕ = caches[ST].ϕ
    g = caches[ST].g
    U = caches[ST].U
    G = caches[ST].G

    # copy x to q
    for k in eachindex(q)
        q[k] = x[k]
    end

    # copy x to λ
    for k in eachindex(λ)
        λ[k] = x[ndims(problem)+k]
    end

    # compute u = λ
    U[1] .= projection(method).RU[1] .* λ
    U[2] .= projection(method).RU[2] .* λ

    # compute g = ∇ϑ(q)⋅λ
    functions(problem).g(g, solstep.t, q, solstep.v, λ)
    G[1] .= projection(method).RG[1] .* g
    G[2] .= projection(method).RG[2] .* g

    # project p
    p .= solstep.p .+ timestep(problem) .* G[2]

    # compute ϕ = ϑ(q) - p
    functions(problem).ϑ(ϕ, solstep.t, q, solstep.v)
    ϕ .-= p
end


# Compute stages of variational partitioned Runge-Kutta methods.
function residual!(
    b::Vector{ST},
    x::Vector{ST},
    solstep::SolutionStep, 
    problem::GeometricProblem,
    method::ProjectedMethod{<:StandardProjection}, 
    caches::CacheDict) where {ST}

    @assert axes(x) == axes(b)

    # compute stages
    components!(x, solstep, problem, method, caches)

    # compute b = q̄ - q - Δt * U
    for k in 1:ndims(problem)
        b[k] = caches[ST].q[k] - solstep.q[k] - timestep(problem) * caches[ST].U[2][k]
    end

    # compute b = ϕ(q) or b = ϕ(q,p) or b = ϕ(...)
    for k in 1:nconstraints(problem)
        b[ndims(problem)+k] = caches[ST].ϕ[k]
    end
end


function preprojection!(
    solstep::SolutionStep{DT},
    problem::GeometricProblem,
    method::StandardProjection, 
    caches::CacheDict) where {DT}

    # TODO: Do not use U and G from cache here but from solstep!
    project!(solstep, problem, method, caches[DT].U[1], caches[DT].G[1], caches[DT].λ)
end

function postprojection!(
    solstep::SolutionStep{DT},
    problem::GeometricProblem,
    method::StandardProjection, 
    caches::CacheDict) where {DT}

    # TODO: Do not use U and G from cache here but from solstep!
    project!(solstep, problem, method, caches[DT].U[2], caches[DT].G[2], caches[DT].λ)
end


function integrate_step!(int::Integrator{<:GeometricProblem, <:ProjectedMethod{<:StandardProjection}})
    # add perturbation for next time step to solution
    # (same vector field as previous time step)
    preprojection!(solstep(int), problem(int), projection(method(int)), caches(int))

    # compute initial guess for parent method
    initial_guess!(subint(int))

    # integrate one step with parent method
    integrate_step!(subint(int))

    # compute initial guess for parent method
    initial_guess!(int)

    # call nonlinear solver for projection
    solve!(nlsolution(int), (b,x) -> residual!(b, x, solstep(int), problem(int), method(int), caches(int)), solver(int))

    # check_jacobian(solver(int))
    # print_jacobian(solver(int))

    # print solver status
    # println(status(solver(int)))

    # check if solution contains NaNs or error bounds are violated
    # println(meets_stopping_criteria(status(solver(int))))

    # TODO: copy λ/U/G to solstep and use λ/U/G from solstep in pre/post projection
    # currently not possible as we are using a SolutionPODE

    # add projection to solution
    postprojection!(solstep(int), problem(int), projection(method(int)), caches(int))

    # copy solver status
    # get_solver_status!(solver(int), solstep(int).internal[:solver])
end
