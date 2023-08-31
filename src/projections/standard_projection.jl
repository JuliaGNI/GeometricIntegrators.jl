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


function Cache{ST}(problem::EquationProblem, method::ProjectedMethod{<:StandardProjection}; kwargs...) where {ST}
    ProjectionCache{ST, ndims(problem), nconstraints(problem), solversize(problem, parent(method))}(; kwargs...)
end

@inline CacheType(ST, problem::EquationProblem, method::ProjectedMethod{<:StandardProjection}) = ProjectionCache{ST, ndims(problem), nconstraints(problem), solversize(problem, parent(method))}


default_solver(::ProjectedMethod{<:StandardProjection}) = Newton()


function initsolver(::NewtonMethod, ::ProjectedMethod{<:StandardProjection}, caches::CacheDict; kwargs...)
    x̄, x̃ = split_nlsolution(cache(caches))
    NewtonSolver(zero(x̃), zero(x̃); kwargs...)
end


const IntegratorStandardProjection{DT,TT} = ProjectionIntegrator{<:EquationProblem{DT,TT}, <:ProjectedMethod{<:StandardProjection}}


# function Base.show(io::IO, int::ProjectedMethod{<:StandardProjection})
#     print(io, "\nProjection method with:\n")
#     print(io, "   Timestep: $(timestep(int))\n")
#     print(io, "   Tableau:  $(description(tableau(int)))\n")
#     print(io, "   $(string(tableau(int)))")
#     # print(io, reference(int.params.tab))
# end


function split_nlsolution(x::AbstractVector, int::IntegratorStandardProjection)
    D = ndims(int)
    M = nconstraints(int)
    N = solversize(problem(int), parent(method(int)))

    x̄ = @view x[1:N]
    x̃ = @view x[N+1:N+D+M]

    return (x̄, x̃)
end


function initial_guess!(int::IntegratorStandardProjection)
    cache(int).x̃[1:ndims(int)] .= solstep(int).q
    cache(int).x̃[ndims(int)+1:end] .= 0
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
    b::AbstractVector{ST},
    x::AbstractVector{ST},
    solstep::SolutionStep, 
    problem::EquationProblem,
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
    problem::EquationProblem,
    method::StandardProjection, 
    caches::CacheDict) where {DT}

    # TODO: Do not use U and G from cache here but from solstep!
    project!(solstep, problem, method, caches[DT].U[1], caches[DT].G[1], caches[DT].λ)
end

function postprojection!(
    solstep::SolutionStep{DT},
    problem::EquationProblem,
    method::StandardProjection, 
    caches::CacheDict) where {DT}

    # TODO: Do not use U and G from cache here but from solstep!
    project!(solstep, problem, method, caches[DT].U[2], caches[DT].G[2], caches[DT].λ)
end


function integrate_step!(int::IntegratorStandardProjection)
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
    x̄, x̃ = split_nlsolution(nlsolution(int), int)
    solve!(x̃, (b,x) -> residual!(b, x, solstep(int), problem(int), method(int), caches(int)), solver(int))

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
