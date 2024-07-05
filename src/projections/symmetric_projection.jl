struct SymmetricProjection{DT} <: ProjectionMethod 
    RU::Vector{DT}
    RG::Vector{DT}
    
    function SymmetricProjection(R∞ = 1)
        DT, RU, RG = _projection_weights([1//2,1//2], [1//2,1//2], R∞)
        new{DT}(RU, RG)
    end
end

SymmetricProjection(method::GeometricMethod) = ProjectedMethod(SymmetricProjection(), method)
SymmetricProjection(method::Union{RKMethod,PRKMethod,VPRKMethod}) = ProjectedMethod(SymmetricProjection(tableau(method).R∞), method)

const SymmetricProjectionIntegrator{PT} = ProjectionIntegrator{<:ProjectedMethod{<:SymmetricProjection, <:GeometricMethod}, PT} where {PT <: AbstractProblem}

function Cache{ST}(problem::EquationProblem, method::ProjectedMethod{<:SymmetricProjection}; kwargs...) where {ST}
    ProjectionCache{ST}(problem, method; kwargs...)
end

@inline CacheType(ST, problem::EquationProblem, method::ProjectedMethod{<:SymmetricProjection}) = 
    ProjectionCache{ST, timetype(problem), typeof(problem), ndims(problem), nconstraints(problem), solversize(problem, parent(method))}


default_solver(::ProjectedMethod{<:SymmetricProjection}) = Newton()
default_iguess(::ProjectedMethod{<:SymmetricProjection}) = HermiteExtrapolation()


# function Base.show(io::IO, int::ProjectedMethod{<:SymmetricProjection})
#     print(io, "\nProjection method with:\n")
#     print(io, "   Timestep: $(timestep(int))\n")
#     print(io, "   Tableau:  $(description(tableau(int)))\n")
#     print(io, "   $(string(tableau(int)))")
#     # print(io, reference(int.params.tab))
# end


function split_nlsolution(x::AbstractVector, int::SymmetricProjectionIntegrator)
    D = ndims(int)
    M = nconstraints(int)
    N = solversize(problem(int), parent(method(int)))

    x̄ = @view x[1:N]
    x̃ = @view x[N+1:N+D+M]

    return (x̄, x̃)
end


function initial_guess!(sol, history, params, int::SymmetricProjectionIntegrator)
    # compute initial guess for parent method
    initial_guess!(sol, history, params, subint(int))

    # copy initial guess for parent method to common solution vector
    cache(int).x̄ .= nlsolution(subint(int))

    # compute initial guess for projected solution
    soltmp = (
        t = (sol.t + history.t[1]) / 2,
        q = cache(int).q̃,
        p = cache(int).p̃,
        v = cache(int).ṽ,
        f = cache(int).f̃,
    )
    solutionstep!(soltmp, history, problem(int), iguess(int))
    # TODO: Fix this!

    # copy initial guess for projected solution to common solution vector
    cache(int).x̃[1:ndims(int)] .= cache(int).q̃

    # set initial guess for Lagrange multiplier to zero
    cache(int).x̃[ndims(int)+1:end] .= 0
end

function components!(x::AbstractVector{ST}, sol, params, int::SymmetricProjectionIntegrator{<:DAEProblem}) where {ST}
    # get cache for internal stages
    local C = cache(int, ST)

    # copy x to q
    for k in eachindex(C.q̃)
        C.q̄[k] = x[k]
    end

    # copy x to λ
    for k in eachindex(C.λ)
        C.λ[k] = x[ndims(int)+k]
    end

    # compute u=λ and g=∇ϑ(q)⋅λ
    equations(int).u(C.U[1], sol.t - timestep(int), sol.q, C.λ, params)
    equations(int).u(C.U[2], sol.t, C.q̄, C.λ, params)

    C.U[1] .*= projection(method(int)).RU[1]
    C.U[2] .*= projection(method(int)).RU[2]
end


function components!(x::AbstractVector{ST}, sol, params, int::SymmetricProjectionIntegrator{<:Union{IODEProblem,LODEProblem}}) where {ST}
    # get cache for internal stages
    local C = cache(int, ST)

    # copy x to q
    for k in eachindex(C.q̃)
        C.q̄[k] = x[k]
        C.ṽ[k] = (C.q̄[k] - sol.q[k]) / timestep(int)
    end

    # copy x to λ
    for k in eachindex(C.λ)
        C.λ[k] = x[ndims(int)+k]
    end

    # compute u = λ
    C.U[1] .= projection(method(int)).RU[1] .* C.λ
    C.U[2] .= projection(method(int)).RU[2] .* C.λ

    # compute g = ∇ϑ(q)⋅λ
    equations(int).g(C.G[1], sol.t - timestep(int), sol.q, C.ṽ, C.λ, params)
    equations(int).g(C.G[2], sol.t, C.q̄, C.ṽ, C.λ, params)

    C.G[1] .*= projection(method(int)).RG[1]
    C.G[2] .*= projection(method(int)).RG[2]
end


function residual!(b::AbstractVector{ST}, sol, params, int::SymmetricProjectionIntegrator{<:DAEProblem}) where {ST}
    # get cache for internal stages
    local C = cache(int, ST)

    # compute b = q̄ - q
    for k in 1:ndims(int)
        b[k] = C.q̄[k] - sol.q[k]
    end

    # compute ϕ = ϕ(q)
    equations(int).ϕ(C.ϕ, sol.t, sol.q, params)

    # compute b = ϕ(q) or b = ϕ(q,p) or b = ϕ(...)
    for k in 1:nconstraints(int)
        b[ndims(int)+k] = C.ϕ[k]
    end
end


function residual!(b::AbstractVector{ST}, sol, params, int::SymmetricProjectionIntegrator{<:Union{IODEProblem,LODEProblem}}) where {ST}
    # get cache for internal stages
    local C = cache(int, ST)

    # compute b = q̄ - q
    for k in 1:ndims(int)
        b[k] = C.q̄[k] - sol.q[k]
    end

    # compute ϑ(q)
    equations(int).ϑ(C.ϑ, sol.t, sol.q, sol.v, params)

    # compute b = ϕ(q) or b = ϕ(q,p) or b = ϕ(...)
    for k in 1:nconstraints(int)
        b[ndims(int)+k] = C.ϑ[k] - sol.p[k]
    end
end


# Compute stages of variational partitioned Runge-Kutta methods.
function residual!(b::AbstractVector{ST}, x::AbstractVector{ST}, sol, params, int::SymmetricProjectionIntegrator) where {ST}
    # check that x and b are compatible
    @assert axes(x) == axes(b)

    # copy previous solution from sol to cache
    # and assign cache to local variable for convenience
    local C = reset!(cache(int, ST), sol)
    local soltemp = current(C)

    # split x and b
    x̄, x̃ = split_nlsolution(x, int)
    b̄, b̃ = split_nlsolution(b, int)

    # compute projection variables
    components!(x̃, soltemp, params, int)

    # compute initial projection (perturbation)
    project!(soltemp, C.U[1], C.G[1], int)

    # call components method of parent integrator
    components!(x̄, soltemp, params, subint(int))

    # compute residual of parent method
    residual!(b̄, soltemp, params, subint(int))

    # update solution with vectorfield of parent integrator
    update!(soltemp, params, subint(int), ST)

    # compute final projection (perturbation)
    project!(soltemp, C.U[2], C.G[2], int)

    # compute residual of projection method
    residual!(b̃, soltemp, params, int)
end


function update!(sol, params, x::AbstractVector{ST}, int::SymmetricProjectionIntegrator) where {ST}
    # split x and b
    x̄, x̃ = split_nlsolution(x, int)

    # assign cache to local variable for convenience
    local C = cache(int, ST)

    # compute projection variables
    components!(x̃, sol, params, int)

    # compute initial projection (perturbation)
    project!(sol, C.U[1], C.G[1], int)

    # compute update of parent integrator
    update!(sol, params, x̄, subint(int))

    # compute final projection (perturbation)
    project!(sol, C.U[2], C.G[2], int)
end


function integrate_step!(sol, history, params, int::SymmetricProjectionIntegrator)
    # call nonlinear solver for projection
    solve!(nlsolution(int), (b,x) -> residual!(b, x, sol, params, int), solver(int))

    # check_jacobian(solver(int))
    # print_jacobian(solver(int))

    # print solver status
    # println(status(solver(int)))

    # check if solution contains NaNs or error bounds are violated
    # println(meets_stopping_criteria(status(solver(int))))

    # copy solver status
    # get_solver_status!(solver(int), solstep(int).internal[:solver])

    # update solution step
    update!(sol, params, nlsolution(int), int)
end
