@doc raw"""
Discrete Euler-Lagrange Method.

We consider discrete Euler-Lagrange equations of the form
```math
D_1 L_d (q_{n}, q_{n+1}) + D_2 L_d (q_{n-1}, q_{n}) = 0,
```
where $q_{n}$ approximates the solution $q(t_{n})$ and $L_d$ is any discrete Lagrangian.
"""
struct DiscreteEulerLagrange <: DELEMethod end

isexplicit(method::DiscreteEulerLagrange) = false
isimplicit(method::DiscreteEulerLagrange) = true
issymmetric(method::DiscreteEulerLagrange) = false
issymplectic(method::DiscreteEulerLagrange) = true


@doc raw"""
Discrete Euler-Lagrange integrator cache.
"""
struct DiscreteEulerLagrangeCache{DT,D} <: DELEIntegratorCache{DT,D}
    x::Vector{DT}
    q::Vector{DT}
    D1Ld::Vector{DT}
    D2Ld::Vector{DT}

    q̃::Vector{DT}
    ṽ::Vector{DT}
    θ̃::Vector{DT}
    f̃::Vector{DT}

    function DiscreteEulerLagrangeCache{DT,D}() where {DT,D}
        x = zeros(DT, D)
        q = zeros(DT, D)
        D1Ld = zeros(DT, D)
        D2Ld = zeros(DT, D)
        new(x, q, D1Ld, D2Ld,
            zeros(DT, D), zeros(DT, D), zeros(DT, D), zeros(DT, D))
    end
end

reset!(::DELEIntegratorCache, t₀, q₀, t₁, q₁) = nothing

nlsolution(cache::DiscreteEulerLagrangeCache) = cache.x

function Cache{ST}(problem::AbstractProblem, method::DiscreteEulerLagrange; kwargs...) where {ST}
    DiscreteEulerLagrangeCache{ST,ndims(problem)}(; kwargs...)
end

@inline CacheType(ST, problem::AbstractProblem, method::DiscreteEulerLagrange) = DiscreteEulerLagrangeCache{ST,ndims(problem)}


solversize(problem::AbstractProblemDELE, ::DiscreteEulerLagrange) = ndims(problem)

default_solver(::DiscreteEulerLagrange) = Newton()
default_iguess(::DiscreteEulerLagrange) = HermiteExtrapolation()


function initial_guess!(sol, history, params, int::GeometricIntegrator{<:DiscreteEulerLagrange})
    # compute initial guess for solution q(n+1)
    # soltmp = (
    #     t=sol.t,
    #     q=cache(int).q̃,
    #     p=cache(int).θ̃,
    #     v=cache(int).ṽ,
    #     f=cache(int).f̃,
    # )
    # solutionstep!(soltmp, history, problem(int), iguess(int))

    # copy initial guess to solution vector
    # nlsolution(int) .= cache(int).q̃
    nlsolution(int) .= sol.q
end

function components!(x::AbstractVector{ST}, sol, history, params, int::GeometricIntegrator{<:DiscreteEulerLagrange}) where {ST}
    local q = cache(int, ST).q
    local D1Ld = cache(int, ST).D1Ld
    local D2Ld = cache(int, ST).D2Ld

    # copy x to q
    q .= x

    # compute D1Ld(t_n, t_n+1, q_n, q_n+1)
    equations(int).D1Ld(D1Ld, history.t[1], sol.t, history.q[1], q, params)

    # compute D2Ld(t_n-1, t_n, q_n-1, q_n)
    equations(int).D2Ld(D2Ld, history.t[2], history.t[1], history.q[2], history.q[1], params)
end


function residual!(b::AbstractVector{ST}, int::GeometricIntegrator{<:DiscreteEulerLagrange}) where {ST}
    # compute b = D1Ld(t_n, t_n+1, q_n, q_n+1) + D2Ld(t_n-1, t_n, q_n-1, q_n)
    b .= cache(int, ST).D1Ld .+ cache(int, ST).D2Ld
end


# Compute stages of discrete Euler-Lagrange methods.
function residual!(b::AbstractVector{ST}, x::AbstractVector{ST}, sol, history, params, int::GeometricIntegrator{<:DiscreteEulerLagrange}) where {ST}
    @assert axes(x) == axes(b)

    # copy previous solution from solstep to cache
    reset!(cache(int, ST), sol..., history.t[1], history.q[1])

    # compute stages from nonlinear solver solution x
    components!(x, sol, history, params, int)

    # compute residual vector
    residual!(b, int)
end


function update!(sol, history, params, x::AbstractVector{DT}, int::GeometricIntegrator{<:DiscreteEulerLagrange}) where {DT}
    # copy previous solution from solstep to cache
    reset!(cache(int, DT), sol..., history.t[1], history.q[1])

    # compute vector field at internal stages
    components!(x, sol, history, params, int)

    # compute final update
    sol.q .= cache(int, DT).q
end


function integrate_step!(sol, history, params, int::GeometricIntegrator{<:DiscreteEulerLagrange,<:AbstractProblemDELE})
    # call nonlinear solver
    solve!(nlsolution(int), solver(int), (sol, history, params, int))

    # print solver status
    # println(status(solver))

    # check if solution contains NaNs or error bounds are violated
    # println(meets_stopping_criteria(status(solver)))

    # compute final update
    update!(sol, history, params, nlsolution(int), int)
end
