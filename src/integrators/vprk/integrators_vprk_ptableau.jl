
using NLsolve
using SimpleSolvers

"Parameters for right-hand side function of projected Gauss-Legendre Runge-Kutta methods."
mutable struct ParametersVPRKpTableau{DT, TT, D, S, ET <: NamedTuple} <: Parameters{DT,TT}
    equs::ET
    tab::CoefficientsPGLRK{TT}
    Δt::TT

    t::TT
    t̄::TT
    q̄::Vector{DT}
    p̄::Vector{DT}

    λ::Vector{DT}
    A::Matrix{DT}
    W::Matrix{DT}
    T::Matrix{DT}

    function ParametersVPRKpTableau{DT,D}(equs::ET, tab::CoefficientsPGLRK{TT}, Δt::TT) where {DT,TT,D,ET}
        S = tab.s
        new{DT, TT, D, S, ET}(equs, tab, Δt, zero(TT), zero(TT), zeros(DT,D), zeros(DT,D), zeros(DT,D), zeros(DT,S,S), zeros(DT,S,S), zeros(DT,S,S))
    end
end


function Integrators.IntegratorCache(params::ParametersVPRKpTableau{DT,TT,D,S}; kwargs...) where {DT,TT,D,S}
    IntegratorCacheVPRK{DT,D,S}(true; kwargs...)
end

function Integrators.IntegratorCache{ST}(params::ParametersVPRKpTableau{DT,TT,D,S}; kwargs...) where {ST,DT,TT,D,S}
    IntegratorCacheVPRK{ST,D,S}(true; kwargs...)
end

@inline Integrators.CacheType(ST, params::ParametersVPRKpTableau{DT,TT,D,S}) where {DT,TT,D,S} = IntegratorCacheVPRK{ST,D,S}


function update_params!(params::ParametersVPRKpTableau, sol::AtomicSolutionPODE)
    # set time for nonlinear solver and copy previous solution
    params.t̄  = sol.t
    params.t  = sol.t + params.Δt
    params.q̄ .= sol.q
    params.p̄ .= sol.p
end


function update_tableau!(params::ParametersVPRKpTableau{DT,TT,D,S}, λ::Vector) where {DT,TT,D,S}
    # copy λ to parameters
    params.λ .= λ

    # compute tableaus
    for k in eachindex(λ)
        params.W[S-k+1, S-k] = +λ[k]
        params.W[S-k, S-k+1] = -λ[k]
    end

    # params.A .= params.tab.P * params.W * params.tab.Q
    mul!(params.T, params.W, params.tab.Q)
    mul!(params.A, params.tab.P, params.T)
end


"""
Projected Variational Gauss-Legendre Runge-Kutta integrator.
"""
struct IntegratorVPRKpTableau{DT, TT, D, S,
                PT <: ParametersVPRKpTableau{DT,TT},
                ST <: NonlinearSolver{DT},
                IT <: InitialGuessIODE{TT}} <: AbstractIntegratorPRK{DT,TT}
    params::PT
    solver::ST
    iguess::IT
    caches::CacheDict{PT}

    function IntegratorVPRKpTableau(params::ParametersVPRKpTableau{DT,TT,D,S}, solver::ST, iguess::IT, caches) where {DT, TT, D, S, ST, IT}
        new{DT, TT, D, S, typeof(params), ST, IT}(params, solver, iguess, caches)
    end

    function IntegratorVPRKpTableau{DT,D}(equations::NamedTuple, tableau::CoefficientsPGLRK{TT}, Δt::TT) where {DT,TT,D}
        # get number of stages
        S = tableau.s

        # create params
        params = ParametersVPRKpTableau{DT,D}(equations, tableau, Δt)

        # create cache dict
        caches = CacheDict(params)

        # create solver
        solver  = create_nonlinear_solver(DT, D*S, params, caches)

        # create initial guess
        iguess = InitialGuessIODE(get_config(:ig_extrapolation), equations[:v̄], equations[:f̄], Δt)

        # create integrator
        IntegratorVPRKpTableau(params, solver, iguess, caches)
    end

    function IntegratorVPRKpTableau(equation::Union{IODE{DT},LODE{DT}}, tableau, Δt; kwargs...) where {DT}
        IntegratorVPRKpTableau{DT, ndims(equation)}(get_functions(equation), tableau, Δt; kwargs...)
    end
end


@inline nstages(integrator::IntegratorVPRKpTableau{DT,TT,D,S}) where {DT,TT,D,S} = S
@inline Base.ndims(int::IntegratorVPRKpTableau{DT,TT,D,S}) where {DT,TT,D,S} = D


function Integrators.initialize!(int::IntegratorVPRKpTableau, sol::AtomicSolutionPODE)
    sol.t̄ = sol.t - timestep(int)

    equation(int, :v̄)(sol.t, sol.q, sol.v)
    equation(int, :f̄)(sol.t, sol.q, sol.v, sol.f)

    initialize!(int.iguess, sol.t, sol.q, sol.p, sol.v, sol.f,
                            sol.t̄, sol.q̄, sol.p̄, sol.v̄, sol.f̄)
end


function initial_guess!(int::IntegratorVPRKpTableau{DT}, sol::AtomicSolutionPODE{DT},
                        cache::IntegratorCacheVPRK{DT}=int.caches[DT]) where {DT}
    for i in eachstage(int)
        evaluate!(int.iguess, sol.q̄, sol.p̄, sol.v̄, sol.f̄,
                              sol.q, sol.p, sol.v, sol.f,
                              cache.q̃, cache.ṽ,
                              tableau(int).c[i])

        for k in eachdim(int)
            int.solver.x[ndims(int)*(i-1)+k] = cache.ṽ[k]
        end
    end
end


function compute_stages!(x::Vector{ST}, cache::IntegratorCacheVPRK{ST}, params::ParametersVPRKpTableau{DT,TT,D,S}) where {ST,DT,TT,D,S}
    compute_stages!(x, cache.Q, cache.V,
                       cache.P, cache.F,
                       cache.Y, cache.Z,
                       cache.q̃, cache.ṽ, cache.p̃,
                       cache.θ̃, cache.y, cache.z,
                       params)
end

function compute_stages!(x::Vector{ST}, Q::Vector{Vector{ST}}, V::Vector{Vector{ST}},
                                        P::Vector{Vector{ST}}, F::Vector{Vector{ST}},
                                        Y::Vector{Vector{ST}}, Z::Vector{Vector{ST}},
                                        q::Vector{ST}, v::Vector{ST}, p::Vector{ST},
                                        θ::Vector{ST}, y::Vector{ST}, z::Vector{ST},
                                        params::ParametersVPRKpTableau{DT,TT,D,S}) where {ST,DT,TT,D,S}

    local tᵢ::TT
    local y₁::ST
    local y₂::ST
    local z₁::ST
    local z₂::ST

    # copy x to Y
    for i in 1:S
        for k in 1:D
            V[i][k] = x[D*(i-1)+k]
        end
    end

    # compute Y
    for i in 1:S
        for k in 1:D
            y₁ = 0
            y₂ = 0
            for j in 1:S
                y₁ += params.tab.a[i,j] * V[j][k]
                y₂ += params.A[i,j] * V[j][k]
            end
            Y[i][k] = y₁ + y₂
        end
    end

    # compute Q=q̄+Δt*Y
    for i in 1:S
        for k in 1:D
            Q[i][k] = params.q̄[k] + params.Δt * Y[i][k]
        end
    end

    # compute P=ϑ(Q,V) and F=f(Q,V)
    for i in 1:S
        tᵢ = params.t̄ + params.Δt * params.tab.c[i]
        params.equs[:ϑ](tᵢ, Q[i], V[i], P[i])
        params.equs[:f](tᵢ, Q[i], V[i], F[i])
    end

    # compute Z
    for i in 1:S
        for k in 1:D
            z₁ = 0
            z₂ = 0
            for j in 1:S
                z₁ += params.tab.a[i,j] * F[j][k]
                z₂ += params.A[i,j] * F[j][k]
            end
            Z[i][k] = z₁ + z₂
        end
    end

    # compute y=B*V and z=B*F
    y .= 0
    z .= 0
    for k in 1:D
        for j in 1:S
            y[k] += params.tab.b[j] * V[j][k]
            z[k] += params.tab.b[j] * F[j][k]
        end
    end


    # compute q=q̄+Δt*y and p=p̄+Δt*z
    q .= params.q̄ .+ params.Δt .* y
    p .= params.p̄ .+ params.Δt .* z

    # compute θ=ϑ(t,q)
    params.equs[:ϑ](params.t, q, v, θ)
end

"Compute stages of projected Gauss-Legendre Runge-Kutta methods."
function Integrators.function_stages!(x::Vector{ST}, b::Vector{ST},
                params::ParametersVPRKpTableau{DT,TT,D,S},
                caches::CacheDict) where {ST,DT,TT,D,S}

    # get cache for internal stages
    cache = caches[ST]

    compute_stages!(x, cache, params)

    # compute b = [P-p-AF]
    for i in 1:S
        for k in 1:D
            b[D*(i-1)+k] = cache.P[i][k] - params.p̄[k] - params.Δt * cache.Z[i][k]
        end
    end
end


function function_dirac_constraint!(λ::Vector, int::IntegratorVPRKpTableau{DT,TT},
                                    cache::IntegratorCacheVPRK{DT}=int.caches[DT]) where {DT,TT}
    # copy λ to integrator parameters
    update_tableau!(int.params, λ)

    # println(λ)

    # call nonlinear solver
    solve!(int.solver)

    # print solver status
    print_solver_status(int.solver.status, int.solver.params)

    # check if solution contains NaNs or error bounds are violated
    check_solver_status(int.solver.status, int.solver.params)

    # compute vector fields at internal stages
    compute_stages!(int.solver.x, cache, int.params)

    # compute and return ϑ(t,q)-p
    return cache.θ̃ .- cache.p̃
end


function Integrators.integrate_step!(int::IntegratorVPRKpTableau{DT,TT}, sol::AtomicSolutionPODE{DT,TT},
                                     cache::IntegratorCacheVPRK{DT}=int.caches[DT]) where {DT,TT}

    # update nonlinear solver parameters from cache
    update_params!(int.params, sol)

    # compute initial guess
    initial_guess!(int, sol, cache)

    # reset solution
    reset!(sol, timestep(int))

    # determine parameter λ
    nlres = nlsolve(λ -> function_dirac_constraint!(λ, int, cache), zero(int.params.λ);
                xtol=SimpleSolvers.get_config(:nls_atol),
                ftol=maximum(int.params.p̄)*SimpleSolvers.get_config(:nls_atol),
                iterations=100)
                #xtol=timestep(int)^nstages(int)*SimpleSolvers.get_config(:nls_atol),

    int.params.λ .= nlres.zero

    # println("x converged = ", nlres.x_converged, ", ",
    #         "f converged = ", nlres.f_converged, ", ",
    #         "residual norm = ", nlres.residual_norm, ", ",
    #         "iterations = ", nlres.iterations)
    #
    # println("λ = ", int.params.λ)
    # println("Δ = ", cache.θ̃ .- cache.p̃)
    # println("p = ", cache.p̃)
    # println()

    # compute final update
    update_solution!(sol.q, cache.V, tableau(int).b, timestep(int))
    update_solution!(sol.p, cache.F, tableau(int).b, timestep(int))

    # copy solution to initial guess
    update_vector_fields!(int.iguess, sol.t, sol.q, sol.p, sol.v, sol.f)
end
