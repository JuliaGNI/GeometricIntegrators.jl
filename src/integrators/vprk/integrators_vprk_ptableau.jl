
using NLsolve

"Parameters for right-hand side function of projected Gauss-Legendre Runge-Kutta methods."
mutable struct ParametersVPRKpTableau{DT,TT,D,S,ET} <: Parameters{DT,TT}
    equ::ET
    tab::CoefficientsPGLRK{TT}
    Δt::TT

    t̅::TT
    t::TT

    q̅::Vector{DT}
    p̅::Vector{DT}

    λ::Vector{DT}
    A::Matrix{DT}
    W::Matrix{DT}

    function ParametersVPRKpTableau{DT,TT,D,S,ET}(equ, tab, Δt) where {DT,TT,D,S,ET}
        new(equ, tab, Δt, zero(TT), zero(TT), zeros(DT,D), zeros(DT,D), zeros(DT,D), zeros(DT,S,S), zeros(DT,S,S))
    end
end


"""
Projected Variational Gauss-Legendre Runge-Kutta integrator.
"""
struct IntegratorVPRKpTableau{DT, TT, PT <: ParametersVPRKpTableau{DT,TT},
                               ST <: NonlinearSolver{DT},
                               IT <: InitialGuessPODE{DT,TT}, D, S} <: IntegratorPRK{DT,TT}
    params::PT
    solver::ST
    iguess::IT
    cache::IntegratorCacheVPRK{DT,D,S}

    function IntegratorVPRKpTableau(params::ParametersVPRKpTableau{DT,TT,D,S,ET}, solver::ST, iguess::IT) where {DT, TT, D, S, ET, ST, IT}
        # create cache
        cache = IntegratorCacheVPRK{DT,D,S}()

        # create integrator
        new{DT, TT, typeof(params), ST, IT, D, S}(params, solver, iguess, cache)
    end
end

function IntegratorVPRKpTableau(equation::IODE{DT,TT,ΑT,FT,GT,HT,VT}, tableau::CoefficientsPGLRK{TT}, Δt::TT) where {DT,TT,ΑT,FT,GT,HT,VT}
    D = equation.d
    M = equation.n
    S = tableau.s

    # create params
    params = ParametersVPRKpTableau{DT,TT,D,S,typeof(equation)}(equation, tableau, Δt)

    # create solver
    solver  = create_nonlinear_solver(DT, D*S, params)

    # create initial guess
    iguess = InitialGuessPODE(get_config(:ig_interpolation), equation, Δt)

    # create integrator
    IntegratorVPRKpTableau(params, solver, iguess)
end

@inline equation(integrator::IntegratorVPRKpTableau) = integrator.params.equ
@inline nstages(integrator::IntegratorVPRKpTableau{DT,TT,PT,ST,IT,D,S}) where {DT,TT,PT,ST,IT,D,S} = S
@inline Base.ndims(int::IntegratorVPRKpTableau{DT,TT,PT,ST,IT,D,S}) where {DT,TT,PT,ST,IT,D,S} = D


function update_params!(params::ParametersVPRKpTableau, sol::AtomicSolutionPODE)
    # set time for nonlinear solver and copy previous solution
    params.t̅  = sol.t
    params.t  = sol.t + params.Δt
    params.q̅ .= sol.q
    params.p̅ .= sol.p
end

function update_tableau!(params::ParametersVPRKpTableau{DT,TT,D,S}, λ::Vector) where {DT,TT,D,S}
    # copy λ to parameters
    params.λ .= λ

    # compute tableaus
    for k in eachindex(λ)
        params.W[S-k+1, S-k] = +λ[k]
        params.W[S-k, S-k+1] = -λ[k]
    end

    # TODO make this more allocation-efficient
    params.A .= params.tab.P * params.W * params.tab.Q
    # simd_mult!(params.T, params.W, params.tab.Q)
    # simd_mult!(params.A, params.tab.P, params.T)
end


function compute_stages!(x::Vector{ST}, cache::IntegratorCacheVPRK{ST}, params::ParametersVPRKpTableau{DT,TT,D,S}) where {ST,DT,TT,D,S}
    compute_stages!(x, cache.Q, cache.V,
                       cache.P, cache.F,
                       cache.Y, cache.Z,
                       cache.q̃, cache.p̃, cache.θ̃,
                       cache.ṽ, cache.f̃,
                       params)
end

function compute_stages!(x::Vector{ST}, Q::Vector{Vector{ST}}, V::Vector{Vector{ST}},
                                        P::Vector{Vector{ST}}, F::Vector{Vector{ST}},
                                        Y::Vector{Vector{ST}}, Z::Vector{Vector{ST}},
                                        q::Vector{ST}, p::Vector{ST}, θ::Vector{ST},
                                        y::Vector{ST}, z::Vector{ST},
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

    # compute Q=q̅+Δt*Y
    for i in 1:S
        for k in 1:D
            Q[i][k] = params.q̅[k] + params.Δt * Y[i][k]
        end
    end

    # compute P=ϑ(Q,V) and F=f(Q,V)
    for i in 1:S
        tᵢ = params.t̅ + params.Δt * params.tab.c[i]
        params.equ.ϑ(tᵢ, Q[i], V[i], P[i])
        params.equ.f(tᵢ, Q[i], V[i], F[i])
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


    # compute q=q̅+Δt*y and p=p̅+Δt*z
    q .= params.q̅ .+ params.Δt .* y
    p .= params.p̅ .+ params.Δt .* z

    # compute θ=ϑ(t,q)
    params.equ.ϑ(params.t, q, θ)
end

"Compute stages of projected Gauss-Legendre Runge-Kutta methods."
@generated function function_stages!(x::Vector{ST}, b::Vector{ST}, params::ParametersVPRKpTableau{DT,TT,D,S}) where {ST,DT,TT,D,S}

    cache = IntegratorCacheVPRK{ST,D,S}()

    quote
        compute_stages!(x, $cache, params)

        # compute b = [P-p-AF]
        for i in 1:S
            for k in 1:D
                b[D*(i-1)+k] = $cache.P[i][k] - params.p̅[k] - params.Δt * $cache.Z[i][k]
            end
        end
    end
end


function function_dirac_constraint!(λ::Vector, int::IntegratorVPRKpTableau{DT,TT}) where {DT,TT}
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
    compute_stages!(int.solver.x, int.cache, int.params)

    # compute and return ϑ(t,q)-p
    return int.cache.θ̃ .- int.cache.p̃
end


function initialize!(int::IntegratorVPRKpTableau, sol::AtomicSolutionPODE)
    sol.t̅ = sol.t - timestep(int)

    equation(int).v(sol.t, sol.q, sol.p, sol.v)
    equation(int).f(sol.t, sol.q, sol.p, sol.f)

    initialize!(int.iguess, sol.t, sol.q, sol.p, sol.v, sol.f,
                            sol.t̅, sol.q̅, sol.p̅, sol.v̅, sol.f̅)
end


function initial_guess!(int::IntegratorVPRKpTableau, sol::AtomicSolutionPODE)
    for i in eachstage(int)
        evaluate!(int.iguess, sol.q, sol.p, sol.v, sol.f,
                              sol.q̅, sol.p̅, sol.v̅, sol.f̅,
                              int.cache.q̃, int.cache.ṽ,
                              tableau(int).c[i])

        for k in eachdim(int)
            int.solver.x[ndims(int)*(i-1)+k] = int.cache.ṽ[k]
        end
    end
end


"Integrate ODE with projected Gauss-Legendre Runge-Kutta integrator."
function integrate_step!(int::IntegratorVPRKpTableau{DT,TT}, sol::AtomicSolutionPODE{DT,TT}) where {DT,TT}
    # update nonlinear solver parameters from cache
    update_params!(int.params, sol)

    # compute initial guess
    initial_guess!(int, sol)

    # reset solution
    reset!(sol, timestep(int))

    # determine parameter λ
    nlres = nlsolve(λ -> function_dirac_constraint!(λ, int), zero(int.params.λ);
                xtol=get_config(:nls_atol),
                ftol=maximum(int.params.p̅)*get_config(:nls_atol),
                iterations=100)
                #xtol=timestep(int)^nstages(int)*get_config(:nls_atol),

    int.params.λ .= nlres.zero

    # println("x converged = ", nlres.x_converged, ", ",
    #         "f converged = ", nlres.f_converged, ", ",
    #         "residual norm = ", nlres.residual_norm, ", ",
    #         "iterations = ", nlres.iterations)
    #
    # println("λ = ", int.params.λ)
    # println("Δ = ", int.cache.θ̃ .- int.cache.p̃)
    # println("p = ", int.cache.p̃)
    # println()

    # compute final update
    update_solution!(sol.q, int.cache.V, tableau(int).b, timestep(int))
    update_solution!(sol.p, int.cache.F, tableau(int).b, timestep(int))

    # copy solution to initial guess
    update!(int.iguess, sol.t, sol.q, sol.p, sol.v, sol.f)
end
