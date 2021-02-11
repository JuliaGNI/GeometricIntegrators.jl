@doc raw"""
`ParametersDGVIP0`: Parameters for right-hand side function of Discontinuous Galerkin Variational Integrator.

### Parameters

* `DT`: data type
* `TT`: parameter type
* `D`: dimension of the system
* `S`: number of basis nodes
* `R`: number of quadrature nodes

### Fields

* `Θ`:  function of the noncanonical one-form (∂L/∂v)
* `f`:  function of the force (∂L/∂q)
* `g`:  function of the projection ∇ϑ(q)⋅v
* `Δt`: time step
* `b`:  quadrature weights
* `c`:  quadrature nodes
* `m`:  mass matrix
* `a`:  derivative matrix
* `r⁻`: reconstruction coefficients, jump lhs value
* `r⁺`: reconstruction coefficients, jump rhs value
* `t`:  current time
* `q`:  current solution of qₙ
* `q⁻`: current solution of qₙ⁻
* `θ`:  one-form ϑ evaluated on qₙ
* `θ⁻`: one-form ϑ evaluated on qₙ⁻
"""
mutable struct ParametersDGVIP0{DT, TT, D, S, R, ET <: NamedTuple} <: Parameters{DT,TT}
    equs::ET
    Δt::TT

    b::Vector{TT}
    c::Vector{TT}
    m::Matrix{TT}
    a::Matrix{TT}
    r⁻::Vector{TT}
    r⁺::Vector{TT}

    t::TT

    q::Vector{DT}
    q⁻::Vector{DT}
    θ::Vector{DT}
    θ⁻::Vector{DT}

    function ParametersDGVIP0{DT,D}(equs::ET, Δt::TT, b, c, m, a, r⁻, r⁺) where {DT, TT, D, ET <: NamedTuple}
        @assert length(b)  == length(c)
        @assert length(r⁻) == length(r⁺)
        new{DT,TT,D,length(r⁻),length(c),ET}(equs, Δt, b, c, m, a, r⁻, r⁺, zero(TT), zeros(DT,D), zeros(DT,D), zeros(DT,D), zeros(DT,D))
    end

    function ParametersDGVIP0{DT,D}(equs::NamedTuple, Δt::TT,
                    basis::Basis{TT}, quadrature::QuadratureRule{TT}) where {DT,TT,D}

        # compute coefficients
        b = weights(quadrature)
        c = nodes(quadrature)
        m = zeros(TT, nnodes(quadrature), nbasis(basis))
        a = zeros(TT, nnodes(quadrature), nbasis(basis))
        r⁻= zeros(TT, nbasis(basis))
        r⁺= zeros(TT, nbasis(basis))

        for i in eachindex(basis)
            for j in eachindex(quadrature)
                m[j,i] = basis[nodes(quadrature)[j], i]
                a[j,i] = basis'[nodes(quadrature)[j], i]
            end
            r⁻[i] = basis[one(TT), i]
            r⁺[i] = basis[zero(TT), i]
        end

        if get_config(:verbosity) > 1
            println()
            println("  Discontinuous Galerkin Variational Integrator")
            println("  =============================================")
            println()
            println("    b = ", b)
            println("    c = ", c)
            println("    m = ", m)
            println("    a = ", a)
            println("    r⁻= ", r⁻)
            println("    r⁺= ", r⁺)
            println()
        end

        ParametersDGVIP0{DT,D}(equs, Δt, b, c, m, a, r⁻, r⁺)
    end
end


# function update_params!(params::ParametersDGVIP0, sol::AtomicSolutionPODE)
#     # set time for nonlinear solver and copy previous solution
#     params.t  = sol.t
#     params.q .= sol.q
#     params.equs[:ϑ] .= sol.p
# end


function IntegratorCache{ST}(params::ParametersDGVIP0{DT,TT,D,S,R}; kwargs...) where {ST,DT,TT,D,S,R}
    IntegratorCacheDGVI{ST,D,S,R}(; kwargs...)
end

@inline CacheType(ST, params::ParametersDGVIP0{DT,TT,D,S,R}) where {DT,TT,D,S,R} = IntegratorCacheDGVI{ST,D,S,R}


@doc raw"""
`IntegratorDGVIP0`: Discontinuous Galerkin Variational Integrator. *EXPERIMENTAL*

### Parameters

### Fields

* `equation`: Implicit Ordinary Differential Equation
* `basis`: piecewise polynomial basis
* `quadrature`: numerical quadrature rule
* `Δt`: time step
* `params`: ParametersDGVIP0
* `solver`: nonlinear solver
* `iguess`: initial guess
* `q`: current solution vector for trajectory
* `p`: current solution vector for one-form
* `cache`: temporary variables for nonlinear solver
"""
struct IntegratorDGVIP0{DT, TT, D, S, R,
                BT <: Basis,
                PT <: ParametersDGVIP0{DT,TT,D,S},
                ST <: NonlinearSolver{DT},
                IT <: InitialGuessODE{TT}} <: IODEIntegrator{DT,TT}
    basis::BT
    quadrature::QuadratureRule{TT,R}

    params::PT
    solver::ST
    iguess::IT
    caches::CacheDict{PT}

    q::Vector{DT}
    q⁻::Vector{DT}
    q⁺::Vector{DT}

    θ::Vector{DT}
    θ⁻::Vector{DT}
    θ⁺::Vector{DT}

    function IntegratorDGVIP0(basis::BT, quadrature::QuadratureRule{TT,R}, params::ParametersDGVIP0{DT,TT,D,S},
                    solver::ST, iguess::IT, caches) where {DT,TT,D,S,R,BT,ST,IT}
        # create solution vectors
        q  = zeros(DT,D)
        q⁻ = zeros(DT,D)
        q⁺ = zeros(DT,D)

        θ  = zeros(DT,D)
        θ⁻ = zeros(DT,D)
        θ⁺ = zeros(DT,D)

        new{DT, TT, D, S, R, BT, typeof(params), ST, IT}(basis, quadrature, params, solver, iguess, caches, q, q⁻, q⁺, θ, θ⁻, θ⁺)
    end

    function IntegratorDGVIP0{DT,D}(equations::NamedTuple, basis::Basis{TT}, quadrature::QuadratureRule{TT,R}, Δt::TT;
                                  interpolation=HermiteInterpolation{DT}) where {DT,TT,D,R}

        # get number of stages
        S = nbasis(basis)

        # create params
        params = ParametersDGVIP0{DT,D}(equations, Δt, basis, quadrature)

        # create cache dict
        caches = CacheDict(params)

        # create nonlinear solver
        solver = create_nonlinear_solver(DT, D*(S+1), params, caches)

        # create initial guess
        iguess = InitialGuessODE(get_config(:ig_interpolation), equations[:v̄], Δt)

        # create integrator
        IntegratorDGVIP0(basis, quadrature, params, solver, iguess, caches)
    end

    function IntegratorDGVIP0(equation::Union{IODE{DT,TT}, LODE{DT,TT}}, basis::Basis{TT}, quadrature::QuadratureRule{TT}, Δt::TT; kwargs...) where {DT,TT}
        IntegratorDGVIP0{DT, ndims(equation)}(get_function_tuple(equation), basis, quadrature, Δt; kwargs...)
    end
end

@inline equation(integrator::IntegratorDGVIP0, i::Symbol) = integrator.params.equs[i]
@inline equations(integrator::IntegratorDGVIP0) = integrator.params.equs
@inline timestep(integrator::IntegratorDGVIP0) = integrator.params.Δt
@inline Base.ndims(::IntegratorDGVIP0{DT,TT,D}) where {DT,TT,D} = D


function update_params!(params::ParametersDGVIP0, int::IntegratorDGVIP0)
    # set time for nonlinear solver and copy previous solution
    params.t  += int.params.Δt
    params.q  .= int.q
    params.q⁻ .= int.q⁻
    params.θ  .= int.θ
    params.θ⁻ .= int.θ⁻
end


function initialize!(int::IntegratorDGVIP0, sol::AtomicSolutionPODE)
    # copy initial conditions from solution
    int.q  .= sol.q
    int.θ  .= sol.p
    int.q⁻ .= int.q
    int.θ⁻ .= int.θ

    sol.t̄ = sol.t - timestep(int)

    equation(int, :v̄)(sol.t, sol.q, sol.v)

    # initialise initial guess
    initialize!(int.iguess, sol.t, sol.q, sol.v,
                            sol.t̄, sol.q̄, sol.v̄)
end


function initial_guess!(int::IntegratorDGVIP0{DT,TT, D, S, R}, sol::AtomicSolutionPODE{DT,TT},
                        cache::IntegratorCacheDGVI{DT}=int.caches[DT]) where {DT,TT,D,S,R}
    if nbasis(int.basis) > 0
        for i in eachindex(int.basis)
            evaluate!(int.iguess, sol.q̄, sol.v̄,
                                  sol.q, sol.v,
                                  cache.q̃,
                                  grid(int.basis)[i])

            for k in 1:D
                int.solver.x[D*(i-1)+k] = cache.q̃[k]
            end
        end
    else
        for i in 1:S
            for k in 1:D
                int.solver.x[D*(i-1)+k] = 0
            end
        end
    end

    evaluate!(int.iguess, sol.q̄, sol.v̄,
                          sol.q, sol.v,
                          cache.q̃,
                          one(TT))

    for k in 1:D
        int.solver.x[D*S+k] = cache.q̃[k]
    end
end


"Compute stages of variational partitioned Runge-Kutta methods."
function function_stages!(x::Vector{ST}, b::Vector{ST}, params::ParametersDGVIP0{DT,TT,D,S,R},
                caches::CacheDict) where {ST,DT,TT,D,S,R}
    @assert length(x) == length(b)

    # get cache for internal stages
    cache = caches[ST]

    # compute stages from nonlinear solver solution x
    compute_stages!(x, cache, params)

    # compute rhs b of nonlinear solver
    compute_rhs!(b, cache, params)
end


function compute_stages!(x, cache::IntegratorCacheDGVI{ST,D,S}, params::ParametersDGVIP0{DT,TT,D,S}) where {ST,DT,TT,D,S}
    # copy x to X
    for i in 1:S
        for k in 1:D
            cache.X[i][k] = x[D*(i-1)+k]
        end
    end

    # copy x to q̄=qₙ+₁
    for k in 1:D
        cache.q̄[k] = x[D*S+k]
    end

    # compute Q, qₙ⁺ and qₙ₊₁⁻
    compute_stages_q!(cache, params)

    # compute V
    compute_stages_v!(cache, params)

    # compute P and F
    compute_stages_p!(cache, params)

    # compute jump
    compute_stages_λ!(cache, params)
end


"Compute solution at quadrature nodes and across jump."
function compute_stages_q!(cache::IntegratorCacheDGVI{ST,D,S,R},
                           params::ParametersDGVIP0{DT,TT,D,S,R}) where {ST,DT,TT,D,S,R}

    local q::ST
    local q⁺::ST
    local q̄⁻::ST

    local X = cache.X

    # copy q and q⁻
    cache.q  .= params.q
    cache.q⁻ .= params.q⁻

    # compute Q
    for i in 1:R
        for k in 1:D
            q = 0
            for j in 1:S
                q += params.m[i,j] * X[j][k]
            end
            cache.Q[i][k] = q
        end
    end

    # compute qₙ⁺ and qₙ₊₁⁻
    for k in 1:D
        q⁺ = 0
        q̄⁻ = 0
        for i in 1:S
            q⁺ += params.r⁺[i] * X[i][k]
            q̄⁻ += params.r⁻[i] * X[i][k]
        end
        cache.q⁺[k] = q⁺
        cache.q̄⁻[k] = q̄⁻
    end
end


"Compute velocities at quadrature nodes."
function compute_stages_v!(cache::IntegratorCacheDGVI{ST,D,S,R},
                           params::ParametersDGVIP0{DT,TT,D,S,R}) where {ST,DT,TT,D,S,R}
    local v::ST

    for i in 1:R
        for k in 1:D
            v = 0
            for j in 1:S
                v += params.a[i,j] * cache.X[j][k]
            end
            cache.V[i][k] = v / params.Δt
        end
    end
end


"Compute one-form and forces at quadrature nodes."
function compute_stages_p!(cache::IntegratorCacheDGVI{ST,D,S,R},
                           params::ParametersDGVIP0{DT,TT,D,S,R}) where {ST,DT,TT,D,S,R}

    local tᵢ::TT

    # copy q and q⁻
    cache.θ  .= params.θ
    cache.θ⁻ .= params.θ⁻

    # compute P=ϑ(Q) and F=f(Q)
    for i in 1:R
        tᵢ = params.t + params.Δt * params.c[i]
        params.equs[:ϑ](tᵢ, cache.Q[i], cache.V[i], cache.P[i])
        params.equs[:f](tᵢ, cache.Q[i], cache.V[i], cache.F[i])
    end
end


function compute_stages_λ!(cache::IntegratorCacheDGVI{ST,D,S,R},
                           params::ParametersDGVIP0{DT,TT,D,S,R}) where {ST,DT,TT,D,S,R}

    local t₀::TT = params.t
    local t₁::TT = params.t + params.Δt

    local Θ̅⁻::ST
    local Θ̅⁺::ST

    # compute ϕ and ϕ̅
    cache.ϕ  .= 0.5 * (cache.q⁻ .+ cache.q⁺)
    cache.ϕ⁻ .= 0.5 * (cache.q⁻ .+ cache.q )
    cache.ϕ⁺ .= 0.5 * (cache.q  .+ cache.q⁺)

    cache.ϕ̅  .= 0.5 * (cache.q̄⁻ .+ cache.q̄⁺)
    cache.ϕ̅⁻ .= 0.5 * (cache.q̄⁻ .+ cache.q̄ )
    cache.ϕ̅⁺ .= 0.5 * (cache.q̄  .+ cache.q̄⁺)

    # compute λ and λ̄
    cache.λ  .= cache.q⁺ .- cache.q⁻
    cache.λ⁻ .= cache.q  .- cache.q⁻
    cache.λ⁺ .= cache.q⁺ .- cache.q

    cache.λ̄ .= cache.q̄⁺ .- cache.q̄⁻
    cache.λ̄⁻ .= cache.q̄  .- cache.q̄⁻
    cache.λ̄⁺ .= cache.q̄⁺ .- cache.q̄

    # compute ϑ
    params.equs[:ϑ](t₀, cache.q,  cache.q,  cache.θ)
    # params.equs[:ϑ](t₀, cache.q⁻, cache.q⁻, cache.θ⁻)
    # params.equs[:ϑ](t₀, cache.q⁺, cache.q⁺, cache.θ⁺)

    params.equs[:ϑ](t₁, cache.q̄,  cache.q̄,  cache.Θ̅)
    # params.equs[:ϑ](t₁, cache.q̄⁻, cache.q̄⁻, cache.Θ̅⁻)
    # params.equs[:ϑ](t₁, cache.q̄⁺, cache.q̄⁺, cache.Θ̅⁺)

    for k in 1:D
        θ⁺ = 0
        Θ̅⁻ = 0
        for i in 1:S
            θ⁺ += params.r⁺[i] * cache.P[i][k]
            Θ̅⁻ += params.r⁻[i] * cache.P[i][k]
        end
        cache.θ⁺[k] = θ⁺
        cache.Θ̅⁻[k] = Θ̅⁻
    end

    # compute projection
    params.equs[:g](t₀, cache.q,  cache.λ,  cache.g)
    params.equs[:g](t₀, cache.q⁻, cache.λ⁻, cache.g⁻)
    params.equs[:g](t₀, cache.q⁺, cache.λ⁺, cache.g⁺)

    params.equs[:g](t₁, cache.q̄,  cache.λ̄,  cache.ḡ)
    params.equs[:g](t₁, cache.q̄⁻, cache.λ̄⁻, cache.ḡ⁻)
    params.equs[:g](t₁, cache.q̄⁺, cache.λ̄⁺, cache.ḡ⁺)

    # # compute ϑ
    # params.equs[:ϑ](t₀, cache.ϕ⁻, cache.ϕ⁻, cache.θ⁻)
    # params.equs[:ϑ](t₀, cache.ϕ⁺, cache.ϕ⁺, cache.θ⁺)
    # params.equs[:ϑ](t₁, cache.ϕ̅⁻, cache.ϕ̅⁻, cache.Θ̅⁻)
    #
    # # compute projection
    # params.equs[:g](t₀, cache.ϕ⁻, cache.λ⁻, cache.g⁻)
    # params.equs[:g](t₀, cache.ϕ⁺, cache.λ⁺, cache.g⁺)
    # params.equs[:g](t₁, cache.ϕ̅⁻, cache.λ̄⁻, cache.ḡ⁻)
end


function compute_rhs!(b::Vector{ST}, cache::IntegratorCacheDGVI{ST,D,S,R},
                params::ParametersDGVIP0{DT,TT,D,S,R}) where {ST,DT,TT,D,S,R}

    local z::ST

    # compute b = - [(P-AF)]
    for i in 1:S
        for k in 1:D
            z = 0
            for j in 1:R
                z += params.b[j] * params.m[j,i] * cache.F[j][k] * params.Δt
                z += params.b[j] * params.a[j,i] * cache.P[j][k]
            end

            z += params.r⁺[i] * 0.5 * ( cache.θ[k] + cache.θ⁺[k] )
            z -= params.r⁻[i] * 0.5 * ( cache.Θ̅[k] + cache.Θ̅⁻[k] )

            z += params.r⁺[i] * 0.5 * cache.g⁺[k]
            z += params.r⁻[i] * 0.5 * cache.ḡ⁻[k]

            b[D*(i-1)+k] = z
        end
    end

    for k in 1:D
        b[D*S+k] = cache.θ⁺[k] - cache.θ⁻[k] - cache.g[k]
    end
end


function update_solution!(int::IntegratorDGVIP0{DT}, cache::IntegratorCacheDGVI{DT}) where {DT}
    int.q  .= cache.q̄
    int.q⁻ .= cache.q̄⁻
    int.q⁺ .= cache.q̄⁺

    int.θ  .= cache.Θ̅
    int.θ⁻ .= cache.Θ̅⁻
    int.θ⁺ .= cache.Θ̅⁺
end


function integrate_step!(int::IntegratorDGVIP0{DT,TT}, sol::AtomicSolutionPODE{DT,TT},
                         cache::IntegratorCacheDGVI{DT}=int.caches[DT]) where {DT,TT}
    # update nonlinear solver parameters from cache
    update_params!(int.params, int)

    # compute initial guess
    initial_guess!(int, sol, cache)

    # reset cache
    reset!(sol, timestep(int))

    # call nonlinear solver
    solve!(int.solver)

    # print solver status
    print_solver_status(int.solver.status, int.solver.params)

    # check if solution contains NaNs or error bounds are violated
    check_solver_status(int.solver.status, int.solver.params)

    # compute vector fields at internal stages
    compute_stages!(int.solver.x, cache, int.params)

    # copy solution from cache to integrator
    update_solution!(int, cache)
    sol.q = int.q
    sol.p = int.θ

    # copy solution to initial guess
    update_vector_fields!(int.iguess, sol.t, sol.q, sol.v)
    # update_vector_fields!(int.iguess, int.q, int.θ)

    # take care of periodic solutions
    # cut_periodic_solution!(int.q,  int.equation.periodicity)
    # cut_periodic_solution!(int.q⁻, int.equation.periodicity)
    # cut_periodic_solution!(int.q⁺, int.equation.periodicity)
end
