@doc raw"""
`ParametersDGVIPI`: Parameters for right-hand side function of Discontinuous
Galerkin Variational Integrator with Path Integral approximation of the jump.

### Parameters

* `DT`: data type
* `TT`: parameter type
* `D`:  dimension of the system
* `S`:  number of basis nodes
* `QR`: number of quadrature nodes
* `FR`: number of quadrature nodes for the discontinuity

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
* `β`:  weights of the quadrature rule for the discontinuity
* `γ`:  nodes of the quadrature rule for the discontinuity
* `μ⁻`: mass vector for the lhs jump value
* `μ⁺`: mass vector for the rhs jump value
* `α⁻`: derivative vector for the discontinuity lhs value
* `α⁺`: derivative vector for the discontinuity rhs value
* `ρ⁻`: reconstruction coefficients for central jump value, lhs value
* `ρ⁺`: reconstruction coefficients for central jump value, rhs value
* `t`:  current time
* `q`:  current solution of qₙ
* `q⁻`: current solution of qₙ⁻
* `q⁺`: current solution of qₙ⁺
"""
mutable struct ParametersDGVIPI{DT, TT, D, S, QR, FR, ET <: NamedTuple} <: Parameters{DT,TT}
    equs::ET
    Δt::TT

    b::Vector{TT}
    c::Vector{TT}
    m::Matrix{TT}
    a::Matrix{TT}
    r⁻::Vector{TT}
    r⁺::Vector{TT}

    β::Vector{TT}
    γ::Vector{TT}
    μ⁻::Vector{TT}
    μ⁺::Vector{TT}
    α⁻::Vector{TT}
    α⁺::Vector{TT}
    ρ⁻::TT
    ρ⁺::TT

    t::TT

    q::Vector{DT}
    q⁻::Vector{DT}
    q⁺::Vector{DT}

    function ParametersDGVIPI{DT,D}(equs::ET, Δt::TT, b, c, m, a, r⁻, r⁺, β, γ, μ⁻, μ⁺, α⁻, α⁺, ρ⁻, ρ⁺) where {DT, TT, D, ET <: NamedTuple}
        @assert length(b)  == length(c)
        @assert length(β)  == length(γ)
        @assert length(r⁻) == length(r⁺)
        new{DT,TT,D,length(r⁻),length(c),length(γ),ET}(equs, Δt, b, c, m, a, r⁻, r⁺, β, γ, μ⁻, μ⁺, α⁻, α⁺, ρ⁻, ρ⁺,
                    zero(TT), zeros(DT,D), zeros(DT,D), zeros(DT,D))
    end

    function ParametersDGVIPI{DT,D}(equs::NamedTuple, Δt::TT,
                    basis::Basis{TT}, quadrature::QuadratureRule{TT}, jump::Discontinuity{TT}) where {DT,TT,D}

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

        β  = weights(jump.quadrature)
        γ  = nodes(jump.quadrature)
        μ⁰ = zeros(TT, 2)
        μ⁻ = zeros(TT, nnodes(jump.quadrature))
        μ⁺ = zeros(TT, nnodes(jump.quadrature))
        α⁻ = zeros(TT, nnodes(jump.quadrature))
        α⁺ = zeros(TT, nnodes(jump.quadrature))

        ρ⁻ = evaluate_l(jump.path, TT(0.5))
        ρ⁺ = evaluate_r(jump.path, TT(0.5))

        for i in 1:nnodes(jump.quadrature)
            μ⁻[i] = evaluate_l(jump.path, nodes(jump.quadrature)[i])
            μ⁺[i] = evaluate_r(jump.path, nodes(jump.quadrature)[i])
            α⁻[i] = derivative_l(jump.path, nodes(jump.quadrature)[i])
            α⁺[i] = derivative_r(jump.path, nodes(jump.quadrature)[i])
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
            println("    β = ", β)
            println("    γ = ", γ)
            println("    μ⁻= ", μ⁻)
            println("    μ⁺= ", μ⁺)
            println("    α⁻= ", α⁻)
            println("    α⁺= ", α⁺)
            println("    ρ⁻= ", ρ⁻)
            println("    ρ⁺= ", ρ⁺)
            println()
        end

        ParametersDGVIPI{DT,D}(equs, Δt, b, c, m, a, r⁻, r⁺,
                    β, γ, μ⁻, μ⁺, α⁻, α⁺, ρ⁻, ρ⁺)
    end
end

# function update_params!(params::ParametersDGVIPI, sol::AtomicSolutionPODE)
#     # set time for nonlinear solver and copy previous solution
#     params.t  = sol.t
#     params.q .= sol.q
#     # params.p .= sol.p
# end


@doc raw"""
Nonlinear function cache for Discontinuous Galerkin Variational Integrator.

### Parameters

* `X`: degrees of freedom
* `Q`: solution at quadrature nodes
* `V`: velocity at quadrature nodes
* `P`: one-form at quadrature nodes
* `F`: forces at quadrature nodes
* `q`:  current solution of qₙ
* `q⁻`: current solution of qₙ⁻
* `q⁺`: current solution of qₙ⁺
* `q̄`:  current solution of qₙ₊₁
* `q̄⁻`: current solution of qₙ₊₁⁻
* `q̄⁺`: current solution of qₙ₊₁⁺
* `λ`:  jump of the solution at tₙ
* `λ̄`:  jump of the solution at tₙ₊₁
* `ϕ`:  solution evaluated across the jump at tₙ
* `ϕ̅`:  solution evaluated across the jump at tₙ₊₁
* `θ`:  one-form evaluated across the jump at tₙ
* `Θ̅`:  one-form evaluated across the jump at tₙ₊₁
* `g`:  projection evaluated across the jump at tₙ
* `ḡ`:  projection evaluated across the jump at tₙ₊₁
"""
struct IntegratorCacheDGVIPI{ST,D,S,QR,FR} <: IODEIntegratorCache{ST,D}
    X::Vector{Vector{ST}}
    Q::Vector{Vector{ST}}
    V::Vector{Vector{ST}}
    P::Vector{Vector{ST}}
    F::Vector{Vector{ST}}

    q̃::Vector{ST}
    p̃::Vector{ST}
    ṽ::Vector{ST}
    f̃::Vector{ST}
    s̃::Vector{ST}

    q::Vector{ST}
    q⁻::Vector{ST}
    q⁺::Vector{ST}
    q̄::Vector{ST}
    q̄⁻::Vector{ST}
    q̄⁺::Vector{ST}

    λ::Vector{Vector{ST}}
    λ̄::Vector{Vector{ST}}
    ϕ::Vector{Vector{ST}}
    ϕ̅::Vector{Vector{ST}}
    θ::Vector{Vector{ST}}
    Θ̅::Vector{Vector{ST}}
    g::Vector{Vector{ST}}
    ḡ::Vector{Vector{ST}}

    function IntegratorCacheDGVIPI{ST,D,S,QR,FR}() where {ST,D,S,QR,FR}
        # create internal stage vectors
        X = create_internal_stage_vector(ST,D,S)
        Q = create_internal_stage_vector(ST,D,QR)
        V = create_internal_stage_vector(ST,D,QR)
        P = create_internal_stage_vector(ST,D,QR)
        F = create_internal_stage_vector(ST,D,QR)

        # create temporary vectors
        q̃ = zeros(ST,D)
        p̃ = zeros(ST,D)
        ṽ = zeros(ST,D)
        f̃ = zeros(ST,D)
        s̃ = zeros(ST,D)

        # create solution vectors
        q  = zeros(ST,D)
        q⁻ = zeros(ST,D)
        q⁺ = zeros(ST,D)
        q̄  = zeros(ST,D)
        q̄⁻ = zeros(ST,D)
        q̄⁺ = zeros(ST,D)

        # create jump vectors
        λ  = create_internal_stage_vector(ST,D,FR)
        λ̄ = create_internal_stage_vector(ST,D,FR)
        ϕ  = create_internal_stage_vector(ST,D,FR)
        ϕ̅  = create_internal_stage_vector(ST,D,FR)
        θ  = create_internal_stage_vector(ST,D,FR)
        Θ̅  = create_internal_stage_vector(ST,D,FR)
        g  = create_internal_stage_vector(ST,D,FR)
        ḡ  = create_internal_stage_vector(ST,D,FR)

        new(X, Q, V, P, F, q̃, p̃, ṽ, f̃, s̃, q, q⁻, q⁺, q̄, q̄⁻, q̄⁺, λ, λ̄, ϕ, ϕ̅, θ, Θ̅, g, ḡ)
    end
end

function IntegratorCache{ST}(params::ParametersDGVIPI{DT,TT,D,S,QR,FR}; kwargs...) where {ST,DT,TT,D,S,QR,FR}
    IntegratorCacheDGVIPI{ST,D,S,QR,FR}(; kwargs...)
end

@inline CacheType(ST, params::ParametersDGVIPI{DT,TT,D,S,QR,FR}) where {DT,TT,D,S,QR,FR} = IntegratorCacheDGVIPI{ST,D,S,QR,FR}


@doc raw"""
`IntegratorDGVIPI`: Discontinuous Galerkin Variational Integrator. *EXPERIMENTAL*

### Parameters

### Fields

* `equation`: Implicit Ordinary Differential Equation
* `basis`: piecewise polynomial basis
* `quadrature`: numerical quadrature rule
* `jump`: jump across discontinuity
* `Δt`: time step
* `params`: ParametersDGVIPI
* `solver`: nonlinear solver
* `iguess`: initial guess
* `q`: current solution vector for trajectory
* `q⁻`: current solution vector for trajectory, lhs of jump
* `q⁺`: current solution vector for trajectory, rhs of jump
* `cache`: temporary variables for nonlinear solver
"""
struct IntegratorDGVIPI{DT, TT, D, S, R,
                BT <: Basis,
                JT <: Discontinuity,
                PT <: ParametersDGVIPI{DT,TT,D,S},
                ST <: NonlinearSolver{DT},
                IT <: InitialGuessODE{TT}} <: IODEIntegrator{DT,TT}
    basis::BT
    quadrature::QuadratureRule{TT,R}
    jump::JT

    params::PT
    solver::ST
    iguess::IT
    caches::CacheDict{PT}

    q::Vector{DT}
    q⁻::Vector{DT}
    q⁺::Vector{DT}

    function IntegratorDGVIPI(basis::BT, quadrature::QuadratureRule{TT,R}, jump::JT,
                    params::ParametersDGVIPI{DT,TT,D,S}, solver::ST, iguess::IT, caches) where {DT,TT,D,S,R,BT,JT,ST,IT}
        # create solution vectors
        q  = zeros(DT,D)
        q⁻ = zeros(DT,D)
        q⁺ = zeros(DT,D)

        new{DT, TT, D, S, R, BT, JT, typeof(params), ST, IT}(basis, quadrature, jump, params, solver, iguess, caches, q, q⁻, q⁺)
    end

    function IntegratorDGVIPI{DT,D}(equations::NamedTuple, basis::Basis{TT}, quadrature::QuadratureRule{TT,R},
                                    jump::Discontinuity{TT}, Δt::TT;
                                    interpolation=HermiteExtrapolation{DT}) where {DT,TT,D,R}

        # get number of stages
        S = nbasis(basis)

        # create params
        params = ParametersDGVIPI{DT,D}(equations, Δt, basis, quadrature, jump)

        # create cache dict
        caches = CacheDict(params)

        # create nonlinear solver
        solver = create_nonlinear_solver(DT, D*(S+1), params, caches)

        # create initial guess
        iguess = InitialGuessODE(get_config(:ig_extrapolation), equations[:v̄], Δt)

        # create integrator
        IntegratorDGVIPI(basis, quadrature, jump, params, solver, iguess, caches)
    end

    function IntegratorDGVIPI(problem::Union{IODEProblem{DT}, LODEProblem{DT}}, basis::Basis, quadrature::QuadratureRule, jump::Discontinuity; kwargs...) where {DT}
        IntegratorDGVIPI{DT, ndims(problem)}(functions(problem), basis, quadrature, jump, timestep(problem); kwargs...)
    end
end

@inline GeometricBase.equation(integrator::IntegratorDGVIPI, i::Symbol) = integrator.params.equs[i]
@inline GeometricBase.equations(integrator::IntegratorDGVIPI) = integrator.params.equs
@inline GeometricBase.timestep(integrator::IntegratorDGVIPI) = integrator.params.Δt
@inline Base.ndims(::IntegratorDGVIPI{DT,TT,D}) where {DT,TT,D} = D


function update_params!(params::ParametersDGVIPI, int::IntegratorDGVIPI)
    # set time for nonlinear solver and copy previous solution
    params.t  += int.params.Δt
    params.q  .= int.q
    params.q⁻ .= int.q⁻
    params.q⁺ .= int.q⁺
end


function initialize!(int::IntegratorDGVIPI, sol::AtomicSolutionPODE)
    # copy initial conditions from solution
    int.q  .= sol.q
    int.q⁻ .= int.q
    int.q⁺ .= q⁺ = (int.q .- int.params.ρ⁻ .* int.q⁻) ./ int.params.ρ⁺

    sol.t̄ = sol.t - timestep(int)

    equation(int, :v̄)(sol.t, sol.q, sol.v)

    # initialise initial guess
    initialize!(int.iguess, sol.t, sol.q, sol.v,
                            sol.t̄, sol.q̄, sol.v̄)
end


function initial_guess!(int::IntegratorDGVIPI{DT,TT,D,S,R}, sol::AtomicSolutionPODE{DT,TT},
                        cache::IntegratorCacheDGVIPI{DT}=int.caches[DT]) where {DT,TT,D,S,R}
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
function function_stages!(x::Vector{ST}, b::Vector{ST}, params::ParametersDGVIPI{DT,TT,D,S,QR,FR},
                caches::CacheDict) where {ST,DT,TT,D,S,QR,FR}
    @assert length(x) == length(b)

    # get cache for internal stages
    cache = caches[ST]

    # copy old solution to cache
    cache.q  .= params.q
    cache.q⁻ .= params.q⁻

    # compute stages from nonlinear solver solution x
    compute_stages!(x, cache, params)

    # compute rhs b of nonlinear solver
    compute_rhs!(b, cache, params)
end


function compute_stages!(x, cache::IntegratorCacheDGVIPI{ST,D,S}, params::ParametersDGVIPI{DT,TT,D,S}) where {ST,DT,TT,D,S}
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

    # compute Q, qₙ⁺, qₙ₊₁⁻ and qₙ₊₁⁺
    compute_stages_q!(cache, params)

    # compute V
    compute_stages_v!(cache, params)

    # compute P and F
    compute_stages_p!(cache, params)

    # compute jump
    compute_stages_λ!(cache, params)
end


"Compute solution at quadrature nodes and across jump."
function compute_stages_q!(cache::IntegratorCacheDGVIPI{ST,D,S,R},
                           params::ParametersDGVIPI{DT,TT,D,S,R}) where {ST,DT,TT,D,S,R}

    local q::ST
    local q⁺::ST
    local q̄⁻::ST

    local X = cache.X
    local Q = cache.Q

    # compute Q
    for i in 1:R
        for k in 1:D
            q = 0
            for j in 1:S
                q += params.m[i,j] * X[j][k]
            end
            Q[i][k] = q
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

    # compute qₙ₊₁⁺
    cache.q̄⁺ .= (cache.q̄ .- params.ρ⁻ .* cache.q̄⁻) ./ params.ρ⁺
end


"Compute velocities at quadrature nodes."
function compute_stages_v!(cache::IntegratorCacheDGVIPI{ST,D,S,R},
                           params::ParametersDGVIPI{DT,TT,D,S,R}) where {ST,DT,TT,D,S,R}
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
function compute_stages_p!(cache::IntegratorCacheDGVIPI{ST,D,S,R},
                           params::ParametersDGVIPI{DT,TT,D,S,R}) where {ST,DT,TT,D,S,R}

    local tᵢ::TT

    # compute P=ϑ(Q) and F=f(Q)
    for i in 1:R
        tᵢ = params.t + params.Δt * params.c[i]
        params.equs[:ϑ](tᵢ, cache.Q[i], cache.V[i], cache.P[i])
        params.equs[:f](tᵢ, cache.Q[i], cache.V[i], cache.F[i])
    end
end


function compute_stages_λ!(cache::IntegratorCacheDGVIPI{ST,D,S,QR,FR},
                           params::ParametersDGVIPI{DT,TT,D,S,QR,FR}) where {ST,DT,TT,D,S,QR,FR}

    local t₀::TT = params.t
    local t₁::TT = params.t + params.Δt

    # compute λ and λ̄
    for i in 1:FR
        cache.λ[i] .= params.α⁻[i] .* cache.q⁻ .+ params.α⁺[i] .* cache.q⁺
        cache.λ̄[i] .= params.α⁻[i] .* cache.q̄⁻ .+ params.α⁺[i] .* cache.q̄⁺
        cache.ϕ[i] .= params.μ⁻[i] .* cache.q⁻ .+ params.μ⁺[i] .* cache.q⁺
        cache.ϕ̅[i] .= params.μ⁻[i] .* cache.q̄⁻ .+ params.μ⁺[i] .* cache.q̄⁺
    end

    for i in 1:FR
        params.equs[:ϑ](t₀, cache.ϕ[i], cache.λ[i], cache.θ[i])
        params.equs[:ϑ](t₁, cache.ϕ̅[i], cache.λ̄[i], cache.Θ̅[i])
        params.equs[:g](t₀, cache.ϕ[i], cache.λ[i], cache.g[i])
        params.equs[:g](t₁, cache.ϕ̅[i], cache.λ̄[i], cache.ḡ[i])
    end
end


function compute_rhs!(b::Vector{ST}, cache::IntegratorCacheDGVIPI{ST,D,S,QR,FR},
                params::ParametersDGVIPI{DT,TT,D,S,QR,FR}) where {ST,DT,TT,D,S,QR,FR}

    local z::ST

    # compute b = - [(P-AF)]
    for i in 1:S
        for k in 1:D
            z = 0
            for j in 1:QR
                z += params.b[j] * params.m[j,i] * cache.F[j][k] * params.Δt
                z += params.b[j] * params.a[j,i] * cache.P[j][k]
            end
            for j in 1:FR
                z += params.β[j] * params.r⁺[i] * params.α⁺[j] * cache.θ[j][k]
                z += params.β[j] * params.r⁻[i] * params.α⁻[j] * cache.Θ̅[j][k]
                z += params.β[j] * params.r⁺[i] * params.μ⁺[j] * cache.g[j][k]
                z += params.β[j] * params.r⁻[i] * params.μ⁻[j] * cache.ḡ[j][k]
            end

            b[D*(i-1)+k] = z
        end
    end

    # compute b = [qₙ - ϕ(0.5; qₙ⁻, qₙ⁺)]
    for k in 1:D
        b[D*S+k] = cache.q[k] - params.ρ⁻ * cache.q⁻[k] - params.ρ⁺ * cache.q⁺[k]
    end
end


function update_solution!(int::IntegratorDGVIPI{DT}, cache::IntegratorCacheDGVIPI{DT}) where {DT}
    int.q  .= cache.q̄
    int.q⁻ .= cache.q̄⁻
    int.q⁺ .= cache.q̄⁺
end


function integrate_step!(int::IntegratorDGVIPI{DT,TT}, sol::AtomicSolutionPODE{DT,TT},
                         cache::IntegratorCacheDGVIPI{DT}=int.caches[DT]) where {DT,TT}
    # update nonlinear solver parameters from cache
    update_params!(int.params, int)

    # compute initial guess
    initial_guess!(int, sol, cache)

    # reset cache
    reset!(sol)

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

    # copy solution to initial guess
    update_vector_fields!(int.iguess, sol.t, sol.q, sol.v)
    # update_vector_fields!(int.iguess, int.q, int.q⁻)

    # take care of periodic solutions
    # cut_periodic_solution!(int.q,  int.equation.periodicity)
    # cut_periodic_solution!(int.q⁻, int.equation.periodicity)

    # copy to solution
    # copy_solution!(sol, int.q, int.q⁻)
end
