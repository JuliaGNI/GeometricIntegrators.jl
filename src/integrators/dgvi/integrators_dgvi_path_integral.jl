"""
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
mutable struct ParametersDGVIPI{DT,TT,D,S,QR,FR,ΘT,FT,GT} <: Parameters{DT,TT}
    Θ::ΘT
    f::FT
    g::GT

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
end

function ParametersDGVIPI(Θ::ΘT, f::FT, g::GT, Δt::TT,
                b::Vector{TT}, c::Vector{TT}, m::Matrix{TT}, a::Matrix{TT},
                r⁻::Vector{TT}, r⁺::Vector{TT}, β::Vector{TT}, γ::Vector{TT},
                μ⁻::Vector{TT}, μ⁺::Vector{TT}, α⁻::Vector{TT}, α⁺::Vector{TT}, ρ⁻::TT, ρ⁺::TT,
                q::Vector{DT}, q⁻::Vector{DT}, q⁺::Vector{DT}) where {DT,TT,ΘT,FT,GT}

    @assert length(q)  == length(q⁻) == length(q⁺)
    @assert length(b)  == length(c)
    @assert length(β)  == length(γ)
    @assert length(r⁻) == length(r⁺)

    D  = length(q)
    S  = length(r⁻)
    QR = length(c)
    FR = length(γ)

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

    ParametersDGVIPI{DT,TT,D,S,QR,FR,ΘT,FT,GT}(
                Θ, f, g, Δt, b, c, m, a, r⁻, r⁺, β, γ, μ⁻, μ⁺, α⁻, α⁺, ρ⁻, ρ⁺, 0, q, q⁻, q⁺)
end


function ParametersDGVIPI(Θ::ΘT, f::FT, g::GT, Δt::TT,
                basis::Basis{TT}, quadrature::Quadrature{TT}, jump::Discontinuity{TT},
                q::Vector{DT}, q⁻::Vector{DT}) where {DT,TT,ΘT,FT,GT}

    # compute coefficients
    m = zeros(TT, nnodes(quadrature), nbasis(basis))
    a = zeros(TT, nnodes(quadrature), nbasis(basis))
    r⁻= zeros(TT, nbasis(basis))
    r⁺= zeros(TT, nbasis(basis))

    for i in 1:nbasis(basis)
        for j in 1:nnodes(quadrature)
            m[j,i] = evaluate(basis, i, nodes(quadrature)[j])
            a[j,i] = derivative(basis, i, nodes(quadrature)[j])
        end
        r⁻[i] = evaluate(basis, i, one(TT))
        r⁺[i] = evaluate(basis, i, zero(TT))
    end

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

    q⁺ = (q .- ρ⁻ .* q⁻) ./ ρ⁺

    ParametersDGVIPI(Θ, f, g, Δt, weights(quadrature), nodes(quadrature), m, a, r⁻, r⁺,
                weights(jump.quadrature), nodes(jump.quadrature), μ⁻, μ⁺, α⁻, α⁺, ρ⁻, ρ⁺,
                q, q⁻, q⁺)
end


"""
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
* `q̅`:  current solution of qₙ₊₁
* `q̅⁻`: current solution of qₙ₊₁⁻
* `q̅⁺`: current solution of qₙ₊₁⁺
* `λ`:  jump of the solution at tₙ
* `λ̅`:  jump of the solution at tₙ₊₁
* `ϕ`:  solution evaluated across the jump at tₙ
* `ϕ̅`:  solution evaluated across the jump at tₙ₊₁
* `θ`:  one-form evaluated across the jump at tₙ
* `Θ̅`:  one-form evaluated across the jump at tₙ₊₁
* `g`:  projection evaluated across the jump at tₙ
* `g̅`:  projection evaluated across the jump at tₙ₊₁
"""
struct NonlinearFunctionCacheDGVIPI{ST,D,S,QR,FR}
    X::Vector{Vector{ST}}
    Q::Vector{Vector{ST}}
    V::Vector{Vector{ST}}
    P::Vector{Vector{ST}}
    F::Vector{Vector{ST}}

    q::Vector{ST}
    q⁻::Vector{ST}
    q⁺::Vector{ST}
    q̅::Vector{ST}
    q̅⁻::Vector{ST}
    q̅⁺::Vector{ST}

    λ::Vector{Vector{ST}}
    λ̅::Vector{Vector{ST}}
    ϕ::Vector{Vector{ST}}
    ϕ̅::Vector{Vector{ST}}
    θ::Vector{Vector{ST}}
    Θ̅::Vector{Vector{ST}}
    g::Vector{Vector{ST}}
    g̅::Vector{Vector{ST}}

    function NonlinearFunctionCacheDGVIPI{ST,D,S,QR,FR}() where {ST,D,S,QR,FR}
        # create internal stage vectors
        X = create_internal_stage_vector(ST,D,S)
        Q = create_internal_stage_vector(ST,D,QR)
        V = create_internal_stage_vector(ST,D,QR)
        P = create_internal_stage_vector(ST,D,QR)
        F = create_internal_stage_vector(ST,D,QR)

        # create solution vectors
        q  = zeros(ST,D)
        q⁻ = zeros(ST,D)
        q⁺ = zeros(ST,D)
        q̅  = zeros(ST,D)
        q̅⁻ = zeros(ST,D)
        q̅⁺ = zeros(ST,D)

        # create jump vectors
        λ  = create_internal_stage_vector(ST,D,FR)
        λ̅  = create_internal_stage_vector(ST,D,FR)
        ϕ  = create_internal_stage_vector(ST,D,FR)
        ϕ̅  = create_internal_stage_vector(ST,D,FR)
        θ  = create_internal_stage_vector(ST,D,FR)
        Θ̅  = create_internal_stage_vector(ST,D,FR)
        g  = create_internal_stage_vector(ST,D,FR)
        g̅  = create_internal_stage_vector(ST,D,FR)

        new(X, Q, V, P, F, q, q⁻, q⁺, q̅, q̅⁻, q̅⁺, λ, λ̅, ϕ, ϕ̅, θ, Θ̅, g, g̅)
    end
end


"Compute stages of variational partitioned Runge-Kutta methods."
@generated function function_stages!(x::Vector{ST}, b::Vector{ST}, params::ParametersDGVIPI{DT,TT,D,S,QR,FR}) where {ST,DT,TT,D,S,QR,FR}
    cache = NonlinearFunctionCacheDGVIPI{ST,D,S,QR,FR}()

    quote
        @assert length(x) == length(b)

        $cache.q  .= params.q
        $cache.q⁻ .= params.q⁻

        compute_stages!(x, $cache, params)

        compute_rhs!(b, $cache, params)
    end
end


function compute_stages!(x, cache::NonlinearFunctionCacheDGVIPI{ST,D,S}, params::ParametersDGVIPI{DT,TT,D,S}) where {ST,DT,TT,D,S}
    # copy x to X
    for i in 1:S
        for k in 1:D
            cache.X[i][k] = x[D*(i-1)+k]
        end
    end

    # copy x to q̅=qₙ+₁
    for k in 1:D
        cache.q̅[k] = x[D*S+k]
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
function compute_stages_q!(cache::NonlinearFunctionCacheDGVIPI{ST,D,S,R},
                           params::ParametersDGVIPI{DT,TT,D,S,R}) where {ST,DT,TT,D,S,R}

    local q::ST
    local q⁺::ST
    local q̅⁻::ST

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
        q̅⁻ = 0
        for i in 1:S
            q⁺ += params.r⁺[i] * X[i][k]
            q̅⁻ += params.r⁻[i] * X[i][k]
        end
        cache.q⁺[k] = q⁺
        cache.q̅⁻[k] = q̅⁻
    end

    # compute qₙ₊₁⁺
    cache.q̅⁺ .= (cache.q̅ .- params.ρ⁻ .* cache.q̅⁻) ./ params.ρ⁺
end


"Compute velocities at quadrature nodes."
function compute_stages_v!(cache::NonlinearFunctionCacheDGVIPI{ST,D,S,R},
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
function compute_stages_p!(cache::NonlinearFunctionCacheDGVIPI{ST,D,S,R},
                           params::ParametersDGVIPI{DT,TT,D,S,R}) where {ST,DT,TT,D,S,R}

    local tᵢ::TT

    # compute P=ϑ(Q) and F=f(Q)
    for i in 1:R
        tᵢ = params.t + params.Δt * params.c[i]
        params.Θ(tᵢ, cache.Q[i], cache.V[i], cache.P[i])
        params.f(tᵢ, cache.Q[i], cache.V[i], cache.F[i])
    end
end


function compute_stages_λ!(cache::NonlinearFunctionCacheDGVIPI{ST,D,S,QR,FR},
                           params::ParametersDGVIPI{DT,TT,D,S,QR,FR}) where {ST,DT,TT,D,S,QR,FR}

    local t₀::TT = params.t
    local t₁::TT = params.t + params.Δt

    # compute λ and λ̅
    for i in 1:FR
        cache.λ[i] .= params.α⁻[i] .* cache.q⁻ .+ params.α⁺[i] .* cache.q⁺
        cache.λ̅[i] .= params.α⁻[i] .* cache.q̅⁻ .+ params.α⁺[i] .* cache.q̅⁺
        cache.ϕ[i] .= params.μ⁻[i] .* cache.q⁻ .+ params.μ⁺[i] .* cache.q⁺
        cache.ϕ̅[i] .= params.μ⁻[i] .* cache.q̅⁻ .+ params.μ⁺[i] .* cache.q̅⁺
    end

    for i in 1:FR
        params.Θ(t₀, cache.ϕ[i], cache.λ[i], cache.θ[i])
        params.Θ(t₁, cache.ϕ̅[i], cache.λ̅[i], cache.Θ̅[i])
        params.g(t₀, cache.ϕ[i], cache.λ[i], cache.g[i])
        params.g(t₁, cache.ϕ̅[i], cache.λ̅[i], cache.g̅[i])
    end
end


function compute_rhs!(b::Vector{ST}, cache::NonlinearFunctionCacheDGVIPI{ST,D,S,QR,FR},
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
                z += params.β[j] * params.r⁻[i] * params.μ⁻[j] * cache.g̅[j][k]
            end

            b[D*(i-1)+k] = z
        end
    end

    # compute b = [qₙ - ϕ(0.5; qₙ⁻, qₙ⁺)]
    for k in 1:D
        b[D*S+k] = cache.q[k] - params.ρ⁻ * cache.q⁻[k] - params.ρ⁺ * cache.q⁺[k]
    end
end


"""
`IntegratorDGVIPI`: Discontinuous Galerkin Variational Integrator.

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
struct IntegratorDGVIPI{DT,TT,D,S,R,ΘT,FT,GT,VT,FPT,ST,IT,BT<:Basis,JT<:Discontinuity} <: DeterministicIntegrator{DT,TT}
    equation::IODE{DT,TT,ΘT,FT,GT,VT}

    basis::BT
    quadrature::Quadrature{TT,R}
    jump::JT

    Δt::TT

    params::FPT
    solver::ST
    iguess::InitialGuessPODE{DT,TT,VT,FT,IT}

    q::Vector{DT}
    q⁻::Vector{DT}
    q⁺::Vector{DT}

    cache::NonlinearFunctionCacheDGVIPI{DT}
end

function IntegratorDGVIPI(equation::IODE{DT,TT,ΘT,FT,GT,VT}, basis::Basis{TT,P},
                quadrature::Quadrature{TT,R}, jump::Discontinuity{TT,PT,QN}, Δt::TT;
                interpolation=HermiteInterpolation{DT}) where {DT,TT,ΘT,FT,GT,VT,P,R,PT,QN}

    D = equation.d
    S = nbasis(basis)

    N = D*(S+1)

    # create solution vector for nonlinear solver
    x = zeros(DT,N)

    # create solution vectors
    q  = zeros(DT,D)
    q⁻ = zeros(DT,D)
    q⁺ = zeros(DT,D)

    # create cache for internal stage vectors and update vectors
    cache = NonlinearFunctionCacheDGVIPI{DT,D,S,R,QN}()

    # create params
    params = ParametersDGVIPI(equation.α, equation.f, equation.g,
                Δt, basis, quadrature, jump, q, q⁻)

    # create rhs function for nonlinear solver
    function_stages = (x,b) -> function_stages!(x, b, params)

    # create nonlinear solver
    solver = get_config(:nls_solver)(x, function_stages)

    # create initial guess
    iguess = InitialGuessPODE(interpolation, equation, Δt)

    # create integrator
    IntegratorDGVIPI{DT, TT, D, S, R, ΘT, FT, GT, VT, typeof(params), typeof(solver),
                typeof(iguess.int), typeof(basis), typeof(jump)}(
                equation, basis, quadrature, jump, Δt, params, solver, iguess,
                q, q⁻, q⁺, cache)
end



function initialize!(int::IntegratorDGVIPI, sol::Union{SolutionPODE, SolutionPDAE}, m::Int)
    @assert m ≥ 1
    @assert m ≤ sol.ni

    # copy initial conditions from solution
    get_initial_conditions!(sol, int.q, int.q⁻, m)

    # initialise initial guess
    initialize!(int.iguess, m, sol.t[0], int.q, int.q⁻)
end


function update_solution!(int::IntegratorDGVIPI{DT,TT}, cache::NonlinearFunctionCacheDGVIPI{DT}) where {DT,TT}
    int.q  .= cache.q̅
    int.q⁻ .= cache.q̅⁻
    int.q⁺ .= cache.q̅⁺
end


@generated function initial_guess!(int::IntegratorDGVIPI{DT,TT,D,S,R}, m::Int) where {DT,TT,D,S,R}
    v = zeros(DT,D)
    y = zeros(DT,D)
    z = zeros(DT,D)

    quote
        # compute initial guess
        if nnodes(int.basis) > 0
            for i in 1:S
                evaluate!(int.iguess, m, $y, $z, $v, nodes(int.basis)[i], nodes(int.basis)[i])
                for k in 1:D
                    int.solver.x[D*(i-1)+k] = $y[k]
                end
            end
        else
            for i in 1:S
                for k in 1:D
                    int.solver.x[D*(i-1)+k] = 0
                end
            end
        end

        evaluate!(int.iguess, m, $y, $z, $v, one(TT), one(TT))
        for k in 1:D
            int.solver.x[D*S+k] = $y[k]
        end
    end
end


function integrate_step!(int::IntegratorDGVIPI{DT,TT}, sol::Union{SolutionPODE{DT,TT}, SolutionPDAE{DT,TT}}, m::Int, n::Int) where {DT,TT}
    # set time for nonlinear solver
    int.params.t = sol.t[0] + (n-1)*int.Δt

    # compute initial guess
    initial_guess!(int, m)

    # call nonlinear solver
    solve!(int.solver)

    # print solver status
    print_solver_status(int.solver.status, int.solver.params)

    # check if solution contains NaNs or error bounds are violated
    check_solver_status(int.solver.status, int.solver.params)

    # compute final update
    compute_stages!(int.solver.x, int.cache, int.params)

    # copy solution from cache to integrator
    update_solution!(int, int.cache)

    # copy solution to initial guess for next time step
    update!(int.iguess, m, sol.t[0] + n*int.Δt, int.q, int.q⁻)

    # take care of periodic solutions
    cut_periodic_solution!(int.q,  int.equation.periodicity)
    cut_periodic_solution!(int.q⁻, int.equation.periodicity)

    # copy to solution
    copy_solution!(sol, int.q, int.q⁻, n, m)
end
