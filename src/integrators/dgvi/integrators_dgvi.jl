@doc raw"""
`ParametersDGVI`: Parameters for right-hand side function of Discontinuous Galerkin Variational Integrator.

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
* `t`:  initial time
* `q`:  solution of q  at time t
* `q⁻`: solution of q⁻ at time t
* `q⁺`: solution of q⁺ at time t
"""
mutable struct ParametersDGVI{DT, TT, D, S, R, ET <: NamedTuple} <: Parameters{DT,TT}
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
    p::Vector{DT}

    q⁻::Vector{DT}
    q⁺::Vector{DT}

    function ParametersDGVI{DT,D}(equs::ET, Δt::TT, b, c, m, a, r⁻, r⁺) where {DT, TT, D, ET <: NamedTuple}
        @assert length(b)  == length(c)
        @assert length(r⁻) == length(r⁺)
        new{DT,TT,D,length(r⁻),length(c),ET}(equs, Δt, b, c, m, a, r⁻, r⁺, zero(TT), zeros(DT,D), zeros(DT,D), zeros(DT,D), zeros(DT,D))
    end

    function ParametersDGVI{DT,D}(equs::NamedTuple, Δt::TT, basis::Basis{TT}, quadrature::QuadratureRule{TT}) where {DT,TT,D}
        # compute coefficients
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

        ParametersDGVI{DT,D}(equs, Δt, weights(quadrature), nodes(quadrature), m, a, r⁻, r⁺)
    end
end


function update_params!(params::ParametersDGVI, sol::AtomicSolutionPODE)
    # set time for nonlinear solver and copy previous solution
    params.t  = sol.t
    params.q .= sol.q
    params.p .= sol.p
end


@doc raw"""
Nonlinear function cache for Discontinuous Galerkin Variational Integrator.

### Parameters

* `ST`: data type
* `D`: number of dimensions
* `S`: number of degrees of freedom
* `R`: number of nodes of quadrature formula

### Fields

* `X`: degrees of freedom
* `Q`: solution at quadrature nodes
* `V`: velocity at quadrature nodes
* `P`: one-form at quadrature nodes
* `F`: forces at quadrature nodes
* `q`:  current solution of ``q_{n}``
* `q⁻`: current solution of ``q_{n}^{-}``
* `q⁺`: current solution of ``q_{n}^{+}``
* `q̅`:  current solution of ``q_{n+1}``
* `q̅⁻`: current solution of ``q_{n+1}^{-}``
* `q̅⁺`: current solution of ``q_{n+1}^{+}``
* `ϕ`:  average of the solution at ``t_{n}``
* `ϕ̅`:  average of the solution at ``t_{n+1}``
* `λ`:  jump of the solution at ``t_{n}``
* `λ̅`:  jump of the solution at ``t_{n+1}``
* `θ`:  one-form evaluated across at ``t_{n}``
* `Θ̅`:  one-form evaluated across at ``t_{n+1}``
* `g`:  projection evaluated across at ``t_{n}``
* `g̅`:  projection evaluated across at ``t_{n+1}``
"""
struct IntegratorCacheDGVI{ST,D,S,R} <: IODEIntegratorCache{ST,D}
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
    q̅::Vector{ST}
    q̅⁻::Vector{ST}
    q̅⁺::Vector{ST}

    p::Vector{ST}
    p̅::Vector{ST}

    ϕ::Vector{ST}
    ϕ⁻::Vector{ST}
    ϕ⁺::Vector{ST}
    ϕ̅::Vector{ST}
    ϕ̅⁻::Vector{ST}
    ϕ̅⁺::Vector{ST}

    λ::Vector{ST}
    λ⁻::Vector{ST}
    λ⁺::Vector{ST}
    λ̅::Vector{ST}
    λ̅⁻::Vector{ST}
    λ̅⁺::Vector{ST}

    θ::Vector{ST}
    θ⁻::Vector{ST}
    θ⁺::Vector{ST}
    Θ̅::Vector{ST}
    Θ̅⁻::Vector{ST}
    Θ̅⁺::Vector{ST}

    g::Vector{ST}
    g⁻::Vector{ST}
    g⁺::Vector{ST}
    g̅::Vector{ST}
    g̅⁻::Vector{ST}
    g̅⁺::Vector{ST}

    h⁺::Vector{ST}
    h̅⁻::Vector{ST}

    function IntegratorCacheDGVI{ST,D,S,R}() where {ST,D,S,R}
        # create internal stage vectors
        X = create_internal_stage_vector(ST,D,S)
        Q = create_internal_stage_vector(ST,D,R)
        V = create_internal_stage_vector(ST,D,R)
        P = create_internal_stage_vector(ST,D,R)
        F = create_internal_stage_vector(ST,D,R)

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
        q̅  = zeros(ST,D)
        q̅⁻ = zeros(ST,D)
        q̅⁺ = zeros(ST,D)

        p = zeros(ST,D)
        p̅ = zeros(ST,D)

        # create jump vectors
        ϕ  = zeros(ST,D)
        ϕ⁻ = zeros(ST,D)
        ϕ⁺ = zeros(ST,D)
        ϕ̅  = zeros(ST,D)
        ϕ̅⁻ = zeros(ST,D)
        ϕ̅⁺ = zeros(ST,D)

        λ  = zeros(ST,D)
        λ⁻ = zeros(ST,D)
        λ⁺ = zeros(ST,D)
        λ̅  = zeros(ST,D)
        λ̅⁻ = zeros(ST,D)
        λ̅⁺ = zeros(ST,D)

        θ  = zeros(ST,D)
        θ⁻ = zeros(ST,D)
        θ⁺ = zeros(ST,D)
        Θ̅  = zeros(ST,D)
        Θ̅⁻ = zeros(ST,D)
        Θ̅⁺ = zeros(ST,D)

        g  = zeros(ST,D)
        g⁻ = zeros(ST,D)
        g⁺ = zeros(ST,D)
        g̅  = zeros(ST,D)
        g̅⁻ = zeros(ST,D)
        g̅⁺ = zeros(ST,D)

        h⁺ = zeros(ST,D)
        h̅⁻ = zeros(ST,D)

        new(X, Q, V, P, F,
            q̃, p̃, ṽ, f̃, s̃,
            q, q⁻, q⁺, q̅, q̅⁻, q̅⁺,
            p, p̅,
            ϕ, ϕ⁻, ϕ⁺, ϕ̅, ϕ̅⁻, ϕ̅⁺,
            λ, λ⁻, λ⁺, λ̅, λ̅⁻, λ̅⁺,
            θ, θ⁻, θ⁺, Θ̅, Θ̅⁻, Θ̅⁺,
            g, g⁻, g⁺, g̅, g̅⁻, g̅⁺,
            h⁺, h̅⁻)
    end
end

function IntegratorCache{ST}(params::ParametersDGVI{DT,TT,D,S,R}; kwargs...) where {ST,DT,TT,D,S,R}
    IntegratorCacheDGVI{ST,D,S,R}(; kwargs...)
end

@inline CacheType(ST, params::ParametersDGVI{DT,TT,D,S,R}) where {DT,TT,D,S,R} = IntegratorCacheDGVI{ST,D,S,R}


@doc raw"""
`IntegratorDGVI`: Discontinuous Galerkin Variational Integrator.


The DGVI integrators arise from the discretization of the action integral
```math
\mathcal{A} [q] = \int \limits_{0}^{T} L(q(t), \dot{q}(t)) \, dt ,
```
with ``L`` a fully degenerate Lagrangian of the form
```math
L(q, \dot{q}) = \vartheta (q) \cdot \dot{q} - H(q) ,
```
where ``\vartheta (q)`` denotes the Cartan one-form and ``H(q)`` the Hamiltonian,
which is usually given by the total energy of the system.


### Discretization

Within each interval ``(t_{n}, t_{n+1})`` a piecewise-polynomial approximation ``q_h``
of the trajectory ``q`` is constructed using ``S`` basis functions ``\varphi_{i}``,
```math
q_h(t) \vert_{(t_{n}, t_{n+1})} = \sum \limits_{i=1}^{S} x_{n,i} \, \bar{\varphi}_{n,i} (t) ,
```
where ``\bar{\varphi}_{n,i} (t)`` is a rescaled basis function, defined by
```math
\bar{\varphi}_{n,i} (t) = \varphi_{i} \bigg( \frac{t - t_{n}}{t_{n+1} - t_{n}} \bigg) ,
```
and it is assumed that ``\varphi_{i}`` is compactly supported on ``[0,1]``.
These approximations ``q_h(t)`` are not assumed to be continuous across interval
boundaries ``t_{n}`` but usuaslly have jumps.

The integral over ``(t_{n}, t_{n+1})`` is approximated by a quadrature rule with
``R`` nodes ``c_i`` and weights ``b_i``.
Denote by ``m`` and ``a`` mass and derivative matrices, respectively, whose elements
 are given by
```math
m_{ij} = \varphi_j (c_i) ,
\qquad
a_{ij} = \varphi_j' (c_i) ,
\qquad
i = 1, ..., R ,
\;
j = 1, ..., S .
```
With that, the solution and its time derivative at the quadrature points can be written as
```math
Q_{n,i} \equiv q_h(t_n + c_i h) = m_{ij} x_{n,j} ,
\qquad
V_{n,i} \equiv \dot{q}_h (t_n + c_i h) = a_{ij} x_{n,j} ,
```
where
```math
x_{n} = ( x_{n,1}, ..., x_{n,S} )^T
```
is the vector containing the degrees of freedom of ``q_h \vert_{[t_{n}, t_{n+1}]}``.
The limits of ``q_h(t)`` at ``t_{n}`` and ``t_{n+1}`` are given by
```math
q_{n}^{+} = \lim \limits_{t \downarrow t_{n}} q_h(t) = \sum \limits_{j=1}^{S} r^{+}_{j} \delta x_{n,j} ,
\hspace{3em}
q_{n+1}^{-} = \lim \limits_{t \uparrow t_{n+1}} q_h(t) = \sum \limits_{j=1}^{S} r^{-}_{j} \delta x_{n,j} .
```

The discrete action reads
```math
\mathcal{A}_d [x_d] = h \sum \limits_{n=0}^{N-1} \bigg[
     \sum \limits_{i=1}^{R} b_i \big[ \vartheta (Q_{n,i}) \cdot V_{n,i} - H(Q_{n,i}) \big] \\
     + \frac{\vartheta (q_n) + \vartheta (q_n^+)}{2} \cdot (q_n^+ - q_n)
     + \frac{\vartheta (q_{n+1}^-) + \vartheta (q_{n+1})}{2} \cdot (q_{n+1} - q_{n+1}^-)
\bigg] ,
```
so that using the relations
```math
\delta Q_{n,i} = m_{ij} \delta x_{n,j} ,
\qquad
\delta V_{n,i} = \frac{a_{ij}}{h} \delta x_{n,j} ,
\qquad
\delta q_{n}^- = r^{-}_{j} \delta x_{n-1,j} ,
\qquad
\delta q_{n}^+ = r^{+}_{j} \delta x_{n,j} ,
```
the discrete action principle leads to the discrete equations of motion,
```math
0 = \sum \limits_{i=1}^{R} b_i \big[ h m_{ij} \nabla \vartheta (Q_{n,i}) \cdot V_{n,i} + a_{ij} \vartheta (Q_{n,i}) - h m_{ij} \nabla H(Q_{n,i}) \big] \\
+ r^{+}_{j} \, \frac{ \vartheta ( q_{n  }   ) + \vartheta( q_{n  }^+ ) }{2}
- r^{-}_{j} \, \frac{ \vartheta ( q_{n+1}^- ) + \vartheta( q_{n+1}   ) }{2}
   \big] \\
+ h r^{+}_{j} \, \nabla \vartheta (q_{n  }^+) \cdot (q_{n  }^+ - q_{n  }  )
+ h r^{-}_{j} \, \nabla \vartheta (q_{n+1}^-) \cdot (q_{n+1}   - q_{n+1}^-) ,
```
and
```math
\vartheta(q_{n}^+) = \vartheta (q_{n}^-) + \nabla \vartheta (q_{n}) \cdot (q_{n}^+ - q_{n}^-) ,
```
for all ``n`` and all ``j``.
Let us introduce the variable ``p_n`` as
```math
p_{n} = \vartheta (q_{n}^-) + \nabla \vartheta (q_{n}) \cdot (q_{n} - q_{n}^-) ,
```
so that
```math
\vartheta(q_{n}^+) = p_{n} + \nabla \vartheta (q_{n}) \cdot (q_{n}^+ - q_{n}) .
```
Then the above equations provide a map ``(q_{n}, p_{n}) \mapsto (q_{n+1}, p_{n+1})``.
In order to solve these equations, initial conditions ``q_{0}`` and
``p_{0} = \vartheta(q_{0})`` have to be prescribed.


### Fields

* `equation`: Implicit Ordinary Differential Equation
* `basis`: piecewise polynomial basis
* `quadrature`: numerical quadrature rule
* `Δt`: time step
* `params`: ParametersDGVI
* `solver`: nonlinear solver
* `iguess`: initial guess
* `q`: current solution vector for trajectory
* `p`: current solution vector for one-form
* `cache`: temporary variables for nonlinear solver
"""
struct IntegratorDGVI{DT, TT, D, S, R,
                BT <: Basis,
                PT <: ParametersDGVI{DT,TT,D,S},
                ST <: NonlinearSolver{DT},
                IT <: InitialGuessODE{TT}} <: IODEIntegrator{DT,TT}
    basis::BT
    quadrature::QuadratureRule{TT,R}

    params::PT
    solver::ST
    iguess::IT
    caches::CacheDict{PT}

    function IntegratorDGVI(basis::BT, quadrature::QuadratureRule{TT,R}, params::ParametersDGVI{DT,TT,D,S},
                    solver::ST, iguess::IT, caches) where {DT,TT,D,S,R,BT,ST,IT}
        new{DT, TT, D, S, R, BT, typeof(params), ST, IT}(basis, quadrature, params, solver, iguess, caches)
    end

    function IntegratorDGVI{DT,D}(equations::NamedTuple, basis::Basis{TT}, quadrature::QuadratureRule{TT,R}, Δt::TT;
                                  interpolation=HermiteInterpolation{DT}) where {DT,TT,D,R}

        # get number of stages
        S = nbasis(basis)

        # create params
        params = ParametersDGVI{DT,D}(equations, Δt, basis, quadrature)

        # create cache dict
        caches = CacheDict(params)

        # create nonlinear solver
        solver = create_nonlinear_solver(DT, D*(S+1), params, caches)

        # create initial guess
        iguess = InitialGuessODE(get_config(:ig_interpolation), equations[:v̄], Δt)

        # create integrator
        IntegratorDGVI(basis, quadrature, params, solver, iguess, caches)
    end

    function IntegratorDGVI(equation::IODE{DT,TT}, basis::Basis{TT}, quadrature::QuadratureRule{TT}, Δt::TT; kwargs...) where {DT,TT}
        IntegratorDGVI{DT, ndims(equation)}(get_function_tuple(equation), basis, quadrature, Δt; kwargs...)
    end
end

@inline equation(integrator::IntegratorDGVI, i::Symbol) = integrator.params.equs[i]
@inline equations(integrator::IntegratorDGVI) = integrator.params.equs
@inline timestep(integrator::IntegratorDGVI) = integrator.params.Δt
@inline Base.ndims(::IntegratorDGVI{DT,TT,D}) where {DT,TT,D} = D


function IntegratorCache(int::IntegratorDGVI{DT,TT}) where {DT,TT}
    IntegratorCacheDGVI{DT, TT, ndims(int), nbasis(int.basis), nnodes(int.quadrature)}()
end


function initialize!(int::IntegratorDGVI, sol::AtomicSolutionPODE)
    sol.t̅ = sol.t - timestep(int)

    equation(int, :v̄)(sol.t, sol.q, sol.v)

    initialize!(int.iguess, sol.t, sol.q, sol.v,
                            sol.t̅, sol.q̅, sol.v̅)
end


function initial_guess!(int::IntegratorDGVI{DT,TT,D,S,R}, sol::AtomicSolutionPODE{DT,TT},
                        cache::IntegratorCacheDGVI{DT}=int.caches[DT]) where {DT,TT,D,S,R}
    if nbasis(int.basis) > 0
        for i in eachindex(int.basis)
            evaluate!(int.iguess, sol.q̅, sol.v̅,
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

    evaluate!(int.iguess, sol.q̅, sol.v̅,
                          sol.q, sol.v,
                          cache.q̃,
                          one(TT))

    for k in 1:D
        int.solver.x[D*S+k] = cache.q̃[k]
    end
end


"Compute stages of variational partitioned Runge-Kutta methods."
function function_stages!(x::Vector{ST}, b::Vector{ST}, params::ParametersDGVI{DT,TT,D,S,R},
                caches::CacheDict) where {ST,DT,TT,D,S,R}
    @assert length(x) == length(b)

    # get cache for internal stages
    cache = caches[ST]

    # compute stages from nonlinear solver solution x
    compute_stages!(x, cache, params)

    # compute rhs b of nonlinear solver
    compute_rhs!(b, cache, params)
end


function compute_stages!(x, cache::IntegratorCacheDGVI{ST,D,S}, params::ParametersDGVI{DT,TT,D,S}) where {ST,DT,TT,D,S}
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
                           params::ParametersDGVI{DT,TT,D,S,R}) where {ST,DT,TT,D,S,R}

    local q::ST
    local q⁺::ST
    local q̅⁻::ST

    local X = cache.X
    local Q = cache.Q

    # copy q and q⁻
    cache.q  .= params.q
    cache.p  .= params.p
    cache.q⁻ .= params.q⁻

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
end


"Compute velocities at quadrature nodes."
function compute_stages_v!(cache::IntegratorCacheDGVI{ST,D,S,R},
                           params::ParametersDGVI{DT,TT,D,S,R}) where {ST,DT,TT,D,S,R}
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
                           params::ParametersDGVI{DT,TT,D,S,R}) where {ST,DT,TT,D,S,R}

    local tᵢ::TT

    # compute P=ϑ(Q) and F=f(Q)
    for i in 1:R
        tᵢ = params.t + params.Δt * params.c[i]
        params.equs[:ϑ](tᵢ, cache.Q[i], cache.V[i], cache.P[i])
        params.equs[:f](tᵢ, cache.Q[i], cache.V[i], cache.F[i])
    end
end


function compute_stages_λ!(cache::IntegratorCacheDGVI{ST,D,S,R},
                           params::ParametersDGVI{DT,TT,D,S,R}) where {ST,DT,TT,D,S,R}

    local t₀::TT = params.t
    local t₁::TT = params.t + params.Δt

    # compute ϕ and ϕ̅
    cache.ϕ  .= 0.5 * (cache.q⁺ .+ cache.q⁻)
    cache.ϕ⁻ .= 0.5 * (cache.q⁻ .+ cache.q )
    cache.ϕ⁺ .= 0.5 * (cache.q  .+ cache.q⁺)
    cache.ϕ̅⁻ .= 0.5 * (cache.q̅⁻ .+ cache.q̅ )

    # compute λ and λ̅
    cache.λ  .= cache.q⁺ .- cache.q⁻
    cache.λ⁻ .= cache.q  .- cache.q⁻
    cache.λ⁺ .= cache.q⁺ .- cache.q
    cache.λ̅⁻ .= cache.q̅  .- cache.q̅⁻

    # compute ϑ
    params.equs[:ϑ](t₀, cache.q,  cache.q,  cache.θ)
    params.equs[:ϑ](t₀, cache.q⁻, cache.q⁻, cache.θ⁻)
    params.equs[:ϑ](t₀, cache.q⁺, cache.q⁺, cache.θ⁺)

    params.equs[:ϑ](t₁, cache.q̅,  cache.q̅,  cache.Θ̅)
    params.equs[:ϑ](t₁, cache.q̅⁻, cache.q̅⁻, cache.Θ̅⁻)

    # compute projection
    params.equs[:g](t₀, cache.q,  cache.λ,  cache.g)
    params.equs[:g](t₀, cache.q⁻, cache.λ⁻, cache.g⁻)
    params.equs[:g](t₀, cache.q⁺, cache.λ⁺, cache.g⁺)
    params.equs[:g](t₁, cache.q̅⁻, cache.λ̅⁻, cache.g̅⁻)
    params.equs[:g](t₀, cache.q,  cache.λ⁺, cache.h⁺)
    params.equs[:g](t₁, cache.q̅,  cache.λ̅⁻, cache.h̅⁻)

    # # compute ϑ
    # params.equs[:ϑ](t₀, cache.ϕ⁻, cache.ϕ⁻, cache.θ⁻)
    # params.equs[:ϑ](t₀, cache.ϕ⁺, cache.ϕ⁺, cache.θ⁺)
    # params.equs[:ϑ](t₁, cache.ϕ̅⁻, cache.ϕ̅⁻, cache.Θ̅⁻)
    #
    # # compute projection
    # params.equs[:g](t₀, cache.ϕ⁻, cache.λ⁻, cache.g⁻)
    # params.equs[:g](t₀, cache.ϕ⁺, cache.λ⁺, cache.g⁺)
    # params.equs[:g](t₁, cache.ϕ̅⁻, cache.λ̅⁻, cache.g̅⁻)

    # compute pₙ₊₁ = ϑ(qₙ₊₁⁻) + ∇ϑ(qₙ₊₁)⋅(qₙ₊₁-qₙ₊₁⁻)
    #    i.e. p̅ = ϑ(q̅⁻) + ∇ϑ(q̅)⋅(q̅-q̅⁻)
    cache.p̅ .= cache.Θ̅⁻ .+ cache.h̅⁻
end


function compute_rhs!(b::Vector{ST}, cache::IntegratorCacheDGVI{ST,D,S,R},
                params::ParametersDGVI{DT,TT,D,S,R}) where {ST,DT,TT,D,S,R}

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
            z += params.r⁻[i] * 0.5 * cache.g̅⁻[k]

            # z += params.r⁺[i] * cache.θ⁺[k]
            # z -= params.r⁻[i] * cache.Θ̅⁻[k]
            #
            # z += params.r⁺[i] * cache.g⁺[k]
            # z += params.r⁻[i] * cache.g̅⁻[k]

            b[D*(i-1)+k] = z
        end
    end

    # compute b = ϑ(qₙ⁺) - pₙ - ∇ϑ(qₙ)⋅(qₙ⁺-qₙ)
    for k in 1:D
        b[D*S+k] = cache.θ⁺[k] - cache.p[k] - cache.h⁺[k]
    end
end


function update_solution!(sol::AtomicSolutionPODE, cache::IntegratorCacheDGVI)
    sol.q .= cache.q̅
    sol.p .= cache.p̅
end


function integrate_step!(int::IntegratorDGVI{DT,TT}, sol::AtomicSolutionPODE{DT,TT},
                         cache::IntegratorCacheDGVI{DT}=int.caches[DT]) where {DT,TT}
    # update nonlinear solver parameters from cache
    update_params!(int.params, sol)

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
    update_solution!(sol, cache)

    # copy solution to initial guess
    update_vector_fields!(int.iguess, sol.t, sol.q, sol.v)
end
