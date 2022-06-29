@doc raw"""
`ParametersCGVI`: Parameters for right-hand side function of continuous Galerkin variational Integrator.

### Parameters

* `Θ`: function of the noncanonical one-form (∂L/∂v)
* `f`: function of the force (∂L/∂q)
* `Δt`: time step
* `b`: weights of the quadrature rule
* `c`: nodes of the quadrature rule
* `x`: nodes of the basis
* `m`: mass matrix
* `a`: derivative matrix
* `r₀`: reconstruction coefficients at the beginning of the interval
* `r₁`: reconstruction coefficients at the end of the interval
* `t`: initial time
* `q`: solution of q at time t
* `p`: solution of p at time t
"""
mutable struct ParametersCGVI{DT, TT, D, S, R, ET <: NamedTuple} <: Parameters{DT,TT}
    equs::ET
    Δt::TT

    b::Vector{TT}
    c::Vector{TT}

    x::Vector{TT}
    m::Matrix{TT}
    a::Matrix{TT}
    r₀::Vector{TT}
    r₁::Vector{TT}

    t::TT
    q::Vector{DT}
    p::Vector{DT}

    function ParametersCGVI{DT,D}(equs::ET, Δt::TT, b, c, x, m, a, r₀, r₁) where {DT, TT, D, ET <: NamedTuple}
        new{DT,TT,D,length(x),length(c),ET}(equs, Δt, b, c, x, m, a, r₀, r₁, zero(TT), zeros(DT,D), zeros(DT,D))
    end

    function ParametersCGVI{DT,D}(equs::NamedTuple, Δt::TT, basis::Basis{TT}, quadrature::QuadratureRule{TT}) where {DT,TT,D}
        # compute coefficients
        r₀ = zeros(TT, nbasis(basis))
        r₁ = zeros(TT, nbasis(basis))
        m  = zeros(TT, nnodes(quadrature), nbasis(basis))
        a  = zeros(TT, nnodes(quadrature), nbasis(basis))

        for i in eachindex(basis)
            r₀[i] = basis[zero(TT), i]
            r₁[i] = basis[one(TT), i]
            for j in eachindex(quadrature)
                m[j,i] = basis[nodes(quadrature)[j], i]
                a[j,i] = basis'[nodes(quadrature)[j], i]
            end
        end

        if get_config(:verbosity) > 1
            println()
            println("  Continuous Galerkin Variational Integrator")
            println("  ==========================================")
            println()
            println("    c  = ", nodes(quadrature))
            println("    b  = ", weights(quadrature))
            println("    x  = ", nodes(basis))
            println("    m  = ", m)
            println("    a  = ", a)
            println("    r₀ = ", r₀)
            println("    r₁ = ", r₁)
            println()
        end

        ParametersCGVI{DT,D}(equs, Δt, weights(quadrature), nodes(quadrature), grid(basis), m, a, r₀, r₁)
    end
end


function update_params!(params::ParametersCGVI, sol::AtomicSolutionPODE)
    # set time for nonlinear solver and copy previous solution
    params.t  = sol.t
    params.q .= sol.q
    params.p .= sol.p
end


struct IntegratorCacheCGVI{ST,D,S,R} <: IODEIntegratorCache{ST,D}
    q̃::Vector{ST}
    p̃::Vector{ST}
    ṽ::Vector{ST}
    f̃::Vector{ST}
    s̃::Vector{ST}

    X::Vector{Vector{ST}}
    Q::Vector{Vector{ST}}
    P::Vector{Vector{ST}}
    V::Vector{Vector{ST}}
    F::Vector{Vector{ST}}


    function IntegratorCacheCGVI{ST,D,S,R}() where {ST,D,S,R}
        # create temporary vectors
        q̃ = zeros(ST,D)
        p̃ = zeros(ST,D)
        ṽ = zeros(ST,D)
        f̃ = zeros(ST,D)
        s̃ = zeros(ST,D)

        # create internal stage vectors
        X = create_internal_stage_vector(ST,D,S)
        Q = create_internal_stage_vector(ST,D,R)
        P = create_internal_stage_vector(ST,D,R)
        V = create_internal_stage_vector(ST,D,R)
        F = create_internal_stage_vector(ST,D,R)

        new(q̃, p̃, ṽ, f̃, s̃, X, Q, P, V, F)
    end
end

function IntegratorCache{ST}(params::ParametersCGVI{DT,TT,D,S,R}; kwargs...) where {ST,DT,TT,D,S,R}
    IntegratorCacheCGVI{ST,D,S,R}(; kwargs...)
end

@inline CacheType(ST, params::ParametersCGVI{DT,TT,D,S,R}) where {DT,TT,D,S,R} = IntegratorCacheCGVI{ST,D,S,R}


"Continuous Galerkin Variational Integrator."
struct IntegratorCGVI{DT, TT, D, S, R, 
                      BT <: Basis,
                      PT <: ParametersCGVI{DT,TT,D,S,R},
                      ST <: NonlinearSolver{DT},
                      IT <: InitialGuessIODE{TT}} <: IODEIntegrator{DT,TT}
    basis::BT
    quadrature::QuadratureRule{TT,R}

    params::PT
    solver::ST
    iguess::IT
    caches::CacheDict{PT}

    function IntegratorCGVI(basis::BT, quadrature::QuadratureRule{TT,R}, params::ParametersCGVI{DT,TT,D,S},
                    solver::ST, iguess::IT, caches) where {DT,TT,D,S,R,BT,ST,IT}
        new{DT, TT, D, S, R, BT, typeof(params), ST, IT}(basis, quadrature, params, solver, iguess, caches)
    end

    function IntegratorCGVI{DT,D}(equations::NamedTuple, basis::Basis{TT}, quadrature::QuadratureRule{TT,R}, Δt::TT;
                                  interpolation=HermiteExtrapolation{DT}) where {DT,TT,D,R}

        # get number of stages
        S = nbasis(basis)

        # create params
        params = ParametersCGVI{DT,D}(equations, Δt, basis, quadrature)

        # create cache dict
        caches = CacheDict(params)

        # create nonlinear solver
        solver = create_nonlinear_solver(DT, D*(S+1), params, caches)

        # create initial guess
        iguess = InitialGuessIODE(get_config(:ig_extrapolation), equations[:v̄], equations[:f̄], Δt)

        # create integrator
        IntegratorCGVI(basis, quadrature, params, solver, iguess, caches)
    end

    function IntegratorCGVI(equation::Union{IODEProblem{DT}, LODEProblem{DT}}, basis::Basis, quadrature::QuadratureRule, Δt=tstep(equation); kwargs...) where {DT}
        IntegratorCGVI{DT, ndims(equation)}(functions(equation), basis, quadrature, Δt; kwargs...)
    end
end

@inline GeometricBase.equation(integrator::IntegratorCGVI, i::Symbol) = integrator.params.equs[i]
@inline GeometricBase.equations(integrator::IntegratorCGVI) = integrator.params.equs
@inline GeometricBase.timestep(integrator::IntegratorCGVI) = integrator.params.Δt
@inline Base.ndims(::IntegratorCGVI{DT,TT,D}) where {DT,TT,D} = D


function IntegratorCache(int::IntegratorCGVI{DT,TT}) where {DT,TT}
    IntegratorCacheCGVI{DT, TT, ndims(int), nbasis(int.basis), nnodes(int.quadrature)}()
end


function initialize!(int::IntegratorCGVI, sol::AtomicSolutionPODE)
    sol.t̄ = sol.t - timestep(int)

    equation(int, :v̄)(sol.t, sol.q, sol.v)
    equation(int, :f̄)(sol.t, sol.q, sol.v, sol.f)

    initialize!(int.iguess, sol.t, sol.q, sol.p, sol.v, sol.f,
                            sol.t̄, sol.q̄, sol.p̄, sol.v̄, sol.f̄)
end


function initial_guess!(int::IntegratorCGVI{DT,TT,D,S,R}, sol::AtomicSolutionPODE{DT,TT},
                        cache::IntegratorCacheCGVI{DT}=int.caches[DT]) where {DT,TT,D,S,R}
    int.solver.x .= 0

    for i in eachindex(int.basis)
        evaluate!(int.iguess, sol.q̄, sol.p̄, sol.v̄, sol.f̄,
                                sol.q, sol.p, sol.v, sol.f,
                                cache.q̃, cache.p̃,
                                cache.ṽ, cache.f̃,
                                grid(int.basis)[i], grid(int.basis)[i])

        for k in 1:D
            int.solver.x[D*(i-1)+k] = cache.q̃[k]
        end
    end

    evaluate!(int.iguess, sol.q̄, sol.p̄, sol.v̄, sol.f̄,
                          sol.q, sol.p, sol.v, sol.f,
                          cache.q̃, cache.p̃,
                          one(TT), one(TT))

    for k in 1:D
        int.solver.x[D*S+k] = cache.p̃[k]
    end
end


"Compute stages of variational partitioned Runge-Kutta methods."
function function_stages!(x::Vector{ST}, b::Vector{ST}, params::ParametersCGVI{DT,TT,D,S,R},
                caches::CacheDict) where {ST,DT,TT,D,S,R}
    @assert length(x) == length(b)

    # get cache for internal stages
    cache = caches[ST]

    # compute stages from nonlinear solver solution x
    compute_stages!(x, cache, params)

    # compute rhs b of nonlinear solver
    compute_rhs!(b, cache.X, cache.Q, cache.P, cache.F, cache.p̃, params)
end


function compute_stages!(x, cache::IntegratorCacheCGVI, params::ParametersCGVI)
    compute_stages!(x, cache.X, cache.Q, cache.V, cache.P, cache.F, cache.q̃, cache.p̃, params)
end

function compute_stages!(x, X, Q, V, P, F, q, p, params::ParametersCGVI{DT,TT,D,S,R}) where {DT,TT,D,S,R}

    # copy x to X
    for i in eachindex(X)
        for k in eachindex(X[i])
            X[i][k] = x[D*(i-1)+k]
        end
    end

    # copy x to p
    for k in eachindex(p)
        p[k] = x[D*S+k]
    end

    # compute Q
    compute_stages_q!(X, Q, q, params)

    # compute V
    compute_stages_v!(X, V, params)

    # compute P and F
    compute_stages_p!(Q, V, P, F, params)
end


function compute_stages_q!(X::Vector{Vector{ST}}, Q::Vector{Vector{ST}}, q::Vector{ST}, params::ParametersCGVI{DT,TT,D,S,R}) where {ST,DT,TT,D,S,R}
    @assert R == length(Q)
    @assert S == length(X)

    local y::ST

    # compute Q
    for i in eachindex(Q)
        @assert D == length(Q[i])
        for k in eachindex(Q[i])
            y = 0
            for j in eachindex(X)
                @assert D == length(X[j])
                y += params.m[i,j] * X[j][k]
            end
            Q[i][k] = y
        end
    end

    # compute q
    for k in eachindex(q)
        y = 0
        for i in eachindex(X)
            y += params.r₁[i] * X[i][k]
        end
        q[k] = y
    end
end


function compute_stages_v!(X::Vector{Vector{ST}}, V::Vector{Vector{ST}}, params::ParametersCGVI{DT,TT,D,S,R}) where {ST,DT,TT,D,S,R}
    @assert R == length(V)
    @assert S == length(X)

    local y::ST

    # compute V
    for i in eachindex(V)
        @assert D == length(V[i])
        for k in eachindex(V[i])
            y = 0
            for j in eachindex(X)
                @assert D == length(X[j])
                y += params.a[i,j] * X[j][k]
            end
            V[i][k] = y / params.Δt
        end
    end
end


function compute_stages_p!(Q::Vector{Vector{ST}}, V::Vector{Vector{ST}},
                           P::Vector{Vector{ST}}, F::Vector{Vector{ST}},
                           params::ParametersCGVI{DT,TT,D,S,R}) where {ST,DT,TT,D,S,R}

    @assert R == length(Q) == length(V) == length(P) == length(F)

    local tᵢ::TT

    # compute P=ϑ(Q,V) and F=f(Q,V)
    for i in eachindex(Q,V,P,F)
        @assert D == length(Q[i]) == length(V[i]) == length(P[i]) == length(F[i])
        tᵢ = params.t + params.Δt * params.c[i]
        params.equs[:ϑ](tᵢ, Q[i], V[i], P[i])
        params.equs[:f](tᵢ, Q[i], V[i], F[i])
    end
end


function compute_rhs!(b::Vector{ST}, X::Vector{Vector{ST}}, Q::Vector{Vector{ST}},
                                     P::Vector{Vector{ST}}, F::Vector{Vector{ST}},
                                     p::Vector{ST}, params::ParametersCGVI{DT,TT,D,S,R}) where {ST,DT,TT,D,S,R}
    local y::ST
    local z::ST

    # compute b = - [(P-AF)]
    for i in 1:S
        for k in 1:D
            z = 0
            for j in eachindex(P,F)
                z += params.b[j] * params.m[j,i] * F[j][k]
                z += params.b[j] * params.a[j,i] * P[j][k] / params.Δt
            end
            b[D*(i-1)+k] = z - (params.r₁[i] * p[k] - params.r₀[i] * params.p[k]) / params.Δt
        end
    end

    # compute b = - [(q-r₀Q)]
    for k in eachindex(params.q)
        y = 0
        for j in eachindex(X)
            y += params.r₀[j] * X[j][k]
        end
        b[D*S+k] = y - params.q[k]
    end
end


function update_solution!(sol::AtomicSolutionPODE, cache::IntegratorCacheCGVI)
    sol.q .= cache.q̃
    sol.p .= cache.p̃
end


function integrate_step!(int::IntegratorCGVI{DT,TT}, sol::AtomicSolutionPODE{DT,TT},
                         cache::IntegratorCacheCGVI{DT}=int.caches[DT]) where {DT,TT}
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

    # compute final update
    update_solution!(sol, cache)

    # copy solution to initial guess
    update_vector_fields!(int.iguess, sol.t, sol.q, sol.p, sol.v, sol.f)
end
