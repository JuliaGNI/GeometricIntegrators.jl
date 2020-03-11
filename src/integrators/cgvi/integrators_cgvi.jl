"""
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
mutable struct ParametersCGVI{DT,TT,D,S,R,ΘT,FT} <: Parameters{DT,TT}
    Θ::ΘT
    f::FT

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

    function ParametersCGVI{DT,TT,D,S,R,ΘT,FT}(Θ::ΘT, f::FT, Δt::TT, b, c, x, m, a, r₀, r₁) where {DT,TT,D,S,R,ΘT,FT}
        new(Θ, f, Δt, b, c, x, m, a, r₀, r₁, zero(TT), zeros(DT,D), zeros(DT,D))
    end
end

function ParametersCGVI(equ::IODE{DT,TT,ΘT,FT}, Δt::TT, b, c, x, m, a, r₀, r₁) where {DT,TT,ΘT,FT}
    ParametersCGVI{DT,TT,ndims(equ),length(x),length(c),ΘT,FT}(equ.ϑ, equ.f, Δt, b, c, x, m, a, r₀, r₁)
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


function update_params!(params::ParametersCGVI, sol::AtomicSolutionPODE)
    # set time for nonlinear solver and copy previous solution
    params.t  = sol.t
    params.q .= sol.q
    params.p .= sol.p
end


"Compute stages of variational partitioned Runge-Kutta methods."
@generated function function_stages!(x::Vector{ST}, b::Vector{ST}, params::ParametersCGVI{DT,TT,D,S,R}) where {ST,DT,TT,D,S,R}
    cache = IntegratorCacheCGVI{ST,D,S,R}()

    quote
        @assert length(x) == length(b)

        compute_stages!(x, $cache, params)

        compute_rhs!(b, $cache.X, $cache.Q, $cache.P, $cache.F, $cache.p̃, params)
    end
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
        params.Θ(tᵢ, Q[i], V[i], P[i])
        params.f(tᵢ, Q[i], V[i], F[i])
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


"Continuous Galerkin Variational Integrator."
struct IntegratorCGVI{DT,TT,D,S,R,ΘT,FT,GT,HT,VT,PT,ST,IT,BT<:Basis} <: DeterministicIntegrator{DT,TT}
    equation::IODE{DT,TT,ΘT,FT,GT,HT,VT}

    basis::BT
    quadrature::Quadrature{TT,R}

    Δt::TT

    params::PT
    solver::ST
    iguess::InitialGuessPODE{DT,TT,VT,FT,IT}
    cache::IntegratorCacheCGVI{DT,D,S,R}
end

function IntegratorCGVI(equation::IODE{DT,TT,ΘT,FT,GT,HT,VT}, basis::Basis{TT,P}, quadrature::Quadrature{TT,R}, Δt::TT;
    interpolation=HermiteInterpolation{DT}) where {DT,TT,ΘT,FT,GT,HT,VT,P,R}
    D = equation.d
    S = nbasis(basis)

    N = D*(S+1)

    # create solution vector for nonlinear solver
    x = zeros(DT,N)

    # compute coefficients
    r₀ = zeros(TT, nbasis(basis))
    r₁ = zeros(TT, nbasis(basis))
    m  = zeros(TT, nnodes(quadrature), S)
    a  = zeros(TT, nnodes(quadrature), S)

    for i in 1:nbasis(basis)
        r₀[i] = evaluate(basis, i, zero(TT))
        r₁[i] = evaluate(basis, i, one(TT))
        for j in 1:nnodes(quadrature)
            m[j,i] = evaluate(basis, i, nodes(quadrature)[j])
            a[j,i] = derivative(basis, i, nodes(quadrature)[j])
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


    # create cache for internal stage vectors and update vectors
    # cache = NonlinearFunctionCacheCGVI{DT}(D,S,R)

    # create params
    params = ParametersCGVI(equation, Δt, weights(quadrature), nodes(quadrature), nodes(basis), m, a, r₀, r₁)

    # create rhs function for nonlinear solver
    function_stages = (x,b) -> function_stages!(x, b, params)

    # create nonlinear solver
    solver = get_config(:nls_solver)(x, function_stages)

    # create initial guess
    iguess = InitialGuessPODE(interpolation, equation, Δt)

    # create cache
    cache = IntegratorCacheCGVI{DT,D,S,R}()

    # create integrator
    IntegratorCGVI{DT, TT, D, S, R, ΘT, FT, GT, HT, VT, typeof(params), typeof(solver), typeof(iguess.int), typeof(basis)}(
                equation, basis, quadrature, Δt, params, solver, iguess, cache)
end

equation(integrator::IntegratorCGVI) = integrator.equation
timestep(integrator::IntegratorCGVI) = integrator.Δt


function IntegratorCache(int::IntegratorCGVI{DT,TT}) where {DT,TT}
    IntegratorCacheCGVI{DT, TT, ndims(int), nbasis(int.basis), nnodes(int.quadrature)}()
end


function initialize!(int::IntegratorCGVI, sol::AtomicSolutionPODE)
    sol.t̅ = sol.t - timestep(int)

    equation(int).v(sol.t, sol.q, sol.p, sol.v)
    equation(int).f(sol.t, sol.q, sol.p, sol.f)

    initialize!(int.iguess, sol.t, sol.q, sol.p, sol.v, sol.f,
                            sol.t̅, sol.q̅, sol.p̅, sol.v̅, sol.f̅)
end


function initial_guess!(int::IntegratorCGVI{DT,TT,D,S,R}, sol::AtomicSolutionPODE{DT,TT}) where {DT,TT,D,S,R}
    if nnodes(int.basis) > 0
        for i in 1:S
            evaluate!(int.iguess, sol.q, sol.p, sol.v, sol.f,
                                  sol.q̅, sol.p̅, sol.v̅, sol.f̅,
                                  int.cache.q̃, int.cache.p̃,
                                  int.cache.ṽ, int.cache.f̃,
                                  nodes(int.basis)[i], nodes(int.basis)[i])

            for k in 1:D
                int.solver.x[D*(i-1)+k] = int.cache.q̃[k]
            end
        end
    else
        for i in 1:S
            for k in 1:D
                int.solver.x[D*(i-1)+k] = 0
            end
        end
    end

    evaluate!(int.iguess, sol.q, sol.p, sol.v, sol.f,
                          sol.q̅, sol.p̅, sol.v̅, sol.f̅,
                          int.cache.q̃, int.cache.p̃,
                          one(TT), one(TT))

    for k in 1:D
        int.solver.x[D*S+k] = int.cache.p̃[k]
    end
end


function update_solution!(sol::AtomicSolutionPODE, cache::IntegratorCacheCGVI)
    sol.q .= cache.q̃
    sol.p .= cache.p̃
end


"Integrate ODE with variational partitioned Runge-Kutta integrator."
function integrate_step!(int::IntegratorCGVI{DT,TT}, sol::AtomicSolutionPODE{DT,TT}) where {DT,TT}
    # update nonlinear solver parameters from cache
    update_params!(int.params, sol)

    # compute initial guess
    initial_guess!(int, sol)

    # reset cache
    reset!(sol, timestep(int))

    # call nonlinear solver
    solve!(int.solver)

    # print solver status
    print_solver_status(int.solver.status, int.solver.params)

    # check if solution contains NaNs or error bounds are violated
    check_solver_status(int.solver.status, int.solver.params)

    # compute vector fields at internal stages
    compute_stages!(int.solver.x, int.cache, int.params)

    # compute final update
    update_solution!(sol, int.cache)

    # copy solution to initial guess
    update_vector_fields!(int.iguess, sol.t, sol.q, sol.p, sol.v, sol.f)
end
