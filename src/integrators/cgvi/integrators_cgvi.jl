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
* `t`: current time
* `q`: current solution of q
* `p`: current solution of p
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
end

function ParametersCGVI(Θ::ΘT, f::FT, Δt::TT, b, c, x, m, a, r₀, r₁, q::Vector{DT}, p::Vector{DT}) where {DT,TT,ΘT,FT}
    @assert length(q) == length(p)
    ParametersCGVI{DT,TT,length(q),length(x),length(c),ΘT,FT}(Θ, f, Δt, b, c, x, m, a, r₀, r₁, 0, q, p)
end


struct NonlinearFunctionCacheCGVI{ST}
    X::Vector{Vector{ST}}
    Q::Vector{Vector{ST}}
    V::Vector{Vector{ST}}
    P::Vector{Vector{ST}}
    F::Vector{Vector{ST}}

    q̅::Vector{ST}
    p̅::Vector{ST}

    function NonlinearFunctionCacheCGVI{ST}(D,S,R) where {ST}
        # create internal stage vectors
        X = create_internal_stage_vector(ST,D,S)
        Q = create_internal_stage_vector(ST,D,R)
        V = create_internal_stage_vector(ST,D,R)
        P = create_internal_stage_vector(ST,D,R)
        F = create_internal_stage_vector(ST,D,R)

        # create solution vectors
        q̅ = zeros(ST,D)
        p̅ = zeros(ST,D)

        new(X, Q, V, P, F, q̅, p̅)
    end
end


"Compute stages of variational partitioned Runge-Kutta methods."
@generated function function_stages!(x::Vector{ST}, b::Vector{ST}, params::ParametersCGVI{DT,TT,D,S,R}) where {ST,DT,TT,D,S,R}
    cache = NonlinearFunctionCacheCGVI{ST}(D, S, R)

    quote
        @assert length(x) == length(b)

        compute_stages!(x, $cache.X, $cache.Q, $cache.V, $cache.P, $cache.F, $cache.q̅, $cache.p̅, params)

        compute_rhs!(b, $cache.X, $cache.Q, $cache.P, $cache.F, $cache.p̅, params)
    end
end


function compute_stages!(x, X, Q, V, P, F, q̅, p̅, params::ParametersCGVI{DT,TT,D,S,R}) where {DT,TT,D,S,R}

    # copy x to X
    for i in eachindex(X)
        for k in eachindex(X[i])
            X[i][k] = x[D*(i-1)+k]
        end
    end

    # copy x to p̅
    for k in eachindex(p̅)
        p̅[k] = x[D*S+k]
    end

    # compute Q
    compute_stages_q!(X, Q, q̅, params)

    # compute V
    compute_stages_v!(X, V, params)

    # compute P and F
    compute_stages_p!(Q, V, P, F, params)
end


function compute_stages_q!(X::Vector{Vector{ST}}, Q::Vector{Vector{ST}}, q̅::Vector{ST}, params::ParametersCGVI{DT,TT,D,S,R}) where {ST,DT,TT,D,S,R}
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

    # compute q̅
    for k in eachindex(q̅)
        y = 0
        for i in eachindex(X)
            y += params.r₁[i] * X[i][k]
        end
        q̅[k] = y
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


function compute_rhs!(b::Vector{ST}, X::Vector{Vector{ST}}, Q::Vector{Vector{ST}}, P::Vector{Vector{ST}}, F::Vector{Vector{ST}}, p̅::Vector{ST}, params::ParametersCGVI{DT,TT,D,S,R}) where {ST,DT,TT,D,S,R}
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
            b[D*(i-1)+k] = z - (params.r₁[i] * p̅[k] - params.r₀[i] * params.p[k]) / params.Δt
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
struct IntegratorCGVI{DT,TT,ΘT,FT,GT,VT,FPT,ST,IT,BT<:Basis,D,S,R} <: DeterministicIntegrator{DT,TT}
    equation::IODE{DT,TT,ΘT,FT,GT,VT}

    basis::BT
    quadrature::Quadrature{TT,R}

    Δt::TT

    params::FPT
    solver::ST
    iguess::InitialGuessPODE{DT,TT,VT,FT,IT}

    q::Vector{DT}
    p::Vector{DT}

    cache::NonlinearFunctionCacheCGVI{DT}
end

function IntegratorCGVI(equation::IODE{DT,TT,ΘT,FT,GT,VT}, basis::Basis{TT,P}, quadrature::Quadrature{TT,R}, Δt::TT;
    interpolation=HermiteInterpolation{DT}) where {DT,TT,ΘT,FT,GT,VT,P,R}
    D = equation.d
    S = nbasis(basis)

    N = D*(S+1)

    # create solution vector for nonlinear solver
    x = zeros(DT,N)

    # create solution vectors
    q = zeros(DT,D)
    p = zeros(DT,D)

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


    # create cache for internal stage vectors and update vectors
    cache = NonlinearFunctionCacheCGVI{DT}(D,S,R)

    # create params
    params = ParametersCGVI(equation.α, equation.f, Δt, weights(quadrature), nodes(quadrature), nodes(basis), m, a, r₀, r₁, q, p)

    # create rhs function for nonlinear solver
    function_stages = (x,b) -> function_stages!(x, b, params)

    # create nonlinear solver
    solver = get_config(:nls_solver)(x, function_stages)

    # create initial guess
    iguess = InitialGuessPODE(interpolation, equation, Δt)

    # create integrator
    IntegratorCGVI{DT, TT, ΘT, FT, GT, VT, typeof(params), typeof(solver), typeof(iguess.int), typeof(basis), D, S, R}(
                equation, basis, quadrature, Δt, params, solver, iguess, q, p, cache)
end



function initialize!(int::IntegratorCGVI, sol::Union{SolutionPODE, SolutionPDAE}, m::Int)
    @assert m ≥ 1
    @assert m ≤ sol.ni

    # copy initial conditions from solution
    get_initial_conditions!(sol, int.q, int.p, m)

    # initialise initial guess
    initialize!(int.iguess, m, sol.t[0], int.q, int.p)
end


function update_solution!(int::IntegratorCGVI{DT,TT}, cache::NonlinearFunctionCacheCGVI{DT}) where {DT,TT}
    for k in eachindex(int.q, cache.q̅)
        int.q[k] = cache.q̅[k]
    end
    for k in eachindex(int.p, cache.p̅)
        int.p[k] = cache.p̅[k]
    end
end


@generated function initial_guess!(int::IntegratorCGVI{DT,TT,ΘT,FT,GT,VT,FPT,ST,IT,BT,D,S,R}, m::Int) where {DT,TT,ΘT,FT,GT,VT,FPT,ST,IT,BT,D,S,R}
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
            int.solver.x[D*S+k] = $z[k]
        end
    end
end


"Integrate ODE with variational partitioned Runge-Kutta integrator."
function integrate_step!(int::IntegratorCGVI{DT,TT}, sol::Union{SolutionPODE{DT,TT}, SolutionPDAE{DT,TT}}, m::Int, n::Int) where {DT,TT}
    # set time for nonlinear solver
    int.params.t = sol.t[0] + (n-1)*int.Δt

    # compute initial guess
    initial_guess!(int, m)

    # call nonlinear solver
    solve!(int.solver)

    # print solver status
    print_solver_status(int.solver.status, int.solver.params, n)

    # check if solution contains NaNs or error bounds are violated
    check_solver_status(int.solver.status, int.solver.params, n)

    # compute final update
    compute_stages!(int.solver.x, int.cache.X, int.cache.Q, int.cache.V, int.cache.P, int.cache.F, int.cache.q̅, int.cache.p̅, int.params)
    update_solution!(int, int.cache)

    # copy solution to initial guess for next time step
    update!(int.iguess, m, sol.t[0] + n*int.Δt, int.q, int.p)

    # take care of periodic solutions
    cut_periodic_solution!(int.q, int.equation.periodicity)

    # copy to solution
    copy_solution!(sol, int.q, int.p, n, m)
end
