
"Holds the coefficients of a projected Gauss-Legendre Runge-Kutta method."
struct CoefficientsPGLRK{T} <: AbstractCoefficients{T}
    @HeaderCoefficientsRK
    @CoefficientsRK

    P::Matrix{T}
    Q::Matrix{T}
    X::Matrix{T}
    W::Matrix{T}
    A::Matrix{T}

    function CoefficientsPGLRK{T}(name,o,s,a,b,c,P,X,W) where {T}
        @assert T <: Real
        @assert isa(name, Symbol)
        @assert isa(o, Integer)
        @assert isa(s, Integer)
        @assert s ≥ 2 "Number of stages must be ≥ 2"
        @assert s==size(a,1)==size(a,2)==length(b)==length(c)
        @assert s==size(P,1)==size(P,2)
        @assert s==size(X,1)==size(X,2)
        @assert s==size(W,1)==size(W,2)

        Q = inv(P)
        A = zero(a)
        B = zero(a)

        simd_mult!(B, W, Q)
        simd_mult!(A, P, B)

        new(name,o,s,a,b,c,P,Q,X,W,A)
    end
end

function CoefficientsPGLRK(name::Symbol, order::Int, a::Matrix{T}, b::Vector{T}, c::Vector{T}, P::Matrix{T}, X::Matrix{T}, W::Matrix{T}) where {T}
    CoefficientsPGLRK{T}(name, order, length(c), a, b, c, P, X, W)
end

Base.hash(tab::CoefficientsPGLRK, h::UInt) = hash(tab.o, hash(tab.s, hash(tab.a, hash(tab.b, hash(tab.c, hash(tab.P, hash(tab.Q, hash(tab.X, hash(tab.W, hash(tab.A, h))))))))))

Base.:(==)(tab1::CoefficientsPGLRK, tab2::CoefficientsPGLRK) = (tab1.o == tab2.o
                                                             && tab1.s == tab2.s
                                                             && tab1.a == tab2.a
                                                             && tab1.b == tab2.b
                                                             && tab1.c == tab2.c
                                                             && tab1.P == tab2.P
                                                             && tab1.Q == tab2.Q
                                                             && tab1.X == tab2.X
                                                             && tab1.W == tab2.W
                                                             && tab1.A == tab2.A)

Base.isequal(tab1::CoefficientsPGLRK{T1}, tab2::CoefficientsPGLRK{T2}) where {T1,T2} = (tab1 == tab2 && T1 == T2 && typeof(tab1) == typeof(tab2))

"Print Runge-Kutta coefficients."
function Base.show(io::IO, tab::CoefficientsPGLRK)
    print(io, "Projected Gauss-Legendre Runge-Kutta Coefficients ", tab.name, " with ", tab.s, " stages and order ", tab.o)
    print(io, "  a = ", tab.a)
    print(io, "  b = ", tab.b)
    print(io, "  c = ", tab.c)
    print(io, "  P = ", tab.P)
    print(io, "  Q = ", tab.Q)
    print(io, "  X = ", tab.X)
    print(io, "  W = ", tab.W)
    print(io, "  A = ", tab.A)
end

function getTableauPGLRK(coeff::CoefficientsPGLRK{T}, λ, a::Matrix{T}) where {T}
    a .= coeff.a .+ λ .* coeff.A
end

function getTableauPGLRK(coeff::CoefficientsPGLRK{T}, λ) where {T}
    coeff.a .+ λ .* coeff.A
end


"Parameters for right-hand side function of projected Gauss-Legendre Runge-Kutta methods."
mutable struct ParametersPGLRK{DT,TT,D,S,ET} <: Parameters{DT,TT}
    equ::ET
    tab::CoefficientsPGLRK{TT}
    Δt::TT

    t̅::TT
    t::TT

    q̅::Vector{DT}

    Q::Vector{Vector{DT}}
    λ::DT

    function ParametersPGLRK{DT,TT,D,S,ET}(equ, tab, Δt) where {DT,TT,D,S,ET}
        new(equ, tab, Δt, zero(TT), zero(TT), zeros(DT,D),
            create_internal_stage_vector(DT,D,S), zero(DT))
    end
end


mutable struct IntegratorCachePGLRK{DT,D,S} <: ODEIntegratorCache{DT,D}
    λ::DT
    λ̅::DT
    h::DT
    h̅::DT

    q̃::Vector{DT}
    ṽ::Vector{DT}
    s̃::Vector{DT}

    Q::Vector{Vector{DT}}
    V::Vector{Vector{DT}}
    Y::Vector{Vector{DT}}

    function IntegratorCachePGLRK{DT,D,S}() where {DT,D,S}
        # create temporary vectors
        q̃ = zeros(DT,D)
        ṽ = zeros(DT,D)
        s̃ = zeros(DT,D)

        # create internal stage vectors
        Q = create_internal_stage_vector(DT, D, S)
        V = create_internal_stage_vector(DT, D, S)
        Y = create_internal_stage_vector(DT, D, S)


        new(0, 0, 0, 0,
            q̃, ṽ, s̃,
            Q, V, Y)
    end
end


"""
Projected Gauss-Legendre Runge-Kutta integrator.

    Reference: LUIGI BRUGNANO, FELICE IAVERNARO, AND DONATO TRIGIANTE.
        ENERGY- AND QUADRATIC INVARIANTS–PRESERVING INTEGRATORS BASED
        UPON GAUSS COLLOCATION FORMULAE.
        SIAM J. NUMER. ANAL. Vol. 50, No. 6, pp. 2897–2916, 2012.
"""
struct IntegratorPGLRK{DT, TT, PT <: ParametersPGLRK{DT,TT},
                               ST <: NonlinearSolver{DT},
                               PST <: NonlinearSolver{DT},
                               IT <: InitialGuessODE{DT,TT}, N, D, S} <: IntegratorPRK{DT,TT}
    params::PT
    solver::ST
    psolver::PST
    iguess::IT
    cache::IntegratorCachePGLRK{DT,D,S}

    x::Vector{DT}
    b::Vector{DT}

    function IntegratorPGLRK(N, params::ParametersPGLRK{DT,TT,D,S,ET}, solver::ST, psolver::PST, iguess::IT, cache,
                    x::Vector{DT}, b::Vector{DT}) where {DT, TT, D, S, ET, ST, PST, IT}
        new{DT, TT, typeof(params), ST, PST, IT, N, D, S}(params, solver, psolver, iguess, cache, x, b)
    end
end

function IntegratorPGLRK(equation::ODE{DT,TT,VT,HT,N}, tableau::CoefficientsPGLRK{TT}, Δt::TT) where {DT,TT,VT,HT,N}
    D = equation.d
    M = equation.n
    S = tableau.s

    # create params
    params = ParametersPGLRK{DT,TT,D,S,typeof(equation)}(equation, tableau, Δt)

    # create solvers
    solver  = create_nonlinear_solver(DT, D*S, params)
    psolver = create_nonlinear_solver(DT, 1,   params; F=function_projection!)

    # create full solution and RHS vector
    x = zeros(DT, D*S+1)
    b = zeros(DT, D*S+1)

    # create initial guess
    iguess = InitialGuessODE(get_config(:ig_interpolation), equation, Δt)

    # create cache
    cache = IntegratorCachePGLRK{DT,D,S}()

    # create integrator
    IntegratorPGLRK(N, params, solver, psolver, iguess, cache, x, b)
end


@inline nstages(integrator::IntegratorPGLRK{DT,TT,PT,ST,IT,N,D,S}) where {DT,TT,PT,ST,IT,N,D,S} = S


function update_params!(params::ParametersPGLRK, sol::AtomicSolutionODE)
    # set time for nonlinear solver and copy previous solution
    params.t̅  = sol.t
    params.t  = params.t̅ + params.Δt
    params.q̅ .= sol.q
end


function compute_stages!(x::Vector{ST}, cache::IntegratorCachePGLRK{ST}, params::ParametersPGLRK{DT,TT,D,S}) where {ST,DT,TT,D,S}
    compute_stages!(x, cache.Q, cache.V, cache.Y, cache.q̃, cache.ṽ, params)
end

function compute_stages!(x::Vector{ST}, Q::Vector{Vector{ST}}, V::Vector{Vector{ST}}, Y::Vector{Vector{ST}},
                                        q::Vector{ST}, y::Vector{ST},
                                        params::ParametersPGLRK{DT,TT,D,S}) where {ST,DT,TT,D,S}

    local tᵢ::TT

    # copy x to Y
    for i in 1:S
        for k in 1:D
            Y[i][k] = x[D*(i-1)+k]
        end
    end

    # compute Q=q̅+Δt*Y
    for i in 1:S
        for k in 1:D
            Q[i][k] = params.q̅[k] + params.Δt * Y[i][k]
        end
    end

    # compute V=v(T,Q)
    for i in 1:S
        tᵢ = params.t̅ + params.Δt * params.tab.c[i]
        params.equ.v(tᵢ, Q[i], V[i])
    end
end

"Compute stages of projected Gauss-Legendre Runge-Kutta methods."
@generated function function_stages!(x::Vector{ST}, b::Vector{ST}, params::ParametersPGLRK{DT,TT,D,S}) where {ST,DT,TT,D,S}

    cache = IntegratorCachePGLRK{ST,D,S}()

    quote
        compute_stages!(x, $cache, params)

        a = getTableauPGLRK(params.tab, params.λ)

        # compute b = [Y-AV]
        for i in 1:S
            for k in 1:D
                b[D*(i-1)+k] = $cache.Y[i][k]
                for j in 1:S
                    b[D*(i-1)+k] -= a[i,j] * $cache.V[j][k]
                end
            end
        end
    end
end


"Compute projection of projected Gauss-Legendre Runge-Kutta methods."
@generated function function_projection!(x::Vector{ST}, b::Vector{ST}, params::ParametersPGLRK{DT,TT,D,S}) where {ST,DT,TT,D,S}

    cache = IntegratorCachePGLRK{ST,D,S}()

    quote
        # copy x to λ
        $cache.λ = x[1]

        # compute full tableau
        a = getTableauPGLRK(params.tab, $cache.λ)
        ainv = inv(a)

        $cache.q̃ .= params.q̅
        for k in 1:D
            for i in 1:S
                for j in 1:S
                    $cache.q̃[k] += params.tab.b[i] * ainv[i,j] * (params.Q[j][k] - params.q̅[k])
                end
            end
        end

        # compute h and h̅
        $cache.h = params.equ.h(params.t, $cache.q̃)
        $cache.h̅ = params.equ.h(params.t̅, params.q̅)

        # compute b = [h̅-h]
        b[1] = $cache.h̅ - $cache.h
    end
end


@generated function function_all!(x::Vector{ST}, b::Vector{ST}, params::ParametersPGLRK{DT,TT,D,S}) where {ST,DT,TT,D,S}

    cache = IntegratorCachePGLRK{ST,D,S}()

    quote
        compute_stages!(x, $cache, params)

        $cache.λ = x[D*S+1]

        a = getTableauPGLRK(params.tab, $cache.λ)

        # compute b = [Y-AV]
        for i in 1:S
            for k in 1:D
                b[D*(i-1)+k] = $cache.Y[i][k]
                for j in 1:S
                    b[D*(i-1)+k] -= a[i,j] * $cache.V[j][k]
                end
            end
        end

        # compute y=B*V
        $cache.ṽ .= 0
        for k in 1:D
            for j in 1:S
                $cache.ṽ[k] += params.tab.b[j] * $cache.V[j][k]
            end
        end

        # compute q=q̅+Δt*y
        $cache.q̃ .= params.q̅ .+ params.Δt .* $cache.ṽ

        # compute h and h̅
        $cache.h = params.equ.h(params.t, $cache.q̃)
        $cache.h̅ = params.equ.h(params.t̅, params.q̅)

        # compute b = [h̅-h]
        b[D*S+1] = $cache.h̅ - $cache.h
    end
end


function compute_residuals(int::IntegratorPGLRK{DT,TT}) where {DT,TT}
    function_all!(int.x, int.b, int.params)

    rₐ²::DT = 0
    for bᵢ in int.b
        rₐ² = max(rₐ², bᵢ^2)
    end
    rₐ = sqrt(rₐ²)

    return rₐ
end


function initialize!(int::IntegratorPGLRK, sol::AtomicSolutionODE)
    sol.t̅ = sol.t - timestep(int)

    equation(int).v(sol.t, sol.q, sol.v)

    initialize!(int.iguess, sol.t, sol.q, sol.v,
                            sol.t̅, sol.q̅, sol.v̅)
end


function initial_guess!(int::IntegratorPGLRK{DT,TT}, sol::AtomicSolutionODE{DT,TT}) where {DT,TT}
    for i in eachstage(int)
        evaluate!(int.iguess, sol.q, sol.v, sol.q̅, sol.v̅,
                              int.cache.q̃, int.cache.ṽ,
                              tableau(int).c[i])

        for k in eachdim(int)
            int.solver.x[ndims(int)*(i-1)+k] = int.cache.ṽ[k]
        end
    end

    int.params.λ = 0
end


"Integrate ODE with projected Gauss-Legendre Runge-Kutta integrator."
function integrate_step!(int::IntegratorPGLRK{DT,TT}, sol::AtomicSolutionODE{DT,TT}) where {DT,TT}
    # update nonlinear solver parameters from cache
    update_params!(int.params, sol)

    # compute initial guess
    initial_guess!(int, sol)

    # reset cache
    reset!(sol, timestep(int))

    for i = 1:get_config(:nls_nmax)
        # call nonlinear solver
        solve!(int.solver)

        # print solver status
        print_solver_status(int.solver.status, int.solver.params)

        # check if solution contains NaNs or error bounds are violated
        check_solver_status(int.solver.status, int.solver.params)

        # compute vector fields at internal stages
        compute_stages!(int.solver.x, int.cache, int.params)

        # copy vector field at internal stages
        for i in eachstage(int)
            int.params.Q[i] .= int.cache.Q[i]
        end

        # call nonlinear solver for projection
        solve!(int.psolver)

        # print solver status
        print_solver_status(int.psolver.status, int.psolver.params)

        # check if solution contains NaNs or error bounds are violated
        check_solver_status(int.psolver.status, int.psolver.params)

        # copy parameter
        int.params.λ = int.psolver.x[1]

        # check convergence for full system
        int.x .= vcat(int.solver.x, int.psolver.x)
        rₐ = compute_residuals(int)

        # println(i, ", ", int.psolver.x, ", ", rₐ)

        rₐ ≤ get_config(:nls_atol) ? break : nothing
    end

    # compute vector fields at internal stages
    compute_stages!(int.solver.x, int.cache, int.params)

    # compute final update
    update_solution!(sol.q, sol.q̃, int.cache.V, tableau(int).b, timestep(int))

    # copy solution to initial guess
    update!(int.iguess, sol.t, sol.q, sol.v)
end
