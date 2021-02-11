
using QuadratureRules

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

        mul!(B, W, Q)
        mul!(A, P, B)

        new(name,o,s,a,b,c,P,Q,X,W,A)
    end
end

function CoefficientsPGLRK(name::Symbol, order::Int, a::Matrix{T}, b::Vector{T}, c::Vector{T}, P::Matrix{T}, X::Matrix{T}, W::Matrix{T}) where {T}
    CoefficientsPGLRK{T}(name, order, length(c), a, b, c, P, X, W)
end

function CoefficientsPGLRK(::Type{T}, s::Int) where {T}
    @assert s ≥ 2

    local ξᵢ::T

    # order
    o = 2s

    # obtain Gauss-Legendre nodes and weights
    gl_quad = GaussLegendreQuadrature(T,s)
    leg_basis = Legendre(T,s)

    a = zeros(T, s, s)
    t = zeros(T, s, s)
    P = zeros(T, s, s)
    X = zeros(T, s, s)
    W = zeros(T, s, s)

    for j in 1:s
        P[:,j] .= leg_basis[nodes(gl_quad), j-1]
    end

    Q = inv(P)

    X[1,1] = 0.5

    for i in 1:s-1
        ξᵢ = one(T) / sqrt(convert(T, 4i^2-1)) / 2
        X[i+1,i] = +ξᵢ
        X[i,i+1] = -ξᵢ
    end

    W[s, s-1] = +1
    W[s-1, s] = -1

    mul!(t, X, Q)
    mul!(a, P, t)

    CoefficientsPGLRK(Symbol("PGLRK", s), o, a, weights(gl_quad), nodes(gl_quad), P, X, W)
end

CoefficientsPGLRK(s::Int) = CoefficientsPGLRK(Float64, s)


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
mutable struct ParametersPGLRK{DT,TT,D,S,ET <: NamedTuple} <: Parameters{DT,TT}
    equs::ET
    tab::CoefficientsPGLRK{TT}
    Δt::TT

    h₀::DT

    t̄::TT
    t::TT

    q̄::Vector{DT}
    λ::DT

    function ParametersPGLRK{DT,D}(equs::ET, tab::CoefficientsPGLRK{TT}, Δt::TT) where {DT, TT, D, ET <: NamedTuple}
        new{DT, TT, D, tab.s, ET}(equs, tab, Δt, zero(TT), zero(TT), zero(DT), zeros(DT,D), zero(DT))
    end
end


mutable struct IntegratorCachePGLRK{DT,D,S} <: ODEIntegratorCache{DT,D}
    h::DT
    λ::DT

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

        new(zero(DT), zero(DT), q̃, ṽ, s̃, Q, V, Y)
    end
end

function IntegratorCache{ST}(params::ParametersPGLRK{DT,TT,D,S}; kwargs...) where {ST,DT,TT,D,S}
    IntegratorCachePGLRK{ST,D,S}(; kwargs...)
end

@inline CacheType(ST, params::ParametersPGLRK{DT,TT,D,S}) where {DT,TT,D,S} = IntegratorCachePGLRK{ST,D,S}


@doc raw"""
Projected Gauss-Legendre Runge-Kutta integrator.

    Reference: LUIGI BRUGNANO, FELICE IAVERNARO, AND DONATO TRIGIANTE.
        ENERGY- AND QUADRATIC INVARIANTS–PRESERVING INTEGRATORS BASED
        UPON GAUSS COLLOCATION FORMULAE.
        SIAM J. NUMER. ANAL. Vol. 50, No. 6, pp. 2897–2916, 2012.
"""
struct IntegratorPGLRK{DT, TT, D, S, PT <: ParametersPGLRK{DT,TT},
                                     ST <: NonlinearSolver{DT},
                                     IT <: InitialGuessODE{TT}} <: AbstractIntegratorRK{DT,TT}
    params::PT
    solver::ST
    iguess::IT
    caches::CacheDict{PT}

    function IntegratorPGLRK(params::ParametersPGLRK{DT,TT,D,S,ET}, solver::ST, iguess::IT, caches) where {DT,TT,D,S,ET,ST,IT}
        new{DT, TT, D, S, typeof(params), ST, IT}(params, solver, iguess, caches)
    end

    function IntegratorPGLRK{DT,D}(equations::NamedTuple, tableau::CoefficientsPGLRK{TT}, Δt::TT) where {D,DT,TT}
        # get number of stages
        S = tableau.s

        # create params
        params = ParametersPGLRK{DT,D}(equations, tableau, Δt)

        # create cache dict
        caches = CacheDict(params)

        # create solver
        solver  = create_nonlinear_solver(DT, D*S, params, caches)

        # create initial guess
        iguess = InitialGuessODE(get_config(:ig_interpolation), equations[:v], Δt)

        # create integrator
        IntegratorPGLRK(params, solver, iguess, caches)
    end

    function IntegratorPGLRK{DT,D}(v::Function, h::Function, tableau::CoefficientsPGLRK{TT}, Δt::TT; kwargs...) where {DT,TT,D}
        IntegratorPGLRK{DT,D}(NamedTuple{(:v,:h)}((v,h)), tableau, Δt; kwargs...)
    end

    function IntegratorPGLRK(equation::ODE{DT,TT}, tableau::CoefficientsPGLRK{TT}, Δt::TT; kwargs...) where {DT,TT}
        IntegratorPGLRK{DT, ndims(equation)}(get_function_tuple(equation), tableau, Δt; kwargs...)
    end
end


@inline Base.ndims(::IntegratorPGLRK{DT,TT,D,S}) where {DT,TT,D,S} = D
@inline nstages(::IntegratorPGLRK{DT,TT,D,S}) where {DT,TT,D,S} = S


function initialize!(int::IntegratorPGLRK, sol::AtomicSolutionODE)
    sol.t̄ = sol.t - timestep(int)

    equations(int)[:v](sol.t, sol.q, sol.v)

    initialize!(int.iguess, sol.t, sol.q, sol.v,
                            sol.t̄, sol.q̄, sol.v̄)

    int.params.h₀ = equations(int)[:h](sol.t, sol.q)
end


function update_params!(params::ParametersPGLRK, sol::AtomicSolutionODE)
    # set time for nonlinear solver and copy previous solution
    params.t̄  = sol.t
    params.t  = sol.t + params.Δt
    params.q̄ .= sol.q
end


function initial_guess!(int::IntegratorPGLRK{DT}, sol::AtomicSolutionODE{DT},
                        cache::IntegratorCachePGLRK{DT}=int.caches[DT]) where {DT}

    for i in eachstage(int)
        evaluate!(int.iguess, sol.q̄, sol.v̄, sol.q, sol.v,
                              cache.q̃, cache.ṽ,
                              tableau(int).c[i])

        for k in eachdim(int)
            int.solver.x[ndims(int)*(i-1)+k] = cache.ṽ[k]
        end
    end
end


function compute_stages!(x::Vector, cache::IntegratorCachePGLRK, params::ParametersPGLRK)
    compute_stages!(x, cache.Q, cache.V, cache.Y, cache.q̃, cache.ṽ, params)
end

function compute_stages!(x::Vector{ST}, Q::Vector{Vector{ST}}, V::Vector{Vector{ST}}, Y::Vector{Vector{ST}},
                         q::Vector{ST}, y::Vector{ST}, params::ParametersPGLRK{DT,TT,D,S}) where {ST,DT,TT,D,S}

    local tᵢ::TT

    # copy x to Y
    for i in 1:S
        for k in 1:D
            Y[i][k] = x[D*(i-1)+k]
        end
    end

    # compute Q=q̄+Δt*Y
    for i in 1:S
        for k in 1:D
            Q[i][k] = params.q̄[k] + params.Δt * Y[i][k]
        end
    end

    # compute V=v(T,Q)
    for i in 1:S
        tᵢ = params.t̄ + params.Δt * params.tab.c[i]
        params.equs[:v](tᵢ, Q[i], V[i])
    end

    # compute y=B*V
    y .= 0
    for k in 1:D
        for j in 1:S
            y[k] += params.tab.b[j] * V[j][k]
        end
    end

    # compute q=q̄+Δt*y
    q .= params.q̄ .+ params.Δt .* y
end

"Compute stages of projected Gauss-Legendre Runge-Kutta methods."
function function_stages!(x::Vector{ST}, b::Vector{ST}, params::ParametersPGLRK{DT,TT,D,S},
                          caches::CacheDict) where {ST,DT,TT,D,S}

    # get cache for internal stages
    cache = caches[ST]

    # compute stages from nonlinear solver solution x
    compute_stages!(x, cache, params)

    # compute tableau
    a = getTableauPGLRK(params.tab, params.λ)

    # compute b = [Y-AV]
    for i in 1:S
        for k in 1:D
            b[D*(i-1)+k] = cache.Y[i][k]
            for j in 1:S
                b[D*(i-1)+k] -= a[i,j] * cache.V[j][k]
            end
        end
    end
end


function bisection(f::Function, λmin::DT, λmax::DT;
                   xtol::AbstractFloat=SimpleSolvers.get_config(:nls_atol),
                   ftol::AbstractFloat=SimpleSolvers.get_config(:nls_atol),
                   maxiter::Integer=SimpleSolvers.get_config(:nls_nmax)) where {DT <: Number}
    a = λmin
    b = λmax
    fa = f(a)
    fb = f(b)

    # fa*fb ≤ 0 || error("Either no or multiple real roots in [λmin,λmax]")

    # local j = 0
    local λ = zero(DT)
    local fλ= zero(DT)

    for i in 1:maxiter
        # j += 1
        λ  = (a+b)/2
        fλ = f(λ)

        !isapprox(fλ, 0, atol=ftol) || break

        if fa*fλ > 0
            a  = λ  # Root is in the right half of [a,b].
            fa = fλ
        else
            b = λ  # Root is in the left half of [a,b].
        end

        abs(b-a) > xtol || break
    end

    # println(j, " bisection iterations, λ=", λ, ", f(λ)=", fλ, ", ftol=", ftol, ", abs(b-a)=", abs(b-a), ", xtol=", xtol)

    # i != maxiter || error("Max iteration number exceeded")

    return λ
end


function function_hamiltonian!(λ::Number, int::IntegratorPGLRK{DT,TT},
                               cache::IntegratorCachePGLRK{DT}=int.caches[DT]) where {DT,TT}
    # copy λ to integrator parameters
    int.params.λ = λ

    # call nonlinear solver
    solve!(int.solver)

    # print solver status
    print_solver_status(int.solver.status, int.solver.params)

    # check if solution contains NaNs or error bounds are violated
    check_solver_status(int.solver.status, int.solver.params)

    # compute vector fields at internal stages
    compute_stages!(int.solver.x, cache, int.params)

    # compute h
    cache.h = int.params.equs[:h](int.params.t, cache.q̃)

    # compute and return h₀-h
    return int.params.h₀ - cache.h
end

function function_hamiltonian!(λ::Vector, int::IntegratorPGLRK{DT,TT},
                               cache::IntegratorCachePGLRK{DT}=int.caches[DT]) where {DT,TT}
    [function_hamiltonian!(λ[1], int, cache)]
end


function integrate_step!(int::IntegratorPGLRK{DT,TT}, sol::AtomicSolutionODE{DT,TT},
                         cache::IntegratorCachePGLRK{DT}=int.caches[DT]) where {DT,TT}
    # update nonlinear solver parameters from atomic solution
    update_params!(int.params, sol)

    # compute initial guess
    initial_guess!(int, sol)

    # reset atomic solution
    reset!(sol, timestep(int))

    # determine parameter λ
    λmin = -0.2^nstages(int)
    λmax = +0.2^nstages(int)

    nls_atol = SimpleSolvers.get_config(:nls_atol)

    int.params.λ = bisection(λ -> function_hamiltonian!(λ, int, cache), λmin, λmax;
                    xtol=abs(λmax-λmin)*nls_atol,
                    ftol=int.params.h₀*nls_atol)
    # int.params.λ = nlsolve(λ -> function_hamiltonian!(λ, int, cache), [zero(DT)];
    #             xtol=abs(λmax-λmin)*nls_atol,
    #             ftol=int.params.h₀*nls_atol).zero[1]
    # println(int.params.λ)

    # compute final update
    update_solution!(sol.q, cache.V, tableau(int).b, timestep(int))

    # copy solution to initial guess
    update_vector_fields!(int.iguess, sol.t, sol.q, sol.v)
end
