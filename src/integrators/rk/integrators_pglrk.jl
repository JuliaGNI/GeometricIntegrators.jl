
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

function getTableauPGLRK(coeff::CoefficientsPGLRK{T}, λ::T, a::Matrix{T}) where {T}
    a .= coeff.a .+ λ .* coeff.A
end


"Parameters for right-hand side function of variational partitioned Runge-Kutta methods."
mutable struct ParametersPGLRK{DT,TT,D,S,ET} <: Parameters{DT,TT}
    equ::ET
    tab::CoefficientsPGLRK{TT}
    Δt::TT

    A_q::Array{TT, 3}
    A_p::Array{TT, 3}
    a_q::Matrix{TT}
    a_p::Matrix{TT}

    t::TT

    q::Vector{DT}
    p::Vector{DT}


    function ParametersPGLRK{DT,TT,D,S,ET}(equ, tab, Δt) where {DT,TT,D,S,ET}
        # create coefficient matrices
        A_q = zeros(TT, S, S, D)
        A_p = zeros(TT, S, S, D)
        a_q = zero(tab.a)
        a_p = zero(tab.a)

        # create solution vectors
        q = zeros(DT,D)
        p = zeros(DT,D)

        new(equ, tab, Δt, A_q, A_p, a_q, a_p, 0, q, p)
    end
end

struct NonlinearFunctionCachePGLRK{DT}
    Q::Vector{Vector{DT}}
    V::Vector{Vector{DT}}
    P::Vector{Vector{DT}}
    F::Vector{Vector{DT}}
    Y::Vector{Vector{DT}}
    Z::Vector{Vector{DT}}

    v::Vector{DT}
    f::Vector{DT}
    y::Vector{DT}
    z::Vector{DT}

    q̅::Vector{DT}
    p̅::Vector{DT}
    θ::Vector{DT}
    λ::Vector{DT}


    function NonlinearFunctionCachePGLRK{DT}(D,S) where {DT}

        # create internal stage vectors
        Q = create_internal_stage_vector(DT, D, S)
        V = create_internal_stage_vector(DT, D, S)
        P = create_internal_stage_vector(DT, D, S)
        F = create_internal_stage_vector(DT, D, S)
        Y = create_internal_stage_vector(DT, D, S)
        Z = create_internal_stage_vector(DT, D, S)

        # create update vectors
        v = zeros(DT,D)
        f = zeros(DT,D)
        y = zeros(DT,D)
        z = zeros(DT,D)

        # create new solution vectors
        q̅ = zeros(DT,D)
        p̅ = zeros(DT,D)
        θ = zeros(DT,D)
        λ = zeros(DT,D)

        new(Q, V, P, F, Y, Z, v, f, y, z, q̅, p̅, θ, λ)
    end
end

function compute_stages!(x::Vector{ST}, Q::Vector{Vector{ST}}, V::Vector{Vector{ST}},
                                        P::Vector{Vector{ST}}, F::Vector{Vector{ST}},
                                        Y::Vector{Vector{ST}}, Z::Vector{Vector{ST}},
                                        q̅::Vector{ST}, p̅::Vector{ST},
                                        θ::Vector{ST}, λ::Vector{ST},
                                        y::Vector{ST}, z::Vector{ST},
                                        params::ParametersPGLRK{DT,TT,D,S}) where {ST,DT,TT,D,S}

    local tᵢ::TT
    local t₀::TT = params.t
    local t₁::TT = params.t + params.Δt

    # copy x to V
    for i in 1:S
        for k in 1:D
            V[i][k] = x[D*(i-1)+k]
        end
    end

    # copy x to λ, q̅ and p̅
    for k in 1:D
        q̅[k] = x[D*(S+0)+k]
        p̅[k] = x[D*(S+1)+k]
        λ[k] = x[D*(S+2)+k]
    end

    # compute tableaus
    for k in 1:D
        getTableauPGLRK(params.tab, λ[k], params.a_q)
        # get_symplectic_conjugate_coefficients(params.a_q, params.tab.b, params.a_p)
        for j=1:S
            for i=1:S
                params.A_q[i,j,k] = params.a_q[i,j]
                params.A_p[i,j,k] = params.a_q[i,j]
            end
        end
    end

    # compute Y and Q
    for i in 1:S
        for k in 1:D
            Y[i][k] = 0
            for j in 1:S
                Y[i][k] += params.A_q[i,j,k] * V[j][k]
            end
            Q[i][k] = params.q[k] + params.Δt * Y[i][k]
        end
    end

    # compute P=α(Q,V) and F=f(Q,V)
    for i in 1:S
        tᵢ = params.t + params.Δt * params.tab.c[i]
        params.equ.α(tᵢ, Q[i], V[i], P[i])
        params.equ.f(tᵢ, Q[i], V[i], F[i])
    end

    # compute Z
    for i in 1:S
        for k in 1:D
            Z[i][k] = 0
            for j in 1:S
                Z[i][k] += params.A_p[i,j,k] * F[j][k]
            end
        end
    end

    # compute y and z
    y .= 0
    z .= 0
    for k in 1:D
        for j in 1:S
            y[k] += params.tab.b[j] * V[j][k]
            z[k] += params.tab.b[j] * F[j][k]
        end
    end

    # compute θ=α(q̅,λ)
    params.equ.α(t₁, q̅, λ, θ)

end

"Compute stages of variational partitioned Runge-Kutta methods."
@generated function function_stages!(x::Vector{ST}, b::Vector{ST}, params::ParametersPGLRK{DT,TT,D,S}) where {ST,DT,TT,D,S}

    cache = NonlinearFunctionCachePGLRK{ST}(D, S)

    quote
        compute_stages!(x, $cache.Q, $cache.V, $cache.P, $cache.F, $cache.Y, $cache.Z,
                           $cache.q̅, $cache.p̅, $cache.θ, $cache.λ, $cache.y, $cache.z, params)

        # compute b = [P-p-AF]
        for i in 1:S
            for k in 1:D
                b[D*(i-1)+k] = $cache.P[i][k] - params.p[k] - params.Δt * $cache.Z[i][k]
            end
        end

        # compute b = [q̅-q-bV, p̅-p-bF, p̅-α(q̅)]
        for k in 1:D
            b[D*(S+0)+k] = $cache.q̅[k] - params.q[k] - params.Δt * $cache.y[k]
            b[D*(S+1)+k] = $cache.p̅[k] - params.p[k] - params.Δt * $cache.z[k]
            b[D*(S+2)+k] = $cache.p̅[k] - $cache.θ[k]
        end
    end
end

"Variational partitioned Runge-Kutta integrator."
struct IntegratorPGLRK{DT, TT, PT <: ParametersPGLRK{DT,TT},
                               ST <: NonlinearSolver{DT},
                               IT <: InitialGuessPODE{DT,TT}, N} <: Integrator{DT,TT}
    params::PT
    solver::ST
    iguess::IT
    fcache::NonlinearFunctionCachePGLRK{DT}

    q::Vector{Vector{TwicePrecision{DT}}}
    p::Vector{Vector{TwicePrecision{DT}}}
end

function IntegratorPGLRK(equation::IODE{DT,TT,ΑT,FT,GT,VT,N}, tableau::CoefficientsPGLRK{TT}, Δt::TT) where {DT,TT,ΑT,FT,GT,VT,N}
    D = equation.d
    M = equation.n
    S = tableau.s

    # create params
    params = ParametersPGLRK{DT,TT,D,S,typeof(equation)}(equation, tableau, Δt)

    # create solver
    solver = create_nonlinear_solver(DT, D*(S+3), params)

    # create initial guess
    iguess = InitialGuessPODE(get_config(:ig_interpolation), equation, Δt)

    # create cache for internal stage vectors and update vectors
    fcache = NonlinearFunctionCachePGLRK{DT}(D, S)

    # create solution vectors
    q = create_solution_vector(DT, D, M)
    p = create_solution_vector(DT, D, M)

    IntegratorPGLRK{DT, TT, typeof(params), typeof(solver), typeof(iguess), N}(
                params, solver, iguess, fcache, q, p)
end

equation(integrator::IntegratorPGLRK) = integrator.params.equ
timestep(integrator::IntegratorPGLRK) = integrator.params.Δt


function initialize!(int::IntegratorPGLRK{DT,TT}, sol::SolutionPDAE, m::Int) where {DT,TT}
    @assert m ≥ 1
    @assert m ≤ sol.ni

    # copy initial conditions from solution
    get_initial_conditions!(sol, int.q[m], int.p[m], m)

    # initialise initial guess
    initialize!(int.iguess, m, sol.t[0], int.q[m], int.p[m])
end

function initial_guess!(int::IntegratorPGLRK{DT,TT}, m::Int) where {DT,TT}
    for i in 1:int.params.tab.s
        evaluate!(int.iguess, m, int.fcache.y, int.fcache.z, int.fcache.v, int.fcache.f, int.params.tab.c[i], int.params.tab.c[i])
        for k in 1:int.params.equ.d
            int.solver.x[int.params.equ.d*(i-1)+k] = int.fcache.v[k]
        end
    end
    evaluate!(int.iguess, m, int.fcache.y, int.fcache.z, int.fcache.v, int.fcache.f, one(TT), one(TT))
    for k in 1:int.params.equ.d
        int.solver.x[int.params.equ.d*(int.params.tab.s+0)+k] = int.fcache.y[k]
    end
    for k in 1:int.params.equ.d
        int.solver.x[int.params.equ.d*(int.params.tab.s+1)+k] = int.fcache.z[k]
    end
    for k in 1:int.params.equ.d
        int.solver.x[int.params.equ.d*(int.params.tab.s+2)+k] = 0
    end
end


"Integrate PODE with variational partitioned Runge-Kutta integrator."
function integrate_step!(int::IntegratorPGLRK{DT,TT}, sol::SolutionPDAE{DT,TT}, m::Int, n::Int) where {DT,TT}
    # set time for nonlinear solver
    int.params.t  = sol.t[0] + (n-1)*int.params.Δt

    # copy previous solution
    int.params.q .= int.q[m]
    int.params.p .= int.p[m]

    # compute initial guess
    initial_guess!(int, m)

    # call nonlinear solver
    solve!(int.solver)

    # print solver status
    print_solver_status(int.solver.status, int.solver.params, n)

    # check if solution contains NaNs or error bounds are violated
    check_solver_status(int.solver.status, int.solver.params, n)

    # compute vector fields at internal stages
    compute_stages!(int.solver.x, int.fcache.Q, int.fcache.V,
                                  int.fcache.P, int.fcache.F,
                                  int.fcache.Y, int.fcache.Z,
                                  int.fcache.q̅, int.fcache.p̅,
                                  int.fcache.θ, int.fcache.λ,
                                  int.fcache.y, int.fcache.z,
                                  int.params)

    # compute final update
    update_solution!(int.q[m], int.fcache.V, int.params.tab.b, int.params.Δt)
    update_solution!(int.p[m], int.fcache.F, int.params.tab.b, int.params.Δt)

    # take care of periodic solutions
    cut_periodic_solution!(int.q[m], int.params.equ.periodicity)

    # copy to solution
    copy_solution!(sol, int.q[m], int.p[m], n, m)
end
