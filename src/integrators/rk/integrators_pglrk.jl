
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
        A = zeros(a)
        B = zeros(a)

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


"Parameters for right-hand side function of variational partitioned Runge-Kutta methods."
mutable struct ParametersPGLRK{DT,TT,ΑT,FT,GT} <: Parameters{DT,TT}
    α::ΑT
    f::FT
    g::GT

    Δt::TT

    d::Int
    s::Int

    tab::CoefficientsPGLRK{TT}
    A_q::Array{TT, 3}
    A_p::Array{TT, 3}
    a_q::Matrix{TT}
    a_p::Matrix{TT}

    t::TT

    q::Vector{DT}
    v::Vector{DT}
    p::Vector{DT}

    q̅::Vector{DT}
    p̅::Vector{DT}
    λ::Vector{DT}

    Q::Matrix{DT}
    V::Matrix{DT}
    P::Matrix{DT}
    F::Matrix{DT}

    Y::Matrix{DT}
    Z::Matrix{DT}

    tQ::Vector{DT}
    tV::Vector{DT}
    tP::Vector{DT}
    tF::Vector{DT}

    function ParametersPGLRK{DT,TT,ΑT,FT,GT}(α, f, g, Δt, d, s, tab) where {DT,TT,ΑT,FT,GT}
        # create coefficient matrices
        A_q = zeros(TT, s, s, d)
        A_p = zeros(TT, s, s, d)
        a_q = zeros(tab.a)
        a_p = zeros(tab.a)

        # create solution vectors
        q = zeros(DT,d)
        v = zeros(DT,d)
        p = zeros(DT,d)

        q̅ = zeros(DT,d)
        p̅ = zeros(DT,d)
        λ = zeros(DT,d)

        # create internal stage vectors
        Q = zeros(DT,d,s)
        V = zeros(DT,d,s)
        P = zeros(DT,d,s)
        F = zeros(DT,d,s)

        Y = zeros(DT,d,s)
        Z = zeros(DT,d,s)

        # create temporary vectors
        tQ = zeros(DT,d)
        tV = zeros(DT,d)
        tP = zeros(DT,d)
        tF = zeros(DT,d)

        new(α, f, g, Δt, d, s, tab, A_q, A_p, a_q, a_p, 0, q, v, p, q̅, p̅, λ, Q, V, P, F, Y, Z, tQ, tV, tP, tF)
    end
end

"Compute stages of variational partitioned Runge-Kutta methods."
function function_stages!(y::Vector{DT}, b::Vector{DT}, params::ParametersPGLRK{DT,TT,ΑT,FT,GT}) where {DT,TT,ΑT,FT,GT}
    local tᵢ::TT
    local t₀::TT = params.t
    local t₁::TT = params.t + params.Δt
    local tf::DT

    # copy y to V
    for i in 1:params.s
        for k in 1:params.d
            params.V[k,i] = y[params.d*(i-1)+k]
        end
    end

    # copy y to λ, q̅ and p̅
    for k in 1:params.d
        params.q̅[k] = y[params.d*(params.s+0)+k]
        params.p̅[k] = y[params.d*(params.s+1)+k]
        params.λ[k] = y[params.d*(params.s+2)+k]
    end

    # compute tableaus
    for k in 1:params.d
        getTableauPGLRK(params.tab, params.λ[k], params.a_q)
        # get_symplectic_conjugate_coefficients(params.a_q, params.tab.b, params.a_p)
        @inbounds for j=1:params.s
            for i=1:params.s
                params.A_q[i,j,k] = params.a_q[i,j]
                params.A_p[i,j,k] = params.a_q[i,j]
            end
        end
    end

    # compute Y and Q
    for i in 1:params.s
        for k in 1:params.d
            params.Y[k,i] = 0
            for j in 1:params.s
                params.Y[k,i] += params.A_q[i,j,k] * params.V[k,j]
            end
            params.Q[k,i] = params.q[k] + params.Δt * params.Y[k,i]
        end
    end

    # compute P=α(Q) and F=f(Q)
    for i in 1:params.s
        tᵢ = params.t + params.Δt * params.tab.c[i]

        simd_copy_xy_first!(params.tQ, params.Q, i)
        simd_copy_xy_first!(params.tV, params.V, i)
        params.α(tᵢ, params.tQ, params.tV, params.tP)
        params.f(tᵢ, params.tQ, params.tV, params.tF)
        simd_copy_yx_first!(params.tP, params.P, i)
        simd_copy_yx_first!(params.tF, params.F, i)
    end

    # compute tP=α(tQ)
    params.α(t₁, params.q̅, params.λ, params.tP)

    # compute b = - [P-AF-U]
    for i in 1:params.s
        for k in 1:params.d
            params.Z[k,i] = 0
            for j in 1:params.s
                params.Z[k,i] += params.A_p[i,j,k] * params.F[k,j]
            end
            b[params.d*(i-1)+k] = - (params.P[k,i] - params.p[k]) + params.Δt * params.Z[k,i]
        end
    end

    # compute b = - [q-bV-U]
    for k in 1:params.d
        params.tV[k] = 0
        for j in 1:params.s
            params.tV[k] += params.tab.b[j] * params.V[k,j]
        end
        b[params.d*(params.s+0)+k] = - (params.q̅[k] - params.q[k]) + params.Δt * params.tV[k]
    end

    # compute b = - [p-bF-G]
    for k in 1:params.d
        params.tF[k] = 0
        for j in 1:params.s
            params.tF[k] += params.tab.b[j] * params.F[k,j]
        end
        b[params.d*(params.s+1)+k] = - (params.p̅[k] - params.p[k]) + params.Δt * params.tF[k]
    end

    # compute b = - [p-α(q)]
    for k in 1:params.d
        b[params.d*(params.s+2)+k] = - params.p̅[k] + params.tP[k]
    end
end


"Variational partitioned Runge-Kutta integrator."
struct IntegratorPGLRK{DT,TT,ΑT,FT,GT,VT,ST,IT} <: Integrator{DT,TT}
    equation::IODE{DT,TT,ΑT,FT,GT,VT}
    tableau::CoefficientsPGLRK{TT}
    Δt::TT

    solver::ST
    iguess::InitialGuessPODE{DT,TT,VT,FT,IT}

    q::Array{DT,1}
    v::Array{DT,1}
    p::Array{DT,1}

    q̅::Array{DT,1}
    p̅::Array{DT,1}
    λ::Array{DT,1}

    qₑᵣᵣ::Vector{DT}
    pₑᵣᵣ::Vector{DT}

    Q::Array{DT,2}
    V::Array{DT,2}

    P::Array{DT,2}
    F::Array{DT,2}

    y::Array{DT,1}
    z::Array{DT,1}
end

function IntegratorPGLRK(equation::IODE{DT,TT,ΑT,FT,GT,VT}, tableau::CoefficientsPGLRK{TT}, Δt::TT) where {DT,TT,ΑT,FT,GT,VT}
    D = equation.d
    S = tableau.s

    N = D*(S+3)

    # create update vectors
    y = zeros(DT,D)
    z = zeros(DT,D)

    # create params
    params = ParametersPGLRK{DT,TT,ΑT,FT,GT}(
                                                equation.α, equation.f, equation.g,
                                                Δt, D, S,
                                                tableau)

    # create solver
    solver = nonlinear_solver(zeros(DT,N), params)

    # create initial guess
    iguess = InitialGuessPODE(get_config(:ig_interpolation), equation, Δt)

    # create integrator
    qₑᵣᵣ = zeros(DT,D)
    pₑᵣᵣ = zeros(DT,D)

    IntegratorPGLRK{DT, TT, ΑT, FT, GT, VT, typeof(solver), typeof(iguess.int)}(
                                        equation, tableau, Δt, solver, iguess,
                                        params.q, params.v, params.p,
                                        params.q̅, params.p̅, params.λ,
                                        qₑᵣᵣ, pₑᵣᵣ,
                                        params.Q, params.V, params.P, params.F,
                                        y, z)
end


"Integrate ODE with variational partitioned Runge-Kutta integrator."
function integrate!(int::IntegratorPGLRK{DT,TT,ΑT,FT,GT,VT}, sol::SolutionPDAE{DT,TT,N}) where {DT,TT,ΑT,FT,GT,VT,N}
    # loop over initial conditions
    for m in 1:sol.ni
        local j::Int
        local tqᵢ::TT
        local tpᵢ::TT

        # copy initial conditions from solution
        get_initial_conditions!(sol, int.q, int.p, m)

        # initialise initial guess
        initialize!(int.iguess, m, sol.t[0], int.q, int.p)

        for n in 1:sol.ntime
            # set time for nonlinear solver
            int.solver.Fparams.t = sol.t[n]

            # copy previous solution to initial guess
            update!(int.iguess, sol.t[n], int.q, int.p)

            # compute initial guess
            for i in 1:int.tableau.s
                evaluate!(int.iguess, int.y, int.z, int.v, int.tableau.c[i], int.tableau.c[i])
                for k in 1:int.equation.d
                    int.solver.x[int.equation.d*(i-1)+k] = int.v[k]
                end
            end
            evaluate!(int.iguess, int.y, int.z, int.v, one(TT), one(TT))
            for k in 1:int.equation.d
                int.solver.x[int.equation.d*(int.tableau.s+0)+k] = int.y[k]
            end
            for k in 1:int.equation.d
                int.solver.x[int.equation.d*(int.tableau.s+1)+k] = int.z[k]
            end
            for k in 1:int.equation.d
                int.solver.x[int.equation.d*(int.tableau.s+2)+k] = 0
            end

            # call nonlinear solver
            solve!(int.solver)

            if !check_solver_status(int.solver.status, int.solver.params)
                println(int.solver.status, ", it=", n)
            end

            if int.solver.status.rₐ == NaN
                break
            end

            # copy solution
            simd_mult!(int.y, int.V, int.tableau.b)
            simd_mult!(int.z, int.F, int.tableau.b)
            simd_axpy!(int.Δt, int.y, int.q, int.qₑᵣᵣ)
            simd_axpy!(int.Δt, int.z, int.p, int.pₑᵣᵣ)

            # take care of periodic solutions
            for k in 1:int.equation.d
                if int.equation.periodicity[k] ≠ 0
                    int.q[k] = mod(int.q[k], int.equation.periodicity[k])
                end
            end

            # copy to solution
            copy_solution!(sol, int.q, int.p, int.λ, n, m)
        end
    end
    nothing
end
