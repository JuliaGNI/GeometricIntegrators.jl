
"Parameters for right-hand side function of variational partitioned Runge-Kutta methods."
type NonlinearFunctionParametersVPRKpSymmetric{DT,TT,ΑT,FT,GT} <: NonlinearFunctionParameters{DT}
    α::ΑT
    f::FT
    g::GT

    Δt::TT

    d::Int
    s::Int

    t_q::CoefficientsRK{TT}
    t_p::CoefficientsRK{TT}
    R∞::Int
    d_v::Vector{TT}
    b_λ::Vector{TT}

    t::TT

    q::Vector{DT}
    v::Vector{DT}
    p::Vector{DT}

    q̅::Vector{DT}
    p̅::Vector{DT}
    λ::Vector{DT}
    μ::Vector{DT}

    Q::Matrix{DT}
    V::Matrix{DT}
    P::Matrix{DT}
    F::Matrix{DT}

    Y::Matrix{DT}
    Z::Matrix{DT}
    U::Matrix{DT}
    G::Matrix{DT}

    tQ::Vector{DT}
    tV::Vector{DT}
    tP::Vector{DT}
    tF::Vector{DT}
    tU::Vector{DT}
    tG::Vector{DT}

    function NonlinearFunctionParametersVPRKpSymmetric(α, f, g, Δt, d, s, t_q, t_p, R∞, d_v)
        # create solution vectors
        q = zeros(DT,d)
        v = zeros(DT,d)
        p = zeros(DT,d)

        q̅ = zeros(DT,d)
        p̅ = zeros(DT,d)
        λ = zeros(DT,d)
        μ = zeros(DT,d)

        # create internal stage vectors
        Q = zeros(DT,d,s)
        V = zeros(DT,d,s)
        P = zeros(DT,d,s)
        F = zeros(DT,d,s)

        Y = zeros(DT,d,s)
        Z = zeros(DT,d,s)
        U = zeros(DT,d,2)
        G = zeros(DT,d,2)

        # create temporary vectors
        tQ = zeros(DT,d)
        tV = zeros(DT,d)
        tP = zeros(DT,d)
        tF = zeros(DT,d)
        tU = zeros(DT,d)
        tG = zeros(DT,d)

        b_λ = float([1, R∞])

        new(α, f, g, Δt, d, s, t_q, t_p, R∞, d_v, b_λ, 0, q, v, p, q̅, p̅, λ, μ, Q, V, P, F, Y, Z, U, G, tQ, tV, tP, tF, tU, tG)
    end
end

"Compute stages of variational partitioned Runge-Kutta methods."
function function_stages!{DT,TT,ΑT,FT,GT}(y::Vector{DT}, b::Vector{DT}, params::NonlinearFunctionParametersVPRKpSymmetric{DT,TT,ΑT,FT,GT})
    local tᵢ::TT
    local t₀::TT = params.t
    local t₁::TT = params.t + params.Δt
    local tf::DT
    local sl::Int = div(params.s+1, 2)

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

    # compute U=λ
    simd_copy_yx_first!(params.λ, params.U, 1)
    simd_copy_yx_first!(params.λ, params.U, 2)

    # compute Y and Q
    for i in 1:params.s
        for k in 1:params.d
            params.Y[k,i] = params.U[k,1]
            for j in 1:params.s
                params.Y[k,i] += params.t_q.a[i,j] * params.V[k,j]
            end
            params.Q[k,i] = params.q[k] + params.Δt * params.Y[k,i]
        end
    end

    # compute P=α(Q) and F=f(Q)
    for i in 1:params.s
        tqᵢ = params.t + params.Δt * params.t_q.c[i]
        tpᵢ = params.t + params.Δt * params.t_p.c[i]

        simd_copy_xy_first!(params.tQ, params.Q, i)
        simd_copy_xy_first!(params.tV, params.V, i)
        params.α(tqᵢ, params.tQ, params.tV, params.tP)
        params.f(tqᵢ, params.tQ, params.tV, params.tF)
        simd_copy_yx_first!(params.tP, params.P, i)
        simd_copy_yx_first!(params.tF, params.F, i)
    end

    # compute G=g(q,λ)
    params.g(t₀, params.q, params.λ, params.tG)
    simd_copy_yx_first!(params.tG, params.G, 1)

    params.g(t₁, params.q̅, params.λ, params.tG)
    simd_copy_yx_first!(params.tG, params.G, 2)

    # compute tP=α(tQ)
    params.α(t₁, params.q̅, params.λ, params.tP)

    # compute b = - [P-AF-U]
    for i in 1:params.s
        for k in 1:params.d
            params.Z[k,i] = params.G[k,1]
            for j in 1:params.s
                params.Z[k,i] += params.t_p.a[i,j] * params.F[k,j]
            end
            b[params.d*(i-1)+k] = - params.P[k,i] + params.p[k] + params.Δt * params.Z[k,i]
        end
    end

    # compute b = - [q-bV-U]
    for k in 1:params.d
        params.tV[k] = params.U[k,1] + params.R∞ * params.U[k,2]
        for j in 1:params.s
            params.tV[k] += params.t_q.b[j] * params.V[k,j]
        end
        b[params.d*(params.s+0)+k] = - (params.q̅[k] - params.q[k]) / params.Δt + params.tV[k]
    end

    # compute b = - [p-α(q)]
    for k in 1:params.d
        b[params.d*(params.s+1)+k] = - params.p̅[k] + params.tP[k]
    end

    # compute b = - [α(q)-bF-G]
    for k in 1:params.d
        params.tF[k] = params.G[k,1] + params.R∞ * params.G[k,2]
        for j in 1:params.s
            params.tF[k] += params.t_p.b[j] * params.F[k,j]
        end
        b[params.d*(params.s+2)+k] = - params.p̅[k] + params.p[k] + params.Δt * params.tF[k]
    end

    if length(params.d_v) > 0
        # compute μ
        for k in 1:params.d
            params.μ[k] = params.t_p.b[sl] / params.d_v[sl] * ( - params.P[k,sl] + params.p[k] + params.Δt * params.Z[k,sl] )
        end

        # replace equation for Pₗ with constraint on V
        for k in 1:params.d
            b[params.d*(sl-1)+k] = 0
            for i in 1:params.s
                b[params.d*(sl-1)+k] += params.V[k,i] * params.d_v[i]
            end
        end

        # modify P₁, ..., Pₛ except for Pₗ
        for i in 1:params.s
            if i ≠ sl
                tf = params.d_v[i] / params.t_p.b[i]
                for k in 1:params.d
                    b[params.d*(i-1)+k] -= tf * params.μ[k]
                end
            end
        end
    end
end


"Variational partitioned Runge-Kutta integrator."
immutable IntegratorVPRKpSymmetric{DT,TT,ΑT,FT,GT,VT,ST,IT} <: Integrator{DT,TT}
    equation::IODE{DT,TT,ΑT,FT,GT,VT}
    tableau::TableauVPRK{TT}
    Δt::TT
    b_λ::Vector{TT}

    solver::ST
    iguess::InitialGuessIODE{DT,TT,VT,FT,IT}

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
    U::Array{DT,2}

    P::Array{DT,2}
    F::Array{DT,2}
    G::Array{DT,2}

    y::Array{DT,1}
    z::Array{DT,1}
    u::Array{DT,1}
    g::Array{DT,1}
end

function IntegratorVPRKpSymmetric{DT,TT,ΑT,FT,GT,VT}(equation::IODE{DT,TT,ΑT,FT,GT,VT}, tableau::TableauVPRK{TT}, Δt::TT;
                                        nonlinear_solver=DEFAULT_NonlinearSolver,
                                        nmax=DEFAULT_nmax, atol=DEFAULT_atol, rtol=DEFAULT_rtol, stol=DEFAULT_stol,
                                        interpolation=HermiteInterpolation{DT})
    D = equation.d
    S = tableau.s

    N = D*(S+3)

    if isdefined(tableau, :d)
        d_v = tableau.d
    else
        d_v = DT[]
    end

    # create update vectors
    y = zeros(DT,D)
    z = zeros(DT,D)
    u = zeros(DT,D)
    g = zeros(DT,D)

    # create params
    params = NonlinearFunctionParametersVPRKpSymmetric{DT,TT,ΑT,FT,GT}(
                                                equation.α, equation.f, equation.g,
                                                Δt, D, S,
                                                tableau.q, tableau.p, tableau.R∞,
                                                d_v)

    # create solver
    solver = nonlinear_solver(zeros(DT,N), params; nmax=nmax, atol=atol, rtol=rtol, stol=stol)

    # create initial guess
    iguess = InitialGuessIODE(interpolation, equation, Δt)

    # create integrator
    qₑᵣᵣ = zeros(DT,D)
    pₑᵣᵣ = zeros(DT,D)

    IntegratorVPRKpSymmetric{DT, TT, ΑT, FT, GT, VT, typeof(solver), typeof(iguess.int)}(
                                        equation, tableau, Δt, params.b_λ, solver, iguess,
                                        params.q, params.v, params.p,
                                        params.q̅, params.p̅, params.λ,
                                        qₑᵣᵣ, pₑᵣᵣ,
                                        params.Q, params.V, params.U,
                                        params.P, params.F, params.G,
                                        y, z, u, g)
end


"Integrate ODE with variational partitioned Runge-Kutta integrator."
function integrate!{DT,TT,ΑT,FT,GT,VT,N}(int::IntegratorVPRKpSymmetric{DT,TT,ΑT,FT,GT,VT}, sol::SolutionPDAE{DT,TT,N})
    # loop over initial conditions
    for m in 1:sol.ni
        local j::Int
        local tqᵢ::TT
        local tpᵢ::TT

        # copy initial conditions from solution
        get_initial_conditions!(sol, int.q, int.p, m)

        # initialise initial guess
        initialize!(int.iguess, sol.t[0], int.q, int.p)

        for n in 1:sol.ntime
            # set time for nonlinear solver
            int.solver.Fparams.t = sol.t[n]

            # copy previous solution to initial guess
            update!(int.iguess, sol.t[n], int.q, int.p)

            # compute initial guess
            for i in 1:int.tableau.s
                evaluate!(int.iguess, int.y, int.z, int.v, int.tableau.q.c[i], int.tableau.p.c[i])
                for k in 1:int.equation.d
                    int.solver.x[int.equation.d*(i-1)+k] = int.v[k]
                end
            end
            # evaluate!(int.iguess, int.y, int.z, int.v, one(TT), one(TT))
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

            if !solverStatusOK(int.solver.status, int.solver.params)
                println(int.solver.status, ", it=", n)
            end

            if int.solver.status.rₐ == NaN
                break
            end

            # copy solution
            simd_mult!(int.y, int.V, int.tableau.q.b)
            simd_mult!(int.z, int.F, int.tableau.p.b)
            simd_axpy!(int.Δt, int.y, int.q, int.qₑᵣᵣ)
            simd_axpy!(int.Δt, int.z, int.p, int.pₑᵣᵣ)

            simd_mult!(int.u, int.U, int.b_λ)
            simd_mult!(int.g, int.G, int.b_λ)
            simd_axpy!(int.Δt, int.u, int.q, int.qₑᵣᵣ)
            simd_axpy!(int.Δt, int.g, int.p, int.pₑᵣᵣ)

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
