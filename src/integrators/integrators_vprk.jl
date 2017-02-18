
"Parameters for right-hand side function of variational partitioned Runge-Kutta methods."
type NonlinearFunctionParametersVPRK{DT,TT,ΑT,FT} <: NonlinearFunctionParameters{DT}
    α::ΑT
    f::FT

    Δt::TT

    d::Int
    s::Int

    t_q::CoefficientsRK{TT}
    t_p::CoefficientsRK{TT}
    d_v::Vector{TT}

    t::TT

    q::Vector{DT}
    v::Vector{DT}
    p::Vector{DT}
    μ::Vector{DT}

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

    function NonlinearFunctionParametersVPRK(α, f, Δt, d, s, t_q, t_p, d_v)
        # create solution vectors
        q = zeros(DT,d)
        v = zeros(DT,d)
        p = zeros(DT,d)
        μ = zeros(DT,d)

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

        new(α, f, Δt, d, s, t_q, t_p, d_v, 0, q, v, p, μ, Q, V, P, F, Y, Z, tQ, tV, tP, tF)
    end
end

"Compute stages of variational partitioned Runge-Kutta methods."
function function_stages!{DT,TT,ΑT,FT}(y::Vector{DT}, b::Vector{DT}, params::NonlinearFunctionParametersVPRK{DT,TT,ΑT,FT})
    local tqᵢ::TT
    local tpᵢ::TT
    local tf::DT
    local sl::Int = div(params.s+1, 2)

    # copy y to V
    for i in 1:params.s
        for k in 1:params.d
            params.V[k,i] = y[params.d*(i-1)+k]
        end
    end

    # compute Y and Q
    for i in 1:params.s
        for k in 1:params.d
            params.Y[k,i] = 0
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

    # compute b = - [(Y-AV), (P-AF), ΛV]
    for i in 1:params.s
        for k in 1:params.d
            tf = 0
            for j in 1:params.s
                tf += params.t_p.a[i,j] * params.F[k,j]
            end
            b[params.d*(i-1)+k] = - (params.P[k,i] - params.p[k]) + params.Δt * tf
        end
    end

    if length(params.d_v) > 0
        # compute μ
        for k in 1:params.d
            tf = 0
            for j in 1:params.s
                tf += params.t_p.a[sl,j] * params.F[k,j]
            end
            params.μ[k] = params.t_p.b[sl] / params.d_v[sl] * ( - params.P[k,sl] + params.p[k] + params.Δt * tf )
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
immutable IntegratorVPRK{DT,TT,ΑT,FT,GT,VT,ST,IT} <: Integrator{DT,TT}
    equation::IODE{DT,TT,ΑT,FT,GT,VT}
    tableau::TableauVPRK{TT}
    Δt::TT

    solver::ST
    iguess::InitialGuessIODE{DT,TT,VT,FT,IT}

    q::Array{DT,1}
    v::Array{DT,1}
    p::Array{DT,1}

    qₑᵣᵣ::Vector{DT}
    pₑᵣᵣ::Vector{DT}

    y::Array{DT,1}
    z::Array{DT,1}

    Q::Array{DT,2}
    V::Array{DT,2}
    P::Array{DT,2}
    F::Array{DT,2}
end

function IntegratorVPRK{DT,TT,ΑT,FT,GT,VT}(equation::IODE{DT,TT,ΑT,FT,GT,VT}, tableau::TableauVPRK{TT}, Δt::TT;
                                        nonlinear_solver=DEFAULT_NonlinearSolver,
                                        nmax=DEFAULT_nmax, atol=DEFAULT_atol, rtol=DEFAULT_rtol, stol=DEFAULT_stol,
                                        interpolation=HermiteInterpolation{DT})
    D = equation.d
    S = tableau.s

    N = D*S

    if isdefined(tableau, :d)
        d_v = tableau.d
    else
        d_v = DT[]
    end

    # create solution vector for internal stages / nonlinear solver
    z = zeros(DT, N)

    # create params
    params = NonlinearFunctionParametersVPRK{DT,TT,ΑT,FT}(
                                                equation.α, equation.f,
                                                Δt, D, S,
                                                tableau.q, tableau.p,
                                                d_v)

    # create solver
    solver = nonlinear_solver(z, params; nmax=nmax, atol=atol, rtol=rtol, stol=stol)

    # create initial guess
    iguess = InitialGuessIODE(interpolation, equation, Δt)

    # create integrator
    qₑᵣᵣ = zeros(DT,D)
    pₑᵣᵣ = zeros(DT,D)

    y = zeros(DT,D)
    z = zeros(DT,D)

    IntegratorVPRK{DT, TT, ΑT, FT, GT, VT, typeof(solver), typeof(iguess.int)}(
                                        equation, tableau, Δt, solver, iguess,
                                        params.q, params.v, params.p,
                                        qₑᵣᵣ, pₑᵣᵣ, y, z,
                                        params.Q, params.V, params.P, params.F)
end


"Integrate ODE with variational partitioned Runge-Kutta integrator."
function integrate!{DT,TT,ΑT,FT,GT,VT,N}(int::IntegratorVPRK{DT,TT,ΑT,FT,GT,VT}, sol::SolutionPDAE{DT,TT,N})
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

            # call nonlinear solver
            solve!(int.solver)

            if !solverStatusOK(int.solver.status, int.solver.params)
                println(int.solver.status, ", it=", n)
            end

            if int.solver.status.rₐ == NaN
                break
            end

            # compute final update
            simd_mult!(int.y, int.V, int.tableau.q.b)
            simd_mult!(int.z, int.F, int.tableau.p.b)
            simd_axpy!(int.Δt, int.y, int.q, int.qₑᵣᵣ)
            simd_axpy!(int.Δt, int.z, int.p, int.pₑᵣᵣ)

            # copy to solution
            copy_solution!(sol, int.q, int.p, n, m)
        end
    end
    nothing
end
