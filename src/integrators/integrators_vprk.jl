
"Parameters for right-hand side function of variational partitioned Runge-Kutta methods."
type NonlinearFunctionParametersVPRK{DT,TT,FT,PT} <: NonlinearFunctionParameters{DT}
    f_f::FT
    f_p::PT

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

    y::Vector{DT}
    z::Vector{DT}

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

    function NonlinearFunctionParametersVPRK(f_f, f_p, Δt, d, s, t_q, t_p, d_v)
        # create solution vectors
        q = zeros(DT,d)
        v = zeros(DT,d)
        p = zeros(DT,d)
        μ = zeros(DT,d)

        y = zeros(DT,d)
        z = zeros(DT,d)

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

        new(f_f, f_p, Δt, d, s, t_q, t_p, d_v, 0, q, v, p, μ, y, z, Q, V, P, F, Y, Z, tQ, tV, tP, tF)
    end
end

"Compute stages of variational partitioned Runge-Kutta methods."
function function_stages!{DT,TT,FT,PT}(y::Vector{DT}, b::Vector{DT}, params::NonlinearFunctionParametersVPRK{DT,TT,FT,PT})
    local tqᵢ::TT
    local tpᵢ::TT

    for i in 1:params.s
        for k in 1:params.d
            # copy y to Y and Z
            params.Y[k,i] = y[2*(params.d*(i-1)+k-1)+1]
            params.Z[k,i] = y[2*(params.d*(i-1)+k-1)+2]

            # compute Q and V
            params.Q[k,i] = params.q[k] + params.Δt * params.Y[k,i]
            params.V[k,i] = params.Z[k,i]
        end

        # compute f(X)
        tqᵢ = params.t + params.Δt * params.t_q.c[i]
        tpᵢ = params.t + params.Δt * params.t_p.c[i]

        simd_copy_xy_first!(params.tQ, params.Q, i)
        simd_copy_xy_first!(params.tV, params.V, i)
        params.f_f(tpᵢ, params.tQ, params.tV, params.tF)
        params.f_p(tqᵢ, params.tQ, params.tV, params.tP)
        simd_copy_yx_first!(params.tP, params.P, i)
        simd_copy_yx_first!(params.tF, params.F, i)
    end

    if length(params.d_v) > 0
        for k in 1:params.d
            params.μ[k] = y[2*params.d*params.s+k]
        end
    end

    # compute b = - [(Y-AV), (P-AF), ΛV]
    for i in 1:params.s
        for k in 1:params.d
            b[2*(params.d*(i-1)+k-1)+1] = - params.Y[k,i]
            b[2*(params.d*(i-1)+k-1)+2] = - params.P[k,i] + params.p[k]
            for j in 1:params.s
                b[2*(params.d*(i-1)+k-1)+1] += params.t_q.a[i,j] * params.V[k,j]
                b[2*(params.d*(i-1)+k-1)+2] += params.t_p.a[i,j] * params.F[k,j] * params.Δt
            end
        end
    end

    if length(params.d_v) > 0
        for i in 1:params.s
            for k in 1:params.d
                b[2*(params.d*(i-1)+k-1)+2] -= params.d_v[i] * params.μ[k]
            end
        end

        for k in 1:params.d
            b[2*params.d*params.s+k] = 0
            for i in 1:params.s
                b[2*params.d*params.s+k] -= params.V[k,i] * params.d_v[i]
            end
        end
    end
end


"Variational partitioned Runge-Kutta integrator."
immutable IntegratorVPRK{DT, TT, FT, PT, VT, ST, IT} <: Integrator{DT, TT}
    equation::IODE{DT,TT,FT,PT,VT}
    tableau::TableauVPRK{TT}
    Δt::TT

    solver::ST
    iguess::InitialGuessIODE{DT, TT, VT, FT, IT}

    q::Array{DT,1}
    v::Array{DT,1}
    p::Array{DT,1}
    y::Array{DT,1}
    z::Array{DT,1}
    Q::Array{DT,2}
    V::Array{DT,2}
    P::Array{DT,2}
    F::Array{DT,2}
end

function IntegratorVPRK{DT,TT,FT,PT,VT}(equation::IODE{DT,TT,FT,PT,VT}, tableau::TableauVPRK{TT}, Δt::TT;
                                        nonlinear_solver=QuasiNewtonSolver,
                                        interpolation=HermiteInterpolation{DT})
    D = equation.d
    S = tableau.s

    N = 2*D*S

    if isdefined(tableau, :d)
        N += D
        d_v = tableau.d
    else
        d_v = DT[]
    end

    # create solution vector for internal stages / nonlinear solver
    z = zeros(DT, N)

    # create params
    params = NonlinearFunctionParametersVPRK{DT,TT,FT,PT}(
                                                equation.f, equation.p,
                                                Δt, D, S,
                                                tableau.q, tableau.p,
                                                d_v)

    # create solver
    solver = nonlinear_solver(z, params)

    # create initial guess
    iguess = InitialGuessIODE(interpolation, equation, Δt)

    # create integrator
    IntegratorVPRK{DT, TT, FT, PT, VT, typeof(solver), typeof(iguess.int)}(
                                        equation, tableau, Δt, solver, iguess,
                                        params.q, params.v, params.p, params.y, params.z,
                                        params.Q, params.V, params.P, params.F)
end


"Integrate ODE with variational partitioned Runge-Kutta integrator."
function integrate!{DT,TT,FT,PT,VT,N}(int::IntegratorVPRK{DT,TT,FT,PT,VT}, sol::SolutionPODE{DT,TT,N})
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
                    int.solver.x[2*(int.equation.d*(i-1)+k-1)+1] = (int.y[k] - int.q[k])/(int.Δt)
                    int.solver.x[2*(int.equation.d*(i-1)+k-1)+2] = int.v[k]
                end
            end
            if isdefined(int.tableau, :d)
                for k in 1:int.equation.d
                    int.solver.x[2*int.equation.d*int.tableau.s+k] = 0
                end
            end

            # call nonlinear solver
            solve!(int.solver)

            if !solverStatusOK(int.solver.status, int.solver.params)
                println(int.solver.status)
            end

            # compute final update
            simd_mult!(int.y, int.V, int.tableau.q.b)
            simd_mult!(int.z, int.F, int.tableau.p.b)
            simd_axpy!(int.Δt, int.y, int.q)
            simd_axpy!(int.Δt, int.z, int.p)

            # copy to solution
            copy_solution!(sol, int.q, int.p, n, m)
        end
    end
    nothing
end
