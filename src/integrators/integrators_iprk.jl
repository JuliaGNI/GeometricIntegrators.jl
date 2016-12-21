
"Parameters for right-hand side function of implicit partitioned Runge-Kutta methods."
type NonlinearFunctionParametersIPRK{DT,TT,FT,GT} <: NonlinearFunctionParameters{DT}
    f::FT
    g::GT
    Δt::TT

    d::Int
    s::Int
    a_q::Matrix{TT}
    a_p::Matrix{TT}
    c_q::Vector{TT}
    c_p::Vector{TT}

    t::TT
    q::Vector{DT}
    p::Vector{DT}
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

    function NonlinearFunctionParametersIPRK(f, g, Δt, d, s, a_q, a_p, c_q, c_p)
        # create solution vectors
        q = zeros(DT,d)
        p = zeros(DT,d)
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

        new(f, g, Δt, d, s, a_q, a_p, c_q, c_p, 0, q, p, y, z, Q, V, P, F, Y, Z, tQ, tV, tP, tF)
    end
end

"Compute stages of implicit partitioned Runge-Kutta methods."
function function_stages!{DT,TT,FT,GT}(y::Vector{DT}, b::Vector{DT}, params::NonlinearFunctionParametersIPRK{DT,TT,FT,GT})
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
        tqᵢ = params.t + params.Δt * params.c_q[i]
        tpᵢ = params.t + params.Δt * params.c_p[i]

        simd_copy_xy_first!(params.tQ, params.Q, i)
        simd_copy_xy_first!(params.tV, params.V, i)
        params.f(tqᵢ, params.tQ, params.tV, params.tP)
        params.g(tpᵢ, params.tQ, params.tV, params.tF)
        simd_copy_yx_first!(params.tP, params.P, i)
        simd_copy_yx_first!(params.tF, params.F, i)
    end

    # compute b = - [(Y-AV), (P-AF)]
    for i in 1:params.s
        for k in 1:params.d
            b[2*(params.d*(i-1)+k-1)+1] = - params.Y[k,i]
            b[2*(params.d*(i-1)+k-1)+2] = - params.P[k,i] + params.p[k]
            for j in 1:params.s
                b[2*(params.d*(i-1)+k-1)+1] += params.a_q[i,j] * params.V[k,j]
                b[2*(params.d*(i-1)+k-1)+2] += params.a_p[i,j] * params.F[k,j] * params.Δt
            end
        end
    end
end


"Implicit partitioned Runge-Kutta integrator."
immutable IntegratorIPRK{DT, TT, FT, GT, ST} <: Integrator{DT, TT}
    equation::IODE{DT,TT,FT,GT}
    tableau::TableauIPRK{TT}
    Δt::TT

    solver::ST

    q::Array{DT,1}
    p::Array{DT,1}
    y::Array{DT,1}
    z::Array{DT,1}
    Q::Array{DT,2}
    V::Array{DT,2}
    P::Array{DT,2}
    F::Array{DT,2}
    Y::Array{DT,2}
    Z::Array{DT,2}
end

function IntegratorIPRK{DT,TT,FT,GT}(equation::IODE{DT,TT,FT,GT}, tableau::TableauIPRK{TT}, Δt::TT)
    D = equation.d
    S = tableau.s

    # create solution vector for internal stages / nonlinear solver
    z = zeros(DT, 2*D*S)

    # create params
    params = NonlinearFunctionParametersIPRK{DT,TT,FT,GT}(
                                                equation.f, equation.p,
                                                Δt, D, S,
                                                tableau.a_q, tableau.a_p,
                                                tableau.c_q, tableau.c_p)

    # create solver
    solver = QuasiNewtonSolver(z, params)
    # TODO allow for other nonlinear solvers based on constructor argument

    # create integrator
    IntegratorIPRK{DT, TT, FT, GT, typeof(solver)}(
                                        equation, tableau, Δt, solver,
                                        params.q, params.p, params.y, params.z,
                                        params.Q, params.V, params.P, params.F,
                                        params.Y, params.Z)
end


"Integrate ODE with implicit partitioned Runge-Kutta integrator."
function integrate!{DT,TT,VT,FT,N}(int::IntegratorIPRK{DT,TT,VT,FT}, sol::SolutionPODE{DT,TT,N})
    # loop over initial conditions
    for m in 1:sol.ni
        local j::Int
        local tqᵢ::TT
        local tpᵢ::TT

        # copy initial conditions from solution
        get_initial_conditions!(sol, int.q, int.p, m)

        for n in 1:sol.ntime
            # set time for nonlinear solver
            int.solver.Fparams.t = sol.t[n]

            # compute initial guess
            # TODO
            for i in 1:int.tableau.s
                for k in 1:int.equation.d
                    int.solver.x[2*(sol.nd*(i-1)+k-1)+1] = int.q[k]
                    # TODO initial guess for velocity
                    # int.solver.x[2*(params.d*(i-1)+k-1)+1] = 0.
                end
            end

            # call nonlinear solver
            solve!(int.solver)

            if !solverStatusOK(int.solver.status, int.solver.params)
                println(int.solver.status)
            end

            # compute final update
            simd_mult!(int.y, int.V, int.tableau.b_q)
            simd_mult!(int.z, int.F, int.tableau.b_p)
            simd_axpy!(int.Δt, int.y, int.q)
            simd_axpy!(int.Δt, int.z, int.p)

            # copy to solution
            copy_solution!(sol, int.q, int.p, n, m)
        end
    end
    nothing
end
