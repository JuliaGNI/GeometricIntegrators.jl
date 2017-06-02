
"Parameters for right-hand side function of implicit partitioned Runge-Kutta methods."
type NonlinearFunctionParametersIPRK{DT,TT,VT,FT} <: NonlinearFunctionParameters{DT}
    f_v::VT
    f_f::FT

    Δt::TT

    d::Int
    s::Int

    t_q::CoefficientsRK{TT}
    t_p::CoefficientsRK{TT}

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

    function NonlinearFunctionParametersIPRK(f_v, f_f, Δt, d, s, t_q, t_p)
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

        new(f_v, f_f, Δt, d, s, t_q, t_p, 0, q, p, y, z, Q, V, P, F, Y, Z, tQ, tV, tP, tF)
    end
end

"Compute stages of implicit partitioned Runge-Kutta methods."
function function_stages!{DT,TT,VT,FT}(y::Vector{DT}, b::Vector{DT}, params::NonlinearFunctionParametersIPRK{DT,TT,VT,FT})
    local tqᵢ::TT
    local tpᵢ::TT

    for i in 1:params.s
        for k in 1:params.d
            # copy y to Y and Z
            params.Y[k,i] = y[2*(params.d*(i-1)+k-1)+1]
            params.Z[k,i] = y[2*(params.d*(i-1)+k-1)+2]

            # compute Q and P
            params.Q[k,i] = params.q[k] + params.Δt * params.Y[k,i]
            params.P[k,i] = params.p[k] + params.Δt * params.Z[k,i]
        end

        # compute f(X)
        tqᵢ = params.t + params.Δt * params.t_q.c[i]
        tpᵢ = params.t + params.Δt * params.t_p.c[i]

        simd_copy_xy_first!(params.tQ, params.Q, i)
        simd_copy_xy_first!(params.tP, params.P, i)
        params.f_v(tqᵢ, params.tQ, params.tP, params.tV)
        params.f_f(tpᵢ, params.tQ, params.tP, params.tF)
        simd_copy_yx_first!(params.tV, params.V, i)
        simd_copy_yx_first!(params.tF, params.F, i)
    end

    # compute b = - [(Y-AV), (Z-AF)]
    for i in 1:params.s
        for k in 1:params.d
            b[2*(params.d*(i-1)+k-1)+1] = - params.Y[k,i]
            b[2*(params.d*(i-1)+k-1)+2] = - params.Z[k,i]
            for j in 1:params.s
                b[2*(params.d*(i-1)+k-1)+1] += params.t_q.a[i,j] * params.V[k,j]
                b[2*(params.d*(i-1)+k-1)+2] += params.t_p.a[i,j] * params.F[k,j]
            end
        end
    end
end


"Implicit partitioned Runge-Kutta integrator."
immutable IntegratorIPRK{DT, TT, VT, FT, SPT, ST} <: Integrator{DT, TT}
    equation::PODE{DT,TT,VT,FT}
    tableau::TableauIPRK{TT}
    Δt::TT

    params::SPT
    solver::ST

    q::Array{DT,1}
    p::Array{DT,1}
    y::Array{DT,1}
    z::Array{DT,1}
    Q::Array{DT,2}
    V::Array{DT,2}
    P::Array{DT,2}
    F::Array{DT,2}
end

function IntegratorIPRK{DT,TT,VT,FT}(equation::PODE{DT,TT,VT,FT}, tableau::TableauIPRK{TT}, Δt::TT;
                                     nonlinear_solver=DEFAULT_NonlinearSolver,
                                     nmax=DEFAULT_nmax, atol=DEFAULT_atol, rtol=DEFAULT_rtol, stol=DEFAULT_stol)
    D = equation.d
    S = tableau.s

    # create solution vector for internal stages / nonlinear solver
    x = zeros(DT, 2*D*S)

    # create params
    params = NonlinearFunctionParametersIPRK{DT,TT,VT,FT}(
                                                equation.v, equation.f,
                                                Δt, D, S,
                                                tableau.q, tableau.p)

    # create rhs function for nonlinear solver
    function_stages = (x,b) -> function_stages!(x, b, params)

    # create solver
    solver = nonlinear_solver(x, function_stages; nmax=nmax, atol=atol, rtol=rtol, stol=stol, autodiff=false)

    # create integrator
    IntegratorIPRK{DT, TT, VT, FT, typeof(params), typeof(solver)}(
                                        equation, tableau, Δt, params, solver,
                                        params.q, params.p, params.y, params.z,
                                        params.Q, params.V, params.P, params.F)
end


function initialize!(int::IntegratorIPRK, sol::SolutionPODE, m::Int)
    @assert m ≥ 1
    @assert m ≤ sol.ni

    # copy initial conditions from solution
    get_initial_conditions!(sol, int.q, int.p, m)
end

"Integrate ODE with implicit partitioned Runge-Kutta integrator."
function integrate_step!{DT,TT,VT,FT,N}(int::IntegratorIPRK{DT,TT,VT,FT}, sol::SolutionPODE{DT,TT,N}, m::Int, n::Int)
    # set time for nonlinear solver
    int.params.t = sol.t[n]

    # compute initial guess
    for i in 1:int.tableau.s
        for k in 1:int.equation.d
            # TODO initial guess
            # int.solver.x[2*(sol.nd*(i-1)+k-1)+1] = int.q[k]
            # int.solver.x[2*(params.d*(i-1)+k-1)+1] = 0.
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
