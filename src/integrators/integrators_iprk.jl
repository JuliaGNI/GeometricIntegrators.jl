
"Parameters for right-hand side function of special partitioned Runge-Kutta methods."
immutable NonlinearFunctionParametersIPRK{T} <: NonlinearFunctionParameters{T}
    f::Function
    g::Function
    Δt::T

    d::Int
    s::Int
    a_q::Matrix{T}
    a_p::Matrix{T}

    q::Vector{T}
    p::Vector{T}
    y::Vector{T}
    z::Vector{T}

    Q::Matrix{T}
    V::Matrix{T}
    P::Matrix{T}
    F::Matrix{T}
    Y::Matrix{T}
    Z::Matrix{T}

    tQ::Vector{T}
    tV::Vector{T}
    tP::Vector{T}
    tF::Vector{T}

    function NonlinearFunctionParametersIPRK(f, g, Δt, d, s, a_q, a_p)
        # create solution vectors
        q = zeros(T,d)
        p = zeros(T,d)
        y = zeros(T,d)
        z = zeros(T,d)

        # create internal stage vectors
        Q = zeros(T,d,s)
        V = zeros(T,d,s)
        P = zeros(T,d,s)
        F = zeros(T,d,s)
        Y = zeros(T,d,s)
        Z = zeros(T,d,s)

        # create temporary vectors
        tQ = zeros(T,d)
        tV = zeros(T,d)
        tP = zeros(T,d)
        tF = zeros(T,d)

        new(f, g, Δt, d, s, a_q, a_p, q, p, y, z, Q, V, P, F, Y, Z, tQ, tV, tP, tF)
    end
end

"Compute stages of special partitioned Runge-Kutta methods."
function function_stages!{T}(y::Vector{T}, b::Vector{T}, params::NonlinearFunctionParametersIPRK{T})
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
        simd_copy_xy_first!(params.tQ, params.Q, i)
        simd_copy_xy_first!(params.tV, params.V, i)
        params.f(params.tQ, params.tV, params.tP)
        params.g(params.tQ, params.tV, params.tF)
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


"Special partitioned Runge-Kutta integrator."
immutable IntegratorIPRK{T, ST} <: Integrator{T}
    equation::SPODE{T}
    tableau::TableauIPRK{T}
    Δt::T

    solver::ST

    q::Array{T,1}
    p::Array{T,1}
    y::Array{T,1}
    z::Array{T,1}
    Q::Array{T,2}
    V::Array{T,2}
    P::Array{T,2}
    F::Array{T,2}
    Y::Array{T,2}
    Z::Array{T,2}
end

function IntegratorIPRK{T}(equation::SPODE{T}, tableau::TableauIPRK{T}, Δt::T)
    D = equation.d
    S = tableau.s

    # create solution vector for internal stages / nonlinear solver
    z = zeros(T, 2*D*S)

    # create params
    params = NonlinearFunctionParametersIPRK{T}(equation.f, equation.g, Δt, D, S, tableau.a_q, tableau.a_p)

    # create solver
    solver = QuasiNewtonSolver(z, params)
    # TODO allow for other nonlinear solvers based on constructor argument

    # create integrator
    IntegratorIPRK{T, typeof(solver)}(equation, tableau, Δt, solver,
                                      params.q, params.p, params.y, params.z,
                                      params.Q, params.V, params.P, params.F,
                                      params.Y, params.Z)
end


"Integrate ODE with special partitioned Runge-Kutta integrator."
function integrate!(int::IntegratorIPRK, sol::SolutionPODE)
    local nt::Int

    # copy initial conditions from solution
    simd_copy_xy_first!(int.q, sol.q, 1)
    simd_copy_xy_first!(int.p, sol.p, 1)

    for n in 1:sol.ntime
        # compute initial guess
        # TODO
        for i in 1:int.tableau.s
            for k in 1:int.equation.d
                int.solver.x[2*(sol.d*(i-1)+k-1)+1] = int.q[k]
                # TODO initial guess for velocity
                # int.solver.x[2*(params.d*(i-1)+k-1)+1] = 0.
            end
        end

        # call nonlinear solver
        solve!(int.solver)

        if (int.solver.status.rₐ > int.solver.params.atol² &&
            int.solver.status.rᵣ > int.solver.params.rtol  &&
            int.solver.status.rₛ > int.solver.params.stol²)||
            int.solver.status.i >= int.solver.params.nmax
            println(int.solver.status.i, ", ", int.solver.status.rₐ,", ",  int.solver.status.rᵣ,", ",  int.solver.status.rₛ)
        end

        # compute final update
        simd_mult!(int.y, int.V, int.tableau.b_q)
        simd_mult!(int.z, int.F, int.tableau.b_p)
        simd_axpy!(int.Δt, int.y, int.q)
        simd_axpy!(int.Δt, int.z, int.p)

        # copy to solution
        if mod(n, sol.nsave) == 0
            nt = div(n, sol.nsave)+1
            simd_copy_yx_first!(int.q, sol.q, nt)
            simd_copy_yx_first!(int.p, sol.p, nt)
            sol.t[nt] = sol.t[1] + n * int.Δt
        end
    end
    nothing
end
