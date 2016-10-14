
"Parameters for right-hand side function of fully implicit Runge-Kutta methods."
immutable NonlinearFunctionParametersFIRK{T} <: NonlinearFunctionParameters{T}
    f::Function
    Δt::T

    d::Int
    s::Int
    a::Matrix{T}

    x::Vector{T}
    y::Vector{T}

    X::Matrix{T}
    Y::Matrix{T}
    F::Matrix{T}

    tX::Vector{T}
    tF::Vector{T}

    function NonlinearFunctionParametersFIRK(f, Δt, d, s, a)
        # create solution vectors
        x = zeros(T,d)
        y = zeros(T,d)

        # create internal stage vectors
        X = zeros(T,d,s)
        Y = zeros(T,d,s)
        F = zeros(T,d,s)

        # create temporary vectors
        tX = zeros(T,d)
        tF = zeros(T,d)

        new(f, Δt, d, s, a, x, y, X, Y, F, tX, tF)
    end
end

"Compute stages of fully implicit Runge-Kutta methods."
function function_stages!{T}(y::Vector{T}, b::Vector{T}, params::NonlinearFunctionParametersFIRK{T})
    for i in 1:params.s
        for k in 1:params.d
            # copy y to Y
            params.Y[k,i] = y[params.d*(i-1)+k]

            # compute X
            params.X[k,i] = params.x[k] + params.Δt * params.Y[k,i]
        end
        # compute f(X)
        simd_copy_xy_first!(params.tX, params.X, i)
        params.f(params.tX, params.tF)
        simd_copy_yx_first!(params.tF, params.F, i)
    end

    # compute b = - (Y-AF)
    for i in 1:params.s
        for k in 1:params.d
            b[params.d*(i-1)+k] = - params.Y[k,i]
            for j in 1:params.s
                b[params.d*(i-1)+k] += params.a[i,j] * params.F[k,j]
            end
        end
    end
end


"Fully implicit Runge-Kutta integrator."
immutable IntegratorFIRK{T, ST} <: Integrator{T}
    equation::ODE{T}
    tableau::TableauFIRK{T}
    Δt::T

    solver::ST

    x::Array{T,1}
    y::Array{T,1}
    X::Array{T,2}
    Y::Array{T,2}
    F::Array{T,2}
end

function IntegratorFIRK{T}(equation::ODE{T}, tableau::TableauFIRK{T}, Δt::T)
    D = equation.d
    S = tableau.s

    # create solution vector for internal stages / nonlinear solver
    z = zeros(T, D*S)

    # create params
    params = NonlinearFunctionParametersFIRK{T}(equation.f, Δt, D, S, tableau.a)

    # create solver
    solver = QuasiNewtonSolver(z, params)
    # TODO allow for other nonlinear solvers based on constructor argument

    # create integrator
    IntegratorFIRK{T, typeof(solver)}(equation, tableau, Δt, solver, params.x, params.y, params.X, params.Y, params.F)
end


"Integrate ODE with fully implicit Runge-Kutta integrator."
function integrate!(int::IntegratorFIRK, sol::SolutionODE)
    # copy initial conditions from solution
    simd_copy_xy_first!(int.x, sol, 0)

    for n in 1:sol.ntime
        # compute initial guess
        # TODO
        for i in 1:int.tableau.s
            for k in 1:int.equation.d
                int.solver.x[int.equation.d*(i-1)+k] = int.x[k]
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
        simd_mult!(int.y, int.F, int.tableau.b)
        simd_axpy!(int.Δt, int.y, int.x)

        # copy to solution
        if mod(n, sol.nsave) == 0
            simd_copy_yx_first!(int.x, sol, div(n, sol.nsave))
            sol.t[div(n, sol.nsave)+1] = sol.t[1] + n * int.Δt
        end
    end
    nothing
end

"Integrate partitioned ODE with fully implicit Runge-Kutta integrator."
function integrate!(int::IntegratorFIRK, s::SolutionPODE)
    # TODO
end
