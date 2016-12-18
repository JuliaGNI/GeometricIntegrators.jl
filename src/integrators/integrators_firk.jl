
"Parameters for right-hand side function of fully implicit Runge-Kutta methods."
type NonlinearFunctionParametersFIRK{DT,TT,FT} <: NonlinearFunctionParameters{DT}
    f::FT
    Δt::TT

    d::Int
    s::Int

    a::Matrix{TT}
    c::Vector{TT}

    t::TT
    x::Vector{DT}
    y::Vector{DT}

    X::Matrix{DT}
    Y::Matrix{DT}
    F::Matrix{DT}

    tX::Vector{DT}
    tF::Vector{DT}

    function NonlinearFunctionParametersFIRK(f, Δt, d, s, a, c)
        # create solution vectors
        x = zeros(DT,d)
        y = zeros(DT,d)

        # create internal stage vectors
        X = zeros(DT,d,s)
        Y = zeros(DT,d,s)
        F = zeros(DT,d,s)

        # create temporary vectors
        tX = zeros(DT,d)
        tF = zeros(DT,d)

        new(f, Δt, d, s, a, c, 0, x, y, X, Y, F, tX, tF)
    end
end

"Compute stages of fully implicit Runge-Kutta methods."
function function_stages!{DT,TT,FT}(y::Vector{DT}, b::Vector{TT}, params::NonlinearFunctionParametersFIRK{DT,TT,FT})
    local tᵢ::TT

    for i in 1:params.s
        tᵢ = params.t + params.Δt * params.c[i]

        for k in 1:params.d
            # copy y to Y
            params.Y[k,i] = y[params.d*(i-1)+k]

            # compute X
            params.X[k,i] = params.x[k] + params.Δt * params.Y[k,i]
        end

        # compute f(X)
        simd_copy_xy_first!(params.tX, params.X, i)
        params.f(tᵢ, params.tX, params.tF)
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
immutable IntegratorFIRK{DT, TT, FT, ST, IT} <: Integrator{DT,TT}
    equation::ODE{DT,TT,FT}
    tableau::TableauFIRK{TT}
    Δt::TT

    solver::ST
    iguess::InitialGuess{DT, TT, FT, IT}

    x::Array{DT,1}
    y::Array{DT,1}
    X::Array{DT,2}
    Y::Array{DT,2}
    F::Array{DT,2}
end

function IntegratorFIRK{DT,TT,FT}(equation::ODE{DT,TT,FT}, tableau::TableauFIRK{TT}, Δt::TT;
                              nonlinear_solver=QuasiNewtonSolver,
                              interpolation=HermiteInterpolation{DT})
    D = equation.d
    S = tableau.s

    # create solution vector for internal stages / nonlinear solver
    z = zeros(DT, D*S)

    # create params
    params = NonlinearFunctionParametersFIRK{DT,TT,FT}(equation.v, Δt, D, S, tableau.a, tableau.c)

    # create solver
    solver = nonlinear_solver(z, params)

    # create initial guess
    iguess = InitialGuess(interpolation, equation, Δt)

    # create integrator
    IntegratorFIRK{DT, TT, FT, typeof(solver), typeof(iguess.int)}(equation, tableau, Δt, solver, iguess, params.x, params.y, params.X, params.Y, params.F)
end


"Integrate ODE with fully implicit Runge-Kutta integrator."
function integrate!{DT, TT, FT, ST, IT}(int::IntegratorFIRK{DT, TT, FT, ST, IT}, sol::SolutionODE{DT,TT})
    # loop over initial conditions
    for m in 1:sol.n0
        # copy initial conditions from solution
        simd_copy_xy_first!(int.x, sol, 0, m)

        # initialise initial guess
        initialize!(int.iguess, sol.t[0], int.x)

        for n in 1:sol.ntime
            # set time for nonlinear solver
            int.solver.Fparams.t = sol.t[n]

            # copy previous solution to initial guess
            update!(int.iguess, sol.t[n], int.x)

            # compute initial guess for internal stages
            for i in 1:int.tableau.s
                evaluate(int.iguess, int.y, int.tableau.c[i])
                for k in 1:int.equation.d
                    int.solver.x[int.equation.d*(i-1)+k] = int.y[k]
                end
            end

            # call nonlinear solver
            solve!(int.solver)

            if !solverStatusOK(int.solver.status, int.solver.params)
                println(int.solver.status)
            end

            # compute final update
            simd_mult!(int.y, int.F, int.tableau.b)
            simd_axpy!(int.Δt, int.y, int.x)

            # copy to solution
            if mod(n, sol.nsave) == 0
                simd_copy_yx_first!(int.x, sol, div(n, sol.nsave), m)
            end
        end
    end
    nothing
end

"Integrate partitioned ODE with fully implicit Runge-Kutta integrator."
function integrate!(int::IntegratorFIRK, s::SolutionPODE)
    # TODO
end
