
"Parameters for right-hand side function of fully implicit Runge-Kutta methods."
type NonlinearFunctionParametersFIRK{T} <: NonlinearFunctionParameters{T}
    f::Function
    Δt::T

    d::Int
    s::Int

    a::Matrix{T}
    c::Vector{T}

    t::T
    x::Vector{T}
    y::Vector{T}

    X::Matrix{T}
    Y::Matrix{T}
    F::Matrix{T}

    tX::Vector{T}
    tF::Vector{T}

    function NonlinearFunctionParametersFIRK(f, Δt, d, s, a, c)
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

        new(f, Δt, d, s, a, c, 0, x, y, X, Y, F, tX, tF)
    end
end

"Compute stages of fully implicit Runge-Kutta methods."
function function_stages!{T}(y::Vector{T}, b::Vector{T}, params::NonlinearFunctionParametersFIRK{T})
    local tᵢ::T

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
immutable IntegratorFIRK{T, ST, IT} <: Integrator{T}
    equation::ODE{T}
    tableau::TableauFIRK{T}
    Δt::T

    solver::ST
    iguess::InitialGuess{T, IT}

    x::Array{T,1}
    y::Array{T,1}
    X::Array{T,2}
    Y::Array{T,2}
    F::Array{T,2}
end

function IntegratorFIRK{T}(equation::ODE{T}, tableau::TableauFIRK{T}, Δt::T;
                           nonlinear_solver=QuasiNewtonSolver,
                           interpolation=HermiteInterpolation{T})
    D = equation.d
    S = tableau.s

    # create solution vector for internal stages / nonlinear solver
    z = zeros(T, D*S)

    # create params
    params = NonlinearFunctionParametersFIRK{T}(equation.f, Δt, D, S, tableau.a, tableau.c)

    # create solver
    solver = nonlinear_solver(z, params)

    # create initial guess
    iguess = InitialGuess(interpolation, equation, Δt)

    # create integrator
    IntegratorFIRK{T, typeof(solver), typeof(iguess.int)}(equation, tableau, Δt, solver, iguess, params.x, params.y, params.X, params.Y, params.F)
end


"Integrate ODE with fully implicit Runge-Kutta integrator."
function integrate!{T}(int::IntegratorFIRK{T}, sol::SolutionODE{T})
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

            if (int.solver.status.rₐ > int.solver.params.atol² &&
                int.solver.status.rᵣ > int.solver.params.rtol  &&
                int.solver.status.rₛ > int.solver.params.stol²)||
                int.solver.status.i >= int.solver.params.nmax
                println(int.solver.status.i, ", ", int.solver.status.rₐ,", ",
                        int.solver.status.rᵣ,", ", int.solver.status.rₛ)
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
