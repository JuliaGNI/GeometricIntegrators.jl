
"Parameters for right-hand side function of fully implicit Runge-Kutta methods."
type NonlinearFunctionParametersFIRK{DT,TT,VT,D,S} <: NonlinearFunctionParameters{DT}
    v::VT
    Δt::TT

    a::Matrix{TT}
    c::Vector{TT}

    t::TT

    q::Vector{DT}
end

function NonlinearFunctionParametersFIRK{DT,TT,VT}(v::VT, Δt::TT, tab, q::Vector{DT})
    NonlinearFunctionParametersFIRK{DT,TT,VT,length(q),tab.s}(v, Δt, tab.a, tab.c, 0, q)
end

immutable NonlinearFunctionCacheFIRK{DT}
    Q::Matrix{DT}
    V::Matrix{DT}
    Y::Matrix{DT}

    tQ::Vector{DT}
    tV::Vector{DT}

    function NonlinearFunctionCacheFIRK(d, s)

        # create internal stage vectors
        Q = zeros(DT,d,s)
        V = zeros(DT,d,s)
        Y = zeros(DT,d,s)

        # create temporary vectors
        tQ = zeros(DT,d)
        tV = zeros(DT,d)

        new(Q, V, Y, tQ, tV)
    end
end

function compute_stages_firk!{DT,TT,VT}(x::Vector{DT}, Q::Matrix{DT}, V::Matrix{DT}, Y::Matrix{DT},
                                           q::Vector, a::Matrix{TT}, c::Vector{TT}, Δt::TT, t::TT, v::VT,
                                           tQ::Vector{DT}, tV::Vector{DT})

    local d::Int = length(q)
    local s::Int = length(c)
    local tᵢ::TT

    @assert d == length(tQ) == length(tV)
    @assert s == size(a,1) == size(a,2)
    @assert d == size(Q,1) == size(V,1) == size(Y,1)
    @assert s == size(Q,2) == size(V,2) == size(Y,2)

    # copy x to Y and compute Q = q + Δt Y
    for i in 1:size(Y,2)
        for k in 1:size(Y,1)
            Y[k,i] = x[d*(i-1)+k]
            Q[k,i] = q[k] + Δt * Y[k,i]
        end
    end

    # compute V = v(Q)
    for i in 1:s
        tᵢ = t + Δt * c[i]
        simd_copy_xy_first!(tQ, Q, i)
        v(tᵢ, tQ, tV)
        simd_copy_yx_first!(tV, V, i)
    end

    nothing
end

"Compute stages of fully implicit Runge-Kutta methods."
@generated function function_stages!{ST,DT,TT,VT,D,S}(x::Vector{ST}, b::Vector{ST}, params::NonlinearFunctionParametersFIRK{DT,TT,VT,D,S})

    cache = NonlinearFunctionCacheFIRK{ST}(D, S)

    function_stages = quote

        compute_stages_firk!(x, $cache.Q, $cache.V, $cache.Y, params.q, params.a, params.c, params.Δt, params.t, params.v, $cache.tQ, $cache.tV)

        # compute b = - (Y-AV)
        for i in 1:S
            for k in 1:D
                b[D*(i-1)+k] = - $cache.Y[k,i]
                for j in 1:S
                    b[D*(i-1)+k] += params.a[i,j] * $cache.V[k,j]
                end
            end
        end
    end

    return function_stages
end


"Fully implicit Runge-Kutta integrator."
immutable IntegratorFIRK{DT, TT, FT, SPT, ST, IT} <: Integrator{DT,TT}
    equation::ODE{DT,TT,FT}
    tableau::TableauFIRK{TT}
    Δt::TT

    Sparams::SPT
    solver::ST
    iguess::InitialGuessODE{DT, TT, FT, IT}

    q::Vector{DT}
    v::Vector{DT}
    y::Vector{DT}

    Q::Matrix{DT}
    V::Matrix{DT}
    Y::Matrix{DT}

    tQ::Vector{DT}
    tV::Vector{DT}
end

function IntegratorFIRK{DT,TT,FT}(equation::ODE{DT,TT,FT}, tableau::TableauFIRK{TT}, Δt::TT;
                                  nonlinear_solver=DEFAULT_NonlinearSolver,
                                  nmax=DEFAULT_nmax, atol=DEFAULT_atol, rtol=DEFAULT_rtol, stol=DEFAULT_stol,
                                  interpolation=HermiteInterpolation{DT})
    D = equation.d
    S = tableau.q.s

    # create solution vector for internal stages / nonlinear solver
    x = zeros(DT, D*S)

    # create solution vectors
    q = zeros(DT,D)
    v = zeros(DT,D)

    # create update vector
    y = zeros(DT,D)

    # create internal stage vectors
    Q = zeros(DT,D,S)
    V = zeros(DT,D,S)
    Y = zeros(DT,D,S)

    # create temporary vectors
    tQ = zeros(DT,D)
    tV = zeros(DT,D)

    # create params
    params = NonlinearFunctionParametersFIRK(equation.v, Δt, tableau.q, q)

    # create rhs function for nonlinear solver
    function_stages = (x,b) -> function_stages!(x, b, params)

    # create solver
    solver = nonlinear_solver(x, function_stages; nmax=nmax, atol=atol, rtol=rtol, stol=stol)

    # create initial guess
    iguess = InitialGuessODE(interpolation, equation, Δt)

    # create integrator
    IntegratorFIRK{DT, TT, FT, typeof(params), typeof(solver), typeof(iguess.int)}(
                                        equation, tableau, Δt, params, solver, iguess,
                                        q, v, y, Q, V, Y, tQ, tV)
end


"Integrate ODE with fully implicit Runge-Kutta integrator."
function integrate!{DT,TT,FT,SPT,ST,IT,N}(int::IntegratorFIRK{DT, TT, FT, SPT, ST, IT}, sol::SolutionODE{DT,TT,N}, m1::Int, m2::Int)
    @assert m1 ≥ 1
    @assert m2 ≤ sol.ni

    # loop over initial conditions
    for m in m1:m2
        # copy initial conditions from solution
        get_initial_conditions!(sol, int.q, m)

        # initialise initial guess
        initialize!(int.iguess, sol.t[0], int.q)

        for n in 1:sol.ntime
            # set time for nonlinear solver
            int.Sparams.t = sol.t[n-1]

            # copy previous solution to initial guess
            update!(int.iguess, sol.t[n], int.q)

            # compute initial guess for internal stages
            for i in 1:int.tableau.q.s
                evaluate!(int.iguess, int.y, int.v, int.tableau.q.c[i])
                for k in 1:int.equation.d
                    int.solver.x[int.equation.d*(i-1)+k] = int.v[k]
                end
            end

            # call nonlinear solver
            solve!(int.solver)

            if !solverStatusOK(int.solver.status, int.solver.params)
                println(int.solver.status)
            end

            compute_stages_firk!(int.solver.x, int.Q, int.V, int.Y, int.q,
                                 int.tableau.q.a, int.tableau.q.c, int.Δt, sol.t[n-1],
                                 int.equation.v, int.tQ, int.tV)

            # compute final update
            simd_mult!(int.y, int.V, int.tableau.q.b)
            simd_axpy!(int.Δt, int.y, int.q)

            # take care of periodic solutions
            for k in 1:int.equation.d
                if int.equation.periodicity[k] ≠ 0
                    int.q[k] = mod(int.q[k], int.equation.periodicity[k])
                end
            end

            # copy to solution
            copy_solution!(sol, int.q, n, m)
        end
    end
    nothing
end

"Integrate partitioned ODE with fully implicit Runge-Kutta integrator."
function integrate!(int::IntegratorFIRK, s::SolutionPODE)
    # TODO
end
