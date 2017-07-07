
"Parameters for right-hand side function of fully implicit Runge-Kutta methods."
mutable struct NonlinearFunctionParametersFIRK{DT,TT,VT,D,S} <: NonlinearFunctionParameters{DT}
    v::VT
    Δt::TT

    a::Matrix{TT}
    â::Matrix{TT}
    c::Vector{TT}

    t::TT

    q::Vector{DT}
end

function NonlinearFunctionParametersFIRK(v::VT, Δt::TT, tab, q::Vector{DT}) where {DT,TT,VT}
    NonlinearFunctionParametersFIRK{DT,TT,VT,length(q),tab.s}(v, Δt, tab.a, tab.â, tab.c, 0, q)
end

struct NonlinearFunctionCacheFIRK{DT}
    Q::Matrix{DT}
    V::Matrix{DT}
    Y::Matrix{DT}

    function NonlinearFunctionCacheFIRK{DT}(d, s) where {DT}

        # create internal stage vectors
        Q = zeros(DT,d,s)
        V = zeros(DT,d,s)
        Y = zeros(DT,d,s)

        new(Q, V, Y)
    end
end

@generated function compute_stages_firk!(x::Vector{ST}, Q::Matrix{ST}, V::Matrix{ST}, Y::Matrix{ST},
                                         params::NonlinearFunctionParametersFIRK{DT,TT,VT,D,S}) where {ST,DT,TT,VT,D,S}

    tQ::Vector{ST} = zeros(ST,D)
    tV::Vector{ST} = zeros(ST,D)

    quote
        local tᵢ::TT

        @assert D == size(Q,1) == size(V,1) == size(Y,1)
        @assert S == size(Q,2) == size(V,2) == size(Y,2)

        # copy x to Y and compute Q = q + Δt Y
        for i in 1:size(Y,2)
            for k in 1:size(Y,1)
                Y[k,i] = x[D*(i-1)+k]
                Q[k,i] = params.q[k] + params.Δt * Y[k,i]
            end
        end

        # compute V = v(Q)
        for i in 1:S
            tᵢ = params.t + params.Δt * params.c[i]
            simd_copy_xy_first!($tQ, Q, i)
            params.v(tᵢ, $tQ, $tV)
            simd_copy_yx_first!($tV, V, i)
        end
    end
end

"Compute stages of fully implicit Runge-Kutta methods."
@generated function function_stages!(x::Vector{ST}, b::Vector{ST}, params::NonlinearFunctionParametersFIRK{DT,TT,VT,D,S}) where {ST,DT,TT,VT,D,S}

    cache = NonlinearFunctionCacheFIRK{ST}(D, S)

    quote
        compute_stages_firk!(x, $cache.Q, $cache.V, $cache.Y, params)

        local y1::ST
        local y2::ST

        # compute b = - (Y-AV)
        for i in 1:S
            for k in 1:D
                y1 = 0
                y2 = 0
                for j in 1:S
                    y1 += params.a[i,j] * $cache.V[k,j]
                    y2 += params.â[i,j] * $cache.V[k,j]
                end
                b[D*(i-1)+k] = - $cache.Y[k,i] + y1 + y2
            end
        end
    end
end


"Fully implicit Runge-Kutta integrator."
struct IntegratorFIRK{DT, TT, FT, SPT, ST, IT} <: Integrator{DT,TT}
    equation::ODE{DT,TT,FT}
    tableau::TableauFIRK{TT}
    Δt::TT

    params::SPT
    solver::ST
    iguess::InitialGuessODE{DT, TT, FT, IT}

    q::Vector{Double{DT}}

    v::Vector{DT}
    y::Vector{DT}

    Q::Matrix{DT}
    V::Matrix{DT}
    Y::Matrix{DT}
end

function IntegratorFIRK(equation::ODE{DT,TT,FT}, tableau::TableauFIRK{TT}, Δt::TT;
                        interpolation=HermiteInterpolation{DT}) where {DT,TT,FT}
    D = equation.d
    S = tableau.q.s

    # create solution vector for internal stages / nonlinear solver
    x = zeros(DT, D*S)

    # create solution vectors
    q = zeros(Double{DT},D)

    # create velocity and update vector
    v = zeros(DT,D)
    y = zeros(DT,D)

    # create internal stage vectors
    Q = zeros(DT,D,S)
    V = zeros(DT,D,S)
    Y = zeros(DT,D,S)

    # create params
    params = NonlinearFunctionParametersFIRK(equation.v, Δt, tableau.q, q)

    # create rhs function for nonlinear solver
    function_stages = (x,b) -> function_stages!(x, b, params)

    # create solver
    solver = get_config(:nls_solver)(x, function_stages)

    # create initial guess
    iguess = InitialGuessODE(interpolation, equation, Δt)

    # create integrator
    IntegratorFIRK{DT, TT, FT, typeof(params), typeof(solver), typeof(iguess.int)}(
                                        equation, tableau, Δt, params, solver, iguess,
                                        q, v, y, Q, V, Y)
end


function initialize!(int::IntegratorFIRK{DT,TT}, sol::SolutionODE, m::Int) where {DT,TT}
    @assert m ≥ 1
    @assert m ≤ sol.ni

    # copy initial conditions from solution
    get_initial_conditions!(sol, int.q, m)

    # initialise initial guess
    initialize!(int.iguess, sol.t[0], int.q)
end

"Integrate ODE with fully implicit Runge-Kutta integrator."
function integrate_step!(int::IntegratorFIRK{DT, TT, FT, SPT, ST, IT}, sol::SolutionODE{DT,TT,N}, m::Int, n::Int) where {DT,TT,FT,SPT,ST,IT,N}
    # set time for nonlinear solver
    int.params.t = sol.t[0] + (n-1)*int.Δt

    # copy previous solution to initial guess
    update!(int.iguess, int.params.t, int.q)

    # compute initial guess for internal stages
    for i in 1:int.tableau.q.s
        evaluate!(int.iguess, int.y, int.v, int.tableau.q.c[i])
        for k in 1:int.equation.d
            int.V[k,i] = int.v[k]
        end
    end
    for i in 1:int.tableau.q.s
        for k in 1:int.equation.d
            int.solver.x[int.equation.d*(i-1)+k] = 0
            for j in 1:int.tableau.q.s
                int.solver.x[int.equation.d*(i-1)+k] += int.tableau.q.a[i,j] * int.V[k,j]
            end
        end
    end

    # call nonlinear solver
    solve!(int.solver)

    # println(int.solver.status)
    if !solverStatusOK(int.solver.status, int.solver.params)
        println(int.solver.status)
    end

    # compute internal stages
    compute_stages_firk!(int.solver.x, int.Q, int.V, int.Y, int.params)

    # compute final update
    update_solution!(int.q, int.V, int.tableau.q.b, int.tableau.q.b̂, int.Δt)

    # take care of periodic solutions
    cut_periodic_solution!(int.q, int.equation.periodicity)

    # copy to solution
    copy_solution!(sol, int.q, n, m)
end
