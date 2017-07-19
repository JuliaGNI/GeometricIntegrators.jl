
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

function NonlinearFunctionParametersFIRK(DT, D, v::VT, Δt::TT, tab) where {TT,VT}
    NonlinearFunctionParametersFIRK{DT,TT,VT,D,tab.s}(v, Δt, tab.a, tab.â, tab.c, 0, zeros(DT,D))
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

@generated function compute_stages!(x::Vector{ST}, Q::Matrix{ST}, V::Matrix{ST}, Y::Matrix{ST},
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
        compute_stages!(x, $cache.Q, $cache.V, $cache.Y, params)

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
                b[D*(i-1)+k] = - $cache.Y[k,i] + (y1 + y2)
            end
        end
    end
end


"Fully implicit Runge-Kutta integrator."
struct IntegratorFIRK{DT, TT, FT, SPT, ST, IT, N} <: Integrator{DT,TT}
    equation::ODE{DT,TT,FT,N}
    tableau::TableauFIRK{TT}
    Δt::TT

    params::SPT
    solver::ST
    iguess::InitialGuessODE{DT, TT, FT, IT}

    q::Vector{Vector{Double{DT}}}

    v::Vector{DT}
    y::Vector{DT}

    Q::Matrix{DT}
    V::Matrix{DT}
    Y::Matrix{DT}
end

function IntegratorFIRK(equation::ODE{DT,TT,FT,N}, tableau::TableauFIRK{TT}, Δt::TT;
                        interpolation=HermiteInterpolation{DT}) where {DT,TT,FT,N}
    D = equation.d
    M = equation.n
    S = tableau.q.s

    # create solution vector for internal stages / nonlinear solver
    x = zeros(DT, D*S)

    # create solution vectors
    q = Array{Vector{Double{DT}}}(M)
    for i in 1:M
        q[i] = zeros(Double{DT},D)
    end

    # create velocity and update vector
    v = zeros(DT,D)
    y = zeros(DT,D)

    # create internal stage vectors
    Q = zeros(DT,D,S)
    V = zeros(DT,D,S)
    Y = zeros(DT,D,S)

    # create params
    params = NonlinearFunctionParametersFIRK(DT, D, equation.v, Δt, tableau.q)

    # create rhs function for nonlinear solver
    function_stages = (x,b) -> function_stages!(x, b, params)

    # create solver
    solver = get_config(:nls_solver)(x, function_stages)

    # create initial guess
    iguess = InitialGuessODE(interpolation, equation, Δt; periodicity=equation.periodicity)

    # create integrator
    IntegratorFIRK{DT, TT, FT, typeof(params), typeof(solver), typeof(iguess.int), N}(
                                        equation, tableau, Δt, params, solver, iguess,
                                        q, v, y, Q, V, Y)
end


function initialize!(int::IntegratorFIRK{DT,TT}, sol::SolutionODE, m::Int) where {DT,TT}
    @assert m ≥ 1
    @assert m ≤ sol.ni

    # copy initial conditions from solution
    get_initial_conditions!(sol, int.q[m], m)

    # initialise initial guess
    initialize!(int.iguess, m, sol.t[0], int.q[m])
end

function initial_guess!(int::IntegratorFIRK, m::Int)
    # compute initial guess for internal stages
    for i in 1:int.tableau.q.s
        evaluate!(int.iguess, m, int.y, int.v, int.tableau.q.c[i])
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
end


"Integrate ODE with fully implicit Runge-Kutta integrator."
function integrate_step!(int::IntegratorFIRK{DT, TT, FT, SPT, ST, IT}, sol::SolutionODE{DT,TT,N}, m::Int, n::Int) where {DT,TT,FT,SPT,ST,IT,N}
    @assert m ≥ 1
    @assert m ≤ sol.ni

    @assert n ≥ 1
    @assert n ≤ sol.ntime

    # set time for nonlinear solver
    int.params.t  = sol.t[0] + (n-1)*int.Δt
    int.params.q .= int.q[m]

    # compute initial guess
    initial_guess!(int, m)

    # call nonlinear solver
    solve!(int.solver)

    # print solver status
    print_solver_status(int.solver.status, int.solver.params, n)

    # check if solution contains NaNs or error bounds are violated
    check_solver_status(int.solver.status, int.solver.params, n)

    # compute vector field at internal stages
    compute_stages!(int.solver.x, int.Q, int.V, int.Y, int.params)

    # compute final update
    update_solution!(int.q[m], int.V, int.tableau.q.b, int.tableau.q.b̂, int.Δt)

    # copy solution to initial guess
    update!(int.iguess, m, sol.t[0] + n*int.Δt, int.q[m])

    # take care of periodic solutions
    cut_periodic_solution!(int.q[m], int.equation.periodicity)

    # copy to solution
    copy_solution!(sol, int.q[m], n, m)
end
