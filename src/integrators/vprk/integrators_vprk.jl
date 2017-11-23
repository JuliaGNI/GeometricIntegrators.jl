
"Parameters for right-hand side function of variational partitioned Runge-Kutta methods."
mutable struct NonlinearFunctionParametersVPRK{DT,TT,ΑT,FT,D,S} <: AbstractNonlinearFunctionParametersVPRK{DT,TT,ΑT,FT,D,S}
    α::ΑT
    f::FT

    Δt::TT

    t_q::CoefficientsRK{TT}
    t_p::CoefficientsRK{TT}
    d_v::Vector{TT}

    t::TT
    q::Vector{DT}
    p::Vector{DT}
end

function NonlinearFunctionParametersVPRK(DT, D, α::AT, f::FT, Δt::TT, t_q, t_p, d_v) where {TT,AT,FT}
    @assert t_q.s == t_p.s
    NonlinearFunctionParametersVPRK{DT,TT,AT,FT,D,t_q.s}(α, f, Δt, t_q, t_p, d_v, 0, zeros(DT,D), zeros(DT,D))
end


"Compute stages of variational partitioned Runge-Kutta methods."
@generated function function_stages!(x::Vector{ST}, b::Vector{ST}, params::NonlinearFunctionParametersVPRK{DT,TT,ΑT,FT,D,S}) where {ST,DT,TT,ΑT,FT,D,S}
    cache = NonlinearFunctionCacheVPRK{ST}(D, S)

    quote
        @assert length(x) == length(b)

        compute_stages!(x, $cache.Q, $cache.V, $cache.P, $cache.F, params)
        compute_rhs_vprk!(b, $cache.P, $cache.F, params)
        compute_rhs_vprk_correction!(b, $cache.V, params)
    end
end


"Variational partitioned Runge-Kutta integrator."
struct IntegratorVPRK{DT,TT,ΑT,FT,GT,VT,FPT,ST,IT} <: AbstractIntegratorVPRK{DT,TT}
    equation::IODE{DT,TT,ΑT,FT,GT,VT}
    tableau::TableauVPRK{TT}
    Δt::TT

    params::FPT
    solver::ST
    iguess::InitialGuessPODE{DT,TT,VT,FT,IT}

    q::Vector{Vector{Double{DT}}}
    p::Vector{Vector{Double{DT}}}

    cache::NonlinearFunctionCacheVPRK{DT}
end

function IntegratorVPRK(equation::IODE{DT,TT,ΑT,FT,GT,VT}, tableau::TableauVPRK{TT}, Δt::TT) where {DT,TT,ΑT,FT,GT,VT}
    D = equation.d
    M = equation.n
    S = tableau.s

    N = D*S

    if isdefined(tableau, :d)
        d_v = tableau.d
    else
        d_v = DT[]
    end

    # create solution vector for nonlinear solver
    x = zeros(DT,N)

    # create solution vectors
    q = Array{Vector{Double{DT}}}(M)
    p = Array{Vector{Double{DT}}}(M)

    for i in 1:M
        q[i] = zeros(Double{DT},D)
        p[i] = zeros(Double{DT},D)
    end

    # create cache for internal stage vectors and update vectors
    cache = NonlinearFunctionCacheVPRK{DT}(D,S)

    # create params
    params = NonlinearFunctionParametersVPRK(DT, D, equation.α, equation.f, Δt,
                                             tableau.q, tableau.p, d_v)

    # create rhs function for nonlinear solver
    function_stages = (x,b) -> function_stages!(x, b, params)

    # create nonlinear solver
    solver = get_config(:nls_solver)(x, function_stages)

    # create initial guess
    iguess = InitialGuessPODE(get_config(:ig_interpolation), equation, Δt)

    # create integrator
    IntegratorVPRK{DT, TT, ΑT, FT, GT, VT, typeof(params), typeof(solver), typeof(iguess.int)}(
                                        equation, tableau, Δt, params, solver, iguess,
                                        q, p, cache)
end


IntegratorVPRKpNone = IntegratorVPRK


function initialize!(int::Union{IntegratorVPRK{DT,TT}, IntegratorVPRKpMidpoint{DT,TT}, IntegratorVPRKpSymmetric{DT,TT}}, sol::SolutionPDAE, m::Int) where {DT,TT}
    @assert m ≥ 1
    @assert m ≤ sol.ni

    # copy initial conditions from solution
    get_initial_conditions!(sol, int.q[m], int.p[m], m)

    # initialise initial guess
    initialize!(int.iguess, m, sol.t[0], int.q[m], int.p[m])
end


function initial_guess!(int::IntegratorVPRK, m::Int)
    for i in 1:int.tableau.s
        evaluate!(int.iguess, m, int.cache.y, int.cache.z, int.cache.v, int.tableau.q.c[i], int.tableau.p.c[i])
        for k in 1:int.equation.d
            int.solver.x[int.equation.d*(i-1)+k] = int.cache.v[k]
        end
    end
end

"Integrate ODE with variational partitioned Runge-Kutta integrator."
function integrate_step!(int::IntegratorVPRK{DT,TT}, sol::SolutionPDAE{DT,TT}, m::Int, n::Int) where {DT,TT}
    @assert m ≥ 1
    @assert m ≤ sol.ni

    @assert n ≥ 1
    @assert n ≤ sol.ntime

    # set time for nonlinear solver
    int.params.t  = sol.t[0] + (n-1)*int.Δt
    int.params.q .= int.q[m]
    int.params.p .= int.p[m]

    # compute initial guess
    initial_guess!(int, m)

    # call nonlinear solver
    solve!(int.solver)

    # print solver status
    print_solver_status(int.solver.status, int.solver.params, n)

    # check if solution contains NaNs or error bounds are violated
    check_solver_status(int.solver.status, int.solver.params, n)

    # compute vector fields at internal stages
    compute_stages!(int.solver.x, int.cache.Q, int.cache.V, int.cache.P, int.cache.F, int.params)

    # compute final update
    update_solution!(int.q[m], int.cache.V, int.tableau.q.b, int.tableau.q.b̂, int.Δt)
    update_solution!(int.p[m], int.cache.F, int.tableau.p.b, int.tableau.p.b̂, int.Δt)

    # copy solution to initial guess for next time step
    update!(int.iguess, m, sol.t[0] + n*int.Δt, int.q[m], int.p[m])

    # take care of periodic solutions
    cut_periodic_solution!(int.q[m], int.equation.periodicity)

    # copy to solution
    copy_solution!(sol, int.q[m], int.p[m], n, m)
end
