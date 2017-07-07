
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

function NonlinearFunctionParametersVPRK(α::AT, f::FT, Δt::TT, t_q, t_p, d_v, q::Vector{DT}, p::Vector{DT}) where {DT,TT,AT,FT}
    @assert t_q.s == t_p.s
    @assert length(q) == length(p)
    NonlinearFunctionParametersVPRK{DT,TT,AT,FT,length(q),t_q.s}(α, f, Δt, t_q, t_p, d_v, 0, q, p)
end


"Compute stages of variational partitioned Runge-Kutta methods."
@generated function function_stages!(x::Vector{ST}, b::Vector{ST}, params::NonlinearFunctionParametersVPRK{DT,TT,ΑT,FT,D,S}) where {ST,DT,TT,ΑT,FT,D,S}
    cache = NonlinearFunctionCacheVPRK{ST}(D, S)

    function_stages = quote
        @assert length(x) == length(b)

        compute_stages_vprk!(x, $cache.Q, $cache.V, $cache.P, $cache.F, params)
        compute_rhs_vprk!(b, $cache.P, $cache.F, params)
        compute_rhs_vprk_correction!(b, $cache.V, params)
    end

    return function_stages
end


"Variational partitioned Runge-Kutta integrator."
struct IntegratorVPRK{DT,TT,ΑT,FT,GT,VT,FPT,ST,IT} <: AbstractIntegratorVPRK{DT,TT}
    equation::IODE{DT,TT,ΑT,FT,GT,VT}
    tableau::TableauVPRK{TT}
    Δt::TT

    params::FPT
    solver::ST
    iguess::InitialGuessIODE{DT,TT,VT,FT,IT}

    q::Array{DT,1}
    p::Array{DT,1}

    qₑᵣᵣ::Vector{DT}
    pₑᵣᵣ::Vector{DT}

    cache::NonlinearFunctionCacheVPRK{DT}
end

function IntegratorVPRK(equation::IODE{DT,TT,ΑT,FT,GT,VT}, tableau::TableauVPRK{TT}, Δt::TT;
                        interpolation=HermiteInterpolation{DT}) where {DT,TT,ΑT,FT,GT,VT}
    D = equation.d
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
    q = zeros(DT,D)
    p = zeros(DT,D)

    # create compensated summation error vectors
    qₑᵣᵣ = zeros(DT,D)
    pₑᵣᵣ = zeros(DT,D)

    # create cache for internal stage vectors and update vectors
    cache = NonlinearFunctionCacheVPRK{DT}(D,S)

    # create params
    params = NonlinearFunctionParametersVPRK(equation.α, equation.f, Δt,
                                             tableau.q, tableau.p, d_v,
                                             q, p)

    # create rhs function for nonlinear solver
    function_stages = (x,b) -> function_stages!(x, b, params)

    # create nonlinear solver
    solver = get_config(:nls_solver)(x, function_stages)

    # create initial guess
    iguess = InitialGuessIODE(interpolation, equation, Δt; periodicity=equation.periodicity)

    # create integrator
    IntegratorVPRK{DT, TT, ΑT, FT, GT, VT, typeof(params), typeof(solver), typeof(iguess.int)}(
                                        equation, tableau, Δt, params, solver, iguess,
                                        q, p, qₑᵣᵣ, pₑᵣᵣ, cache)
end



function initialize!(int::Union{IntegratorVPRK{DT,TT}, IntegratorVPRKpMidpoint{DT,TT}, IntegratorVPRKpSymmetric{DT,TT}}, sol::SolutionPDAE, m::Int) where {DT,TT}
    @assert m ≥ 1
    @assert m ≤ sol.ni

    # copy initial conditions from solution
    get_initial_conditions!(sol, int.q, int.p, m)

    # initialise initial guess
    initialize!(int.iguess, sol.t[0], int.q, int.p)

    # reset compensated summation error
    int.qₑᵣᵣ .= 0
    int.pₑᵣᵣ .= 0
end


"Integrate ODE with variational partitioned Runge-Kutta integrator."
function integrate_step!(int::IntegratorVPRK{DT,TT,ΑT,FT,GT,VT}, sol::SolutionPDAE{DT,TT,N}, m::Int, n::Int) where {DT,TT,ΑT,FT,GT,VT,N}
    # set time for nonlinear solver
    int.params.t = sol.t[0] + (n-1)*int.Δt

    # compute initial guess
    for i in 1:int.tableau.s
        evaluate!(int.iguess, int.cache.y, int.cache.z, int.cache.v, int.tableau.q.c[i], int.tableau.p.c[i])
        for k in 1:int.equation.d
            int.solver.x[int.equation.d*(i-1)+k] = int.cache.v[k]
        end
    end

    # call nonlinear solver
    solve!(int.solver)

    # print solver status
    printSolverStatus(int.solver.status, int.solver.params, n)

    # if isnan(int.solver.status.rₐ)
    #     println("WARNING: Detected NaN in it=", n)
    #     break
    # end

    # compute final update
    compute_stages_vprk!(int.solver.x, int.cache.Q, int.cache.V, int.cache.P, int.cache.F, int.params)
    update_solution!(int, int.cache)

    # copy solution to initial guess for next time step
    update!(int.iguess, sol.t[0] + n*int.Δt, int.q, int.p)

    # take care of periodic solutions
    cut_periodic_solution!(int)

    # copy to solution
    copy_solution!(sol, int.q, int.p, n, m)
end
