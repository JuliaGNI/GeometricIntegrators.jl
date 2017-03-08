
"Parameters for right-hand side function of variational partitioned Runge-Kutta methods."
type NonlinearFunctionParametersVPRK{DT,TT,ΑT,FT,D,S} <: AbstractNonlinearFunctionParametersVPRK{DT,TT,ΑT,FT,D,S}
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

function NonlinearFunctionParametersVPRK{DT,TT,AT,FT}(α::AT, f::FT, Δt::TT, t_q, t_p, d_v, q::Vector{DT}, p::Vector{DT})
    @assert t_q.s == t_p.s
    @assert length(q) == length(p)
    NonlinearFunctionParametersVPRK{DT,TT,AT,FT,length(q),t_q.s}(α, f, Δt, t_q, t_p, d_v, 0, q, p)
end


"Compute stages of variational partitioned Runge-Kutta methods."
@generated function function_stages!{ST,DT,TT,ΑT,FT,D,S}(x::Vector{ST}, b::Vector{ST}, params::NonlinearFunctionParametersVPRK{DT,TT,ΑT,FT,D,S})
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
immutable IntegratorVPRK{DT,TT,ΑT,FT,GT,VT,FPT,ST,IT} <: AbstractIntegratorVPRK{DT,TT}
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

function IntegratorVPRK{DT,TT,ΑT,FT,GT,VT}(equation::IODE{DT,TT,ΑT,FT,GT,VT}, tableau::TableauVPRK{TT}, Δt::TT;
                                        nonlinear_solver=DEFAULT_NonlinearSolver,
                                        nmax=DEFAULT_nmax, atol=DEFAULT_atol, rtol=DEFAULT_rtol, stol=DEFAULT_stol,
                                        interpolation=HermiteInterpolation{DT})
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
    solver = nonlinear_solver(x, function_stages; nmax=nmax, atol=atol, rtol=rtol, stol=stol)

    # create initial guess
    iguess = InitialGuessIODE(interpolation, equation, Δt)

    # create integrator
    IntegratorVPRK{DT, TT, ΑT, FT, GT, VT, typeof(params), typeof(solver), typeof(iguess.int)}(
                                        equation, tableau, Δt, params, solver, iguess,
                                        q, p, qₑᵣᵣ, pₑᵣᵣ, cache)
end


"Integrate ODE with variational partitioned Runge-Kutta integrator."
function integrate!{DT,TT,ΑT,FT,GT,VT,N}(int::IntegratorVPRK{DT,TT,ΑT,FT,GT,VT}, sol::SolutionPDAE{DT,TT,N})
    # loop over initial conditions
    for m in 1:sol.ni
        # copy initial conditions from solution
        get_initial_conditions!(sol, int.q, int.p, m)

        # initialise initial guess
        initialize!(int.iguess, sol.t[0], int.q, int.p)

        for n in 1:sol.ntime
            # set time for nonlinear solver
            int.params.t = sol.t[0] + (n-1)*int.Δt

            # copy previous solution to initial guess
            update!(int.iguess, sol.t[0] + n*int.Δt, int.q, int.p)

            # compute initial guess
            for i in 1:int.tableau.s
                evaluate!(int.iguess, int.cache.y, int.cache.z, int.cache.v, int.tableau.q.c[i], int.tableau.p.c[i])
                for k in 1:int.equation.d
                    int.solver.x[int.equation.d*(i-1)+k] = int.cache.v[k]
                end
            end

            # call nonlinear solver
            solve!(int.solver)

            if !solverStatusOK(int.solver.status, int.solver.params)
                println(int.solver.status, ", it=", n)
            end

            if isnan(int.solver.status.rₐ)
                break
            end

            # compute final update
            compute_stages_vprk!(int.solver.x, int.cache.Q, int.cache.V, int.cache.P, int.cache.F, int.params)
            update_solution!(int, int.cache)
            cut_periodic_solution!(int)

            # copy to solution
            copy_solution!(sol, int.q, int.p, n, m)
        end
    end
end
