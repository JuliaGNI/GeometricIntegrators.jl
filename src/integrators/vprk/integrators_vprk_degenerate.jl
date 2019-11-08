
# using Printf

"Parameters for right-hand side function of variational partitioned Runge-Kutta methods."
mutable struct ParametersVPRKdegenerate{DT, TT, ET <: IODE{DT,TT}, D, S} <: AbstractParametersVPRK{DT,TT,ET,D,S}
    equ::ET
    tab::TableauVPRK{TT}
    Δt::TT
    cache::NonlinearFunctionCacheVPRK{DT}

    t::TT
    q::Vector{DT}
    p::Vector{DT}
end

function ParametersVPRKdegenerate(equ::ET, tab::TableauVPRK{TT}, Δt::TT,
                                  cache::NonlinearFunctionCacheVPRK{DT}) where {DT, TT, ET <: IODE{DT,TT}}
    q = zeros(DT, equ.d)
    p = zeros(DT, equ.d)

    ParametersVPRKdegenerate{DT, TT, ET, equ.d, tab.s}(equ, tab, Δt, cache, 0, q, p)
end


"Compute solution of degenerate symplectic partitioned Runge-Kutta methods."
@generated function function_stages!(x::Vector{ST}, b::Vector{ST},
                params::ParametersVPRKdegenerate{DT,TT,ET,D,S}
            ) where {ST,DT,TT,ET,D,S}

    # create temporary vectors
    q̅ = zeros(ST,D)
    p̅ = zeros(ST,D)

    quote
        @assert length(x) == length(b)

        compute_solution!(x, $q̅, $p̅, params)

        # compute b = - [q̅-q-BV]
        for k in 1:div(D,2)
            b[k] = params.q[k] - $q̅[k]
            for i in 1:S
                b[k] += params.Δt * params.tab.q.b[i] * params.cache.V[k,i]
            end
        end

        # compute b = - [p̅-p-BF]
        for k in 1:div(D,2)
            b[div(D,2)+k] = params.p[k] - $p̅[k]
            for i in 1:S
                b[div(D,2)+k] += params.Δt * params.tab.p.b[i] * params.cache.F[k,i]
            end
        end
    end
end

function compute_solution!(
                x::Vector{ST}, q̅::Vector{ST}, p̅::Vector{ST},
                params::ParametersVPRKdegenerate{DT,TT,ET,D,S}
            ) where {ST,DT,TT,ET,D,S}

    @assert length(x) == length(q̅) == length(p̅)

    # copy x to q
    q̅ .= x

    # compute p̅ = ϑ(q)
    params.equ.ϑ(params.t + params.Δt, q̅, p̅)

end

"Variational partitioned Runge-Kutta integrator."
mutable struct IntegratorVPRKdegenerate{DT, TT,
                SPT <: ParametersVPRK{DT,TT},
                PPT <: ParametersVPRKdegenerate{DT,TT},
                SST <: NonlinearSolver{DT},
                STP <: NonlinearSolver{DT},
                IT  <: InitialGuessPODE{DT,TT}} <: AbstractIntegratorVPRK{DT,TT}

    sparams::SPT
    pparams::PPT

    solver::SST
    projector::STP
    iguess::IT

    cache::NonlinearFunctionCacheVPRK{DT}

    q::Vector{Vector{TwicePrecision{DT}}}
    p::Vector{Vector{TwicePrecision{DT}}}
end

function IntegratorVPRKdegenerate(equation::ET, tableau::TableauVPRK{TT}, Δt::TT) where {DT, TT, ET <: IODE{DT,TT}}
    D = equation.d
    M = equation.n
    S = tableau.s

    # create cache for internal stage vectors and update vectors
    cache = NonlinearFunctionCacheVPRK{DT}(D,S)

    # create solver params
    sparams = ParametersVPRK(equation, tableau, Δt)

    # create projector params
    pparams = ParametersVPRKdegenerate(equation, tableau, Δt, cache)

    # create nonlinear solver
    solver = create_nonlinear_solver(DT, D*S, sparams)

    # create projector
    projector = create_nonlinear_solver(DT, D, pparams)

    # create initial guess
    iguess = InitialGuessPODE(get_config(:ig_interpolation), equation, Δt)

    # create solution vectors
    q = create_solution_vector(DT, D, M)
    p = create_solution_vector(DT, D, M)

    # create integrator
    IntegratorVPRKdegenerate{DT, TT, typeof(sparams), typeof(pparams), typeof(solver), typeof(projector), typeof(iguess)}(
                sparams, pparams, solver, projector, iguess, cache, q, p)
end


function initial_guess!(int::IntegratorVPRKdegenerate, m::Int)
    for i in 1:int.sparams.tab.s
        evaluate!(int.iguess, m, int.cache.y, int.cache.z, int.cache.v, int.sparams.tab.q.c[i], int.sparams.tab.p.c[i])
        for k in 1:int.sparams.equ.d
            int.solver.x[int.sparams.equ.d*(i-1)+k] = int.cache.v[k]
        end
    end
end

"Integrate ODE with variational partitioned Runge-Kutta integrator."
function integrate_step!(int::IntegratorVPRKdegenerate{DT,TT}, sol::SolutionPDAE{DT,TT}, m::Int, n::Int) where {DT,TT}
        # check if m and n are compatible with solution dimensions
        check_solution_dimension_asserts(sol, m, n)

        # set time and solution for nonlinear solver
        int.sparams.t = sol.t[0] + (n-1)*int.sparams.Δt
        int.sparams.q .= int.q[m]
        int.sparams.p .= int.p[m]

        # set time and solution for solution solver
        int.pparams.t = sol.t[0] + (n-1)*int.pparams.Δt
        int.pparams.q .= int.q[m]
        int.pparams.p .= int.p[m]

        # compute initial guess
        initial_guess!(int, m)

        # call nonlinear solver
        solve!(int.solver)

        # print solver status
        print_solver_status(int.solver.status, int.solver.params, n)

        # check if solution contains NaNs or error bounds are violated
        check_solver_status(int.solver.status, int.solver.params, n)

        # compute vector fields at internal stages
        compute_stages!(int.solver.x, int.cache.Q, int.cache.V, int.cache.P, int.cache.F, int.sparams)

        # compute internal stage solution
        update_solution!(int.q[m], int.cache.V, int.sparams.tab.q.b, int.sparams.tab.q.b̂, int.sparams.Δt)
        update_solution!(int.p[m], int.cache.F, int.sparams.tab.p.b, int.sparams.tab.p.b̂, int.sparams.Δt)

        # set initial guess for solution
        int.projector.x .= int.q[m]

        # call solution solver
        solve!(int.projector)

        # print solver status
        print_solver_status(int.projector.status, int.projector.params, n)

        # check if solution contains NaNs or error bounds are violated
        check_solver_status(int.projector.status, int.projector.params, n)

        # copy solution to integrator
        compute_solution!(int.projector.x, int.cache.q̅, int.cache.p̅, int.pparams)
        int.q[m] .= int.cache.q̅
        int.p[m] .= int.cache.p̅

        # copy solution to initial guess for next time step
        update!(int.iguess, m, sol.t[0] + n*int.sparams.Δt, int.q[m], int.p[m])

        # take care of periodic solutions
        cut_periodic_solution!(int.q[m], int.sparams.equ.periodicity)

        # copy to solution
        copy_solution!(sol, int.q[m], int.p[m], n, m)
end
