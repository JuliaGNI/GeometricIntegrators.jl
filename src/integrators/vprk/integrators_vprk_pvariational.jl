
"Parameters for right-hand side function of variational partitioned Runge-Kutta methods."
type NonlinearFunctionParametersVPRKpVariational{DT, TT, ET <: IODE{DT,TT}, D, S} <: AbstractNonlinearFunctionParametersVPRK{DT,TT,ET,D,S}
    equ::ET
    tab::TableauVPRK{TT}
    Δt::TT

    R::Vector{TT}
    R1::Vector{TT}
    R2::Vector{TT}

    t::TT
    q::Vector{DT}
    p::Vector{DT}
end

function NonlinearFunctionParametersVPRKpVariational(equ::ET, tab::TableauVPRK{TT}, Δt::TT, R::Vector) where {DT, TT, ET <: IODE{DT,TT}}
    R  = convert(Vector{TT}, R)
    R1 = [one(TT), zero(TT)]
    R2 = [zero(TT), one(TT)]

    q = zeros(DT, equ.d)
    p = zeros(DT, equ.d)

    NonlinearFunctionParametersVPRKpVariational{DT, TT, ET, equ.d, tab.s}(equ, tab, Δt, R, R1, R2, 0, q, p)
end


@generated function compute_projection!(
                x::Vector{ST}, q̅::Vector{ST}, p̅::Vector{ST},
                λ::Vector{ST}, U::Matrix{ST}, G::Matrix{ST},
                params::NonlinearFunctionParametersVPRKpVariational{DT,TT,ET,D,S}
            ) where {ST,DT,TT,ET,D,S}

    # create temporary vectors
    tG = zeros(ST,D)

    quote
        @assert length(q̅) == length(p̅) == length(λ) == size(U,1) == size(G,1)
        @assert length(x) ≥ 2length(q̅)

        # copy x to q, λ
        for k in 1:D
            q̅[k] = x[0*D+k]
            λ[k] = x[1*D+k]
        end

        # compute u=λ and g=∇α(q)⋅λ
        # simd_copy_yx_first!(λ, U, 1)
        simd_copy_yx_first!(λ, U, 2)

        params.equ.g(params.t + params.Δt, q̅, λ, $tG)
        simd_copy_yx_first!($tG, G, 1)
        # simd_copy_yx_first!($tG, G, 2)

        # compute p̅=α(q)
        params.equ.α(params.t + params.Δt, q̅, λ, p̅)
    end
end

"Compute stages of projected variational partitioned Runge-Kutta methods."
@generated function function_stages!(x::Vector{ST}, b::Vector{ST},
                params::NonlinearFunctionParametersVPRKpVariational{DT,TT,ET,D,S}
            ) where {ST,DT,TT,ET,D,S}

    cache = NonlinearFunctionCacheVPRKprojection{ST}(D)

    function_stages = quote
        @assert length(x) == length(b)

        compute_projection!(x, $cache.q̅, $cache.p̅, $cache.λ, $cache.U, $cache.G, params)

        # # compute b = - [q̅-q-U]
        for k in 1:D
            b[0*D+k] = - ($cache.q̅[k] - params.q[k]) + params.Δt * params.R[2] * $cache.U[k,2]
        end

        # compute b = - [p̅-p-G]
        for k in 1:D
            b[1*D+k] = - ($cache.p̅[k] - params.p[k])
        end
    end

    return function_stages
end

"Variational partitioned Runge-Kutta integrator."
immutable IntegratorVPRKpVariational{DT, TT,
                SPT <: NonlinearFunctionParametersVPRK{DT,TT},
                PPT <: NonlinearFunctionParametersVPRKpVariational{DT,TT},
                SST <: NonlinearSolver{DT},
                STP <: NonlinearSolver{DT},
                IT <: InitialGuessPODE{DT,TT}} <: AbstractIntegratorVPRK{DT,TT}

    sparams::SPT
    pparams::PPT

    solver::SST
    projector::STP
    iguess::IT

    scache::NonlinearFunctionCacheVPRK{DT}
    pcache::NonlinearFunctionCacheVPRKprojection{DT}

    q::Vector{Vector{Double{DT}}}
    p::Vector{Vector{Double{DT}}}
end

function IntegratorVPRKpVariational(equation::ET, tableau::TableauVPRK{TT}, Δt::TT;
                                    R=[1,1]) where {DT, TT, ET <: IODE{DT,TT}}

    D = equation.d
    S = tableau.s

    # create solver params
    sparams = NonlinearFunctionParametersVPRK(equation, tableau, Δt)

    # create projector params
    pparams = NonlinearFunctionParametersVPRKpVariational(equation, tableau, Δt)

    # create nonlinear solver
    solver = create_nonlinear_solver(DT, D*S, sparams)

    # create projector
    projector = create_nonlinear_solver(DT, 2*D, pparams)

    # create initial guess
    iguess = InitialGuessPODE(get_config(:ig_interpolation), equation, Δt)

    # create cache for internal stage vectors and update vectors
    scache = NonlinearFunctionCacheVPRK{DT}(D,S)
    pcache = NonlinearFunctionCacheVPRKprojection{DT}(D)

    # create solution vectors
    q = create_solution_vector_double_double(DT, D, M)
    p = create_solution_vector_double_double(DT, D, M)

    # create integrator
    IntegratorVPRKpVariational{DT, TT, ET, typeof(sparams), typeof(pparams), typeof(solver), typeof(projector), typeof(iguess)}(
            sparams, pparams, solver, projector, iguess, scache, pcache, q, p)
end

equation(integrator::IntegratorVPRKpVariational) = integrator.sparams.equ
timestep(integrator::IntegratorVPRKpVariational) = integrator.sparams.Δt
tableau(integrator::IntegratorVPRKpVariational) = integrator.sparams.tab
dims(integrator::IntegratorVPRKpVariational) = integrator.sparams.equ.d


function initialize!(int::IntegratorVPRKpVariational{DT,TT}, sol::SolutionPDAE, m::Int) where {DT,TT}
    @assert m ≥ 1
    @assert m ≤ sol.ni

    local tq = zeros(DT, int.equation.d)
    local tG = zeros(DT, int.equation.d)

    # copy initial conditions from solution
    get_initial_conditions!(sol, int.q[m], int.p[m], int.pcache.λ, m)

    # initialise initial guess
    initialize!(int.iguess, m, sol.t[0], int.q[m], int.p[m])

    # initialise projector
    tq .= int.q[m]
    simd_copy_yx_first!(int.pcache.λ, int.pcache.U, 1)
    int.equation.g(sol.t[0], tq, int.pcache.λ, tG)
    simd_copy_yx_first!(tG, int.pcache.G, 1)

    # add perturbation for first time step to solution
    update_solution!(int.q[m], int.pcache.U, int.pparams.R1, int.Δt)
    update_solution!(int.p[m], int.pcache.G, int.pparams.R1, int.Δt)
end


function initial_guess!(int::IntegratorVPRKpVariational{DT,TT}, m::Int) where {DT,TT}
    for i in 1:params.tab.s
        evaluate!(int.iguess, int.scache.y, int.scache.z, int.scache.v, params.tab.q.c[i], params.tab.p.c[i])
        for k in 1:int.params.equ.d
            int.solver.x[int.params.equ.d*(i-1)+k] = int.scache.v[k]
        end
    end
end


"Integrate ODE with variational partitioned Runge-Kutta integrator."
function integrate!(int::IntegratorVPRKpVariational{DT,TT}, sol::SolutionPDAE{DT,TT}, m::Int, n::Int) where {DT,TT}
    # check if m and n are compatible with solution dimensions
    check_solution_dimension_asserts(sol, m, n)

    # set time for nonlinear and projection solver
    int.sparams.t = sol.t[0] + (n-1)*int.params.Δt
    int.sparams.q .= int.q[m]
    int.sparams.p .= int.p[m]

    int.pparams.t = sol.t[0] + (n-1)*int.params.Δt
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

    # compute unprojected solution
    compute_stages!(int.solver.x, int.scache.Q, int.scache.V, int.scache.P, int.scache.F, int.sparams)

    # compute intermediate update
    update_solution!(int.q[m], int.scache.V, int.tableau.q.b, int.tableau.q.b̂, int.Δt)
    update_solution!(int.p[m], int.scache.F, int.tableau.p.b, int.tableau.p.b̂, int.Δt)

    # set initial guess for projection
    for k in 1:int.params.equ.d
        int.projector.x[0*int.params.equ.d+k] = int.q[k]
        int.projector.x[1*int.params.equ.d+k] = 0
    end

    # call projection solver
    solve!(int.projector)

    # print projector status
    print_solver_status(int.projector.status, int.projector.params, n)

    # check if solution contains NaNs or error bounds are violated
    check_solver_status(int.projector.status, int.projector.params, n)

    # compute projection vector fields
    compute_projection!(int.projector.x, int.pcache.q̅, int.pcache.p̅, int.pcache.λ, int.pcache.U, int.pcache.G, int.pparams)

    # compute final update
    update_solution!(int.q[m], int.pcache.U, int.pparams.R2, int.Δt)
    update_solution!(int.p[m], int.pcache.G, int.pparams.R2, int.Δt)

    # copy solution to initial guess for next time step
    update!(int.iguess, m, sol.t[0] + n*int.params.Δt, int.q[m], int.p[m])

    # take care of periodic solutions
    cut_periodic_solution!(int.q[m], int.params.equ.periodicity)

    # copy to solution
    copy_solution!(sol, int.q[m], int.p[m], int.pcache.λ, n, m)

    # add perturbation to solution
    # (same vector field as previous time step)
    # project_solution!(int, int.pcache, int.pparams.R1)
    update_solution!(int.q[m], int.pcache.U, int.pparams.R1, int.Δt)
    update_solution!(int.p[m], int.pcache.G, int.pparams.R1, int.Δt)
end
