
"Parameters for right-hand side function of variational partitioned Runge-Kutta methods."
mutable struct ParametersVPRKpStandard{DT, TT, ET <: IODE{DT,TT}, D, S} <: AbstractParametersVPRK{DT,TT,ET,D,S}
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

function ParametersVPRKpStandard(equ::ET, tab::TableauVPRK{TT}, Δt::TT, R::Vector) where {DT, TT, ET <: IODE{DT,TT}}
    R  = Vector{TT}(R)
    R1 = [R[1], zero(TT)]
    R2 = [zero(TT), tab.R∞ * R[2]]

    # TODO make R a matrix with different components for q and p, i.e.,
    #      different components for U1,U2,G1,G2
    #      -> this will account for the cases of variational projections

    q = zeros(DT, equ.d)
    p = zeros(DT, equ.d)

    ParametersVPRKpStandard{DT, TT, ET, equ.d, tab.s}(equ, tab, Δt, R, R1, R2, 0, q, p)
end


@generated function compute_projection!(
                x::Vector{ST}, q̅::Vector{ST}, p̅::Vector{ST},
                λ::Vector{ST}, U::Matrix{ST}, G::Matrix{ST},
                params::ParametersVPRKpStandard{DT,TT,ET,D,S}
            ) where {ST,DT,TT,ET,D,S}

    # create temporary vectors
    tG = zeros(ST,D)

    compute_projection_vprk = quote
        @assert length(q̅) == length(p̅) == length(λ) == size(U,1) == size(G,1)
        # @assert length(x) ≥ 2length(q̅)

        # copy x to q, λ
        for k in 1:D
            q̅[k] = x[0*D+k]
            λ[k] = x[1*D+k]
        end
        # for k in 1:D
        #     λ[k] = x[0*D+k]
        #     q̅[k] = params.q[k] + params.Δt * params.R[2] * λ[k]
        # end

        # compute u=λ and g=∇α(q)⋅λ
        simd_copy_yx_first!(λ, U, 1)
        simd_copy_yx_first!(λ, U, 2)

        params.equ.g(params.t + params.Δt, q̅, λ, $tG)
        simd_copy_yx_first!($tG, G, 1)
        simd_copy_yx_first!($tG, G, 2)

        # compute p̅=α(q)
        params.equ.α(params.t + params.Δt, q̅, λ, p̅)
    end

    return compute_projection_vprk
end

"Compute stages of projected variational partitioned Runge-Kutta methods."
@generated function function_stages!(x::Vector{ST}, b::Vector{ST},
                params::ParametersVPRKpStandard{DT,TT,ET,D,S}
            ) where {ST,DT,TT,ET,D,S}

    cache = NonlinearFunctionCacheVPRKprojection{ST}(D,S)

    quote
        @assert length(x) == length(b)

        compute_projection!(x, $cache.q̅, $cache.p̅, $cache.λ, $cache.U, $cache.G, params)

        # compute b = - [q̅-q-U]
        for k in 1:D
            b[0*D+k] = - ($cache.q̅[k] - params.q[k]) + params.Δt * params.R[2] * $cache.U[k,2]
        end

        # compute b = - [p̅-p-G]
        for k in 1:D
            b[1*D+k] = - ($cache.p̅[k] - params.p[k]) + params.Δt * params.R[2] * $cache.G[k,2]
            # b[0*D+k] = - ($cache.p̅[k] - params.p[k]) + params.Δt * params.R[2] * $cache.G[k,2]
        end
    end
end

"Variational partitioned Runge-Kutta integrator."
mutable struct IntegratorVPRKpStandard{DT, TT,
                SPT <: ParametersVPRK{DT,TT},
                PPT <: ParametersVPRKpStandard{DT,TT},
                SST <: NonlinearSolver{DT},
                STP <: NonlinearSolver{DT},
                IT  <: InitialGuessPODE{DT,TT}} <: AbstractIntegratorVPRK{DT,TT}

    sparams::SPT
    pparams::PPT

    solver::SST
    projector::STP
    iguess::IT

    scache::NonlinearFunctionCacheVPRK{DT}
    pcache::NonlinearFunctionCacheVPRKprojection{DT}

    q::Vector{Vector{TwicePrecision{DT}}}
    p::Vector{Vector{TwicePrecision{DT}}}
end

function IntegratorVPRKpSymplectic(args...; kwargs...)
    IntegratorVPRKpStandard(args...; kwargs..., R=[1,1])
end

# function IntegratorVPRKpVariationalQ(args...; kwargs...)
#     IntegratorVPRKpStandard(args...; kwargs..., R=[...])
# end
#
# function IntegratorVPRKpVariationalP(args...; kwargs...)
#     IntegratorVPRKpStandard(args...; kwargs..., R=[...])
# end

function IntegratorVPRKpStandard(equation::ET, tableau::TableauVPRK{TT}, Δt::TT;
                                 R=[0,1]) where {DT, TT, ET <: IODE{DT,TT}}
    D = equation.d
    M = equation.n
    S = tableau.s

    # create solver params
    sparams = ParametersVPRK(equation, tableau, Δt)

    # create projector params
    pparams = ParametersVPRKpStandard(equation, tableau, Δt, R)

    # create nonlinear solver
    solver = create_nonlinear_solver(DT, D*S, sparams)

    # create projector
    projector = create_nonlinear_solver(DT, 2*D, pparams)

    # create initial guess
    iguess = InitialGuessPODE(get_config(:ig_interpolation), equation, Δt)

    # create cache for internal stage vectors and update vectors
    scache = NonlinearFunctionCacheVPRK{DT}(D,S)
    pcache = NonlinearFunctionCacheVPRKprojection{DT}(D,S)

    # create solution vectors
    q = create_solution_vector_double_double(DT, D, M)
    p = create_solution_vector_double_double(DT, D, M)

    # create integrator
    IntegratorVPRKpStandard{DT, TT, typeof(sparams), typeof(pparams), typeof(solver), typeof(projector), typeof(iguess)}(
                sparams, pparams, solver, projector, iguess, scache, pcache, q, p)
end

equation(integrator::IntegratorVPRKpStandard) = integrator.sparams.equ
timestep(integrator::IntegratorVPRKpStandard) = integrator.sparams.Δt
tableau(integrator::IntegratorVPRKpStandard) = integrator.sparams.tab
dims(integrator::IntegratorVPRKpStandard) = integrator.sparams.equ.d


function initialize!(int::IntegratorVPRKpStandard{DT,TT}, sol::SolutionPDAE, m::Int) where {DT,TT}
    @assert m ≥ 1
    @assert m ≤ sol.ni

    # copy initial conditions from solution
    get_initial_conditions!(sol, int.q[m], int.p[m], int.pcache.λ, m)

    # initialise initial guess
    initialize!(int.iguess, m, sol.t[0], int.q[m], int.p[m])

    local tq = zeros(DT, int.pparams.equ.d)
    local tG = zeros(DT, int.pparams.equ.d)

    # initialise projector
    tq .= int.q[m]
    simd_copy_yx_first!(int.pcache.λ, int.pcache.U, 1)
    int.pparams.equ.g(sol.t[0], tq, int.pcache.λ, tG)
    simd_copy_yx_first!(tG, int.pcache.G, 1)

    # add perturbation for first time step to solution
    update_solution!(int.q[m], int.pcache.U, int.pparams.R1, int.pparams.Δt)
    update_solution!(int.p[m], int.pcache.G, int.pparams.R1, int.pparams.Δt)
end


function initial_guess!(int::IntegratorVPRKpStandard, m::Int)
    for i in 1:int.sparams.tab.s
        evaluate!(int.iguess, m, int.scache.y, int.scache.z, int.scache.v, int.sparams.tab.q.c[i], int.sparams.tab.p.c[i])
        for k in 1:int.sparams.equ.d
            int.solver.x[int.sparams.equ.d*(i-1)+k] = int.scache.v[k]
        end
    end
end

"Integrate ODE with variational partitioned Runge-Kutta integrator."
function integrate_step!(int::IntegratorVPRKpStandard{DT,TT}, sol::SolutionPDAE{DT,TT}, m::Int, n::Int) where {DT,TT}

    # check if m and n are compatible with solution dimensions
    check_solution_dimension_asserts(sol, m, n)

    # set time and solution for nonlinear solver
    int.sparams.t = sol.t[0] + (n-1)*int.sparams.Δt
    int.sparams.q .= int.q[m]
    int.sparams.p .= int.p[m]

    # compute initial guess
    initial_guess!(int, m)

    # call nonlinear solver
    solve!(int.solver)

    # print solver status
    print_solver_status(int.solver.status, int.solver.params, n)

    # check if solution contains NaNs or error bounds are violated
    check_solver_status(int.solver.status, int.solver.params, n)

    # compute vector fields at internal stages
    compute_stages!(int.solver.x, int.scache.Q, int.scache.V, int.scache.P, int.scache.F, int.sparams)

    # compute unprojected solution
    update_solution!(int.q[m], int.scache.V, int.sparams.tab.q.b, int.sparams.tab.q.b̂, int.sparams.Δt)
    update_solution!(int.p[m], int.scache.F, int.sparams.tab.p.b, int.sparams.tab.p.b̂, int.sparams.Δt)

    # set time and solution for projection solver
    int.pparams.t = sol.t[0] + (n-1)*int.pparams.Δt
    int.pparams.q .= int.q[m]
    int.pparams.p .= int.p[m]

    # set initial guess for projection
    for k in 1:int.pparams.equ.d
        int.projector.x[0*int.pparams.equ.d+k] = int.pparams.q[k]
        int.projector.x[1*int.pparams.equ.d+k] = 0
        # int.projector.x[0*int.pparams.equ.d+k] = 0
    end

    # call projection solver
    solve!(int.projector)

    # print solver status
    print_solver_status(int.projector.status, int.projector.params, n)

    # check if solution contains NaNs or error bounds are violated
    check_solver_status(int.projector.status, int.projector.params, n)

    # compute projection vector fields
    compute_projection!(int.projector.x, int.pcache.q̅, int.pcache.p̅, int.pcache.λ, int.pcache.U, int.pcache.G, int.pparams)

    # add projection to solution
    update_solution!(int.q[m], int.pcache.U, int.pparams.R2, int.pparams.Δt)
    update_solution!(int.p[m], int.pcache.G, int.pparams.R2, int.pparams.Δt)

    # copy solution to initial guess for next time step
    update!(int.iguess, m, sol.t[0] + n*int.sparams.Δt, int.q[m], int.p[m])

    # take care of periodic solutions
    cut_periodic_solution!(int.q[m], int.sparams.equ.periodicity)

    # copy to solution
    copy_solution!(sol, int.q[m], int.p[m], int.pcache.λ, n, m)

    # add perturbation for next time step to solution
    # (same vector field as previous time step)
    update_solution!(int.q[m], int.pcache.U, int.pparams.R1, int.pparams.Δt)
    update_solution!(int.p[m], int.pcache.G, int.pparams.R1, int.pparams.Δt)
end
