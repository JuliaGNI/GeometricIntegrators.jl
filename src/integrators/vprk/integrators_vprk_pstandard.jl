
"Parameters for right-hand side function of variational partitioned Runge-Kutta methods."
mutable struct NonlinearFunctionParametersVPRKpStandard{DT,TT,AT,GT,D,S} <: NonlinearFunctionParameters{DT}
    α::AT
    g::GT

    Δt::TT

    o::Int
    R::Vector{TT}
    R1::Vector{TT}
    R2::Vector{TT}

    t::TT

    q::Vector{DT}
    p::Vector{DT}
end

function NonlinearFunctionParametersVPRKpStandard(DT, α::AT, g::GT, Δt::TT, o::Int, R::Vector, D::Int, S::Int) where {TT,AT,GT}
    R = convert(Vector{TT}, R)
    R1 = [R[1], zero(TT)]
    R2 = [zero(TT), R[2]]
    # R  = [one(TT), one(TT)]
    # R1 = [one(TT), zero(TT)]
    # R2 = [zero(TT), one(TT)]
    # println("R  = ", R)
    # println("R1 = ", R1)
    # println("R2 = ", R2)
    NonlinearFunctionParametersVPRKpStandard{DT,TT,AT,GT,D,S}(α, g, Δt, o, R, R1, R2, 0, zeros(DT,D), zeros(DT,D))
end


@generated function compute_projection!(x::Vector{ST}, q̅::Vector{ST}, p̅::Vector{ST}, λ::Vector{ST}, U::Matrix{ST}, G::Matrix{ST}, params::NonlinearFunctionParametersVPRKpStandard{DT,TT,ΑT,GT,D,S}) where {ST,DT,TT,ΑT,GT,D,S}
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

        params.g(params.t + params.Δt, q̅, λ, $tG)
        simd_copy_yx_first!($tG, G, 1)
        simd_copy_yx_first!($tG, G, 2)

        # compute p̅=α(q)
        params.α(params.t + params.Δt, q̅, λ, p̅)
    end

    return compute_projection_vprk
end

"Compute stages of projected variational partitioned Runge-Kutta methods."
@generated function function_stages!(x::Vector{ST}, b::Vector{ST}, params::NonlinearFunctionParametersVPRKpStandard{DT,TT,ΑT,GT,D,S}) where {ST,DT,TT,ΑT,GT,D,S}
    cache = NonlinearFunctionCacheVPRKprojection{ST}(D,S)

    function_stages = quote
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

    return function_stages
end

"Variational partitioned Runge-Kutta integrator."
struct IntegratorVPRKpStandard{DT, TT, ΑT, FT, GT, VT, SPT, PPT, SST, STP, IT} <: AbstractIntegratorVPRK{DT, TT}
    equation::IODE{DT,TT,ΑT,FT,GT,VT}
    tableau::TableauVPRK{TT}
    Δt::TT

    sparams::SPT
    pparams::PPT

    solver::SST
    projector::STP

    scache::NonlinearFunctionCacheVPRK{DT}
    pcache::NonlinearFunctionCacheVPRKprojection{DT}

    iguess::InitialGuessPODE{DT, TT, VT, FT, IT}

    q::Vector{Vector{Double{DT}}}
    p::Vector{Vector{Double{DT}}}
    # q::Vector{Vector{DT}}
    # p::Vector{Vector{DT}}
end

function IntegratorVPRKpSymplectic(args...; kwargs...)
    IntegratorVPRKpStandard(args...; kwargs..., R=[1,1])
end

function IntegratorVPRKpStandard(equation::IODE{DT,TT,ΑT,FT,GT,VT}, tableau::TableauVPRK{TT}, Δt::TT;
                                 R=[0,1], interpolation=HermiteInterpolation{DT}) where {DT,TT,ΑT,FT,GT,VT}
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
    x̃ = zeros(DT,N)

    # create solution vector for projector
    x = zeros(DT,2*D)
    # x = zeros(DT,1*D)

    # create solution vectors
    q = Array{Vector{Double{DT}}}(M)
    p = Array{Vector{Double{DT}}}(M)
    # q = Array{Vector{DT}}(M)
    # p = Array{Vector{DT}}(M)

    for i in 1:M
        q[i] = zeros(Double{DT},D)
        p[i] = zeros(Double{DT},D)
        # q[i] = zeros(DT,D)
        # p[i] = zeros(DT,D)
    end

    # create cache for internal stage vectors and update vectors
    scache = NonlinearFunctionCacheVPRK{DT}(D,S)
    pcache = NonlinearFunctionCacheVPRKprojection{DT}(D,S)

    # create solver params
    sparams = NonlinearFunctionParametersVPRK(DT, D, equation.α, equation.f, Δt,
                                              tableau.q, tableau.p, d_v)

    # create rhs function for nonlinear solver
    function_stages_solver = (x,b) -> function_stages!(x, b, sparams)

    # create nonlinear solver
    solver = get_config(:nls_solver)(x̃, function_stages_solver)

    # create projector params
    R[2] *= tableau.R∞
    pparams = NonlinearFunctionParametersVPRKpStandard(DT, equation.α, equation.g, Δt, tableau.o, R, D, S)

    # create rhs function for projector
    function_stages_projector = (x,b) -> function_stages!(x, b, pparams)

    # create projector
    projector = get_config(:nls_solver)(x, function_stages_projector)

    # create initial guess
    iguess = InitialGuessPODE(interpolation, equation, Δt)


    IntegratorVPRKpStandard{DT, TT, ΑT, FT, GT, VT, typeof(sparams), typeof(pparams), typeof(solver), typeof(projector), typeof(iguess.int)}(
                                        equation, tableau, Δt, sparams, pparams, solver, projector, scache, pcache, iguess, q, p)
end


function initialize!(int::IntegratorVPRKpStandard{DT,TT}, sol::SolutionPDAE, m::Int) where {DT,TT}
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


function initial_guess!(int::IntegratorVPRKpStandard, m::Int)
    for i in 1:int.tableau.s
        evaluate!(int.iguess, m, int.scache.y, int.scache.z, int.scache.v, int.tableau.q.c[i], int.tableau.p.c[i])
        for k in 1:int.equation.d
            int.solver.x[int.equation.d*(i-1)+k] = int.scache.v[k]
        end
    end
end

"Integrate ODE with variational partitioned Runge-Kutta integrator."
function integrate_step!(int::IntegratorVPRKpStandard{DT,TT,ΑT,FT,GT,VT}, sol::SolutionPDAE{DT,TT,N}, m::Int, n::Int) where {DT,TT,ΑT,FT,GT,VT,N}
    @assert m ≥ 1
    @assert m ≤ sol.ni

    @assert n ≥ 1
    @assert n ≤ sol.ntime

    # set time and solution for nonlinear solver
    int.sparams.t = sol.t[0] + (n-1)*int.Δt
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
    update_solution!(int.q[m], int.scache.V, int.tableau.q.b, int.tableau.q.b̂, int.Δt)
    update_solution!(int.p[m], int.scache.F, int.tableau.p.b, int.tableau.p.b̂, int.Δt)

    # set time and solution for projection solver
    int.pparams.t = sol.t[0] + (n-1)*int.Δt
    int.pparams.q .= int.q[m]
    int.pparams.p .= int.p[m]

    # set initial guess for projection
    for k in 1:int.equation.d
        int.projector.x[0*int.equation.d+k] = int.pparams.q[k]
        int.projector.x[1*int.equation.d+k] = 0
        # int.projector.x[0*int.equation.d+k] = 0
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
    update_solution!(int.q[m], int.pcache.U, int.pparams.R2, int.Δt)
    update_solution!(int.p[m], int.pcache.G, int.pparams.R2, int.Δt)

    # copy solution to initial guess for next time step
    update!(int.iguess, m, sol.t[0] + n*int.Δt, int.q[m], int.p[m])

    # take care of periodic solutions
    cut_periodic_solution!(int.q[m], int.equation.periodicity)

    # copy to solution
    copy_solution!(sol, int.q[m], int.p[m], int.pcache.λ, n, m)

    # add perturbation for next time step to solution
    # (same vector field as previous time step)
    update_solution!(int.q[m], int.pcache.U, int.pparams.R1, int.Δt)
    update_solution!(int.p[m], int.pcache.G, int.pparams.R1, int.Δt)
end
