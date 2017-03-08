
"Parameters for right-hand side function of variational partitioned Runge-Kutta methods."
type NonlinearFunctionParametersVPRKpStandard{DT,TT,AT,GT,D} <: NonlinearFunctionParameters{DT}
    α::AT
    g::GT

    Δt::TT

    R::Vector{TT}
    R1::Vector{TT}
    R2::Vector{TT}

    t::TT

    q::Vector{DT}
    p::Vector{DT}
end

function NonlinearFunctionParametersVPRKpStandard{DT,TT,AT,GT}(α::AT, g::GT, Δt::TT, R::Vector, q::Vector{DT}, p::Vector{DT})
    @assert length(q)==length(p)
    R = convert(Vector{TT}, R)
    R1 = [R[1], zero(TT)]
    R2 = [zero(TT), R[2]]
    NonlinearFunctionParametersVPRKpStandard{DT,TT,AT,GT,length(q)}(α, g, Δt, R, R1, R2, 0, q, p)
end


@generated function compute_projection_vprk!{ST,DT,TT,ΑT,GT,D}(x::Vector{ST}, q̅::Vector{ST}, p̅::Vector{ST}, λ::Vector{ST}, U::Matrix{ST}, G::Matrix{ST}, params::NonlinearFunctionParametersVPRKpStandard{DT,TT,ΑT,GT,D})
    # create temporary vectors
    tG = zeros(ST,D)

    compute_projection_vprk = quote
        @assert length(q̅) == length(p̅) == length(λ) == size(U,1) == size(G,1)
        @assert length(x) ≥ 2length(q̅)

        # copy x to q, λ
        for k in 1:D
            q̅[k] = x[0*D+k]
            λ[k] = x[1*D+k]
        end

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
@generated function function_stages!{ST,DT,TT,ΑT,GT,D}(x::Vector{ST}, b::Vector{ST}, params::NonlinearFunctionParametersVPRKpStandard{DT,TT,ΑT,GT,D})
    cache = NonlinearFunctionCacheVPRKprojection{ST}(D)

    function_stages = quote
        @assert length(x) == length(b)

        compute_projection_vprk!(x, $cache.q̅, $cache.p̅, $cache.λ, $cache.U, $cache.G, params)

        # compute b = - [q̅-q-U]
        for k in 1:D
            b[0*D+k] = - ($cache.q̅[k] - params.q[k]) + params.Δt * params.R[2] * $cache.U[k,2]
        end

        # compute b = - [p̅-p-G]
        for k in 1:D
            b[1*D+k] = - ($cache.p̅[k] - params.p[k]) + params.Δt * params.R[2] * $cache.G[k,2]
        end
    end

    return function_stages
end

"Variational partitioned Runge-Kutta integrator."
immutable IntegratorVPRKpStandard{DT, TT, ΑT, FT, GT, VT, SPT, PPT, SST, STP, IT} <: AbstractIntegratorVPRK{DT, TT}
    equation::IODE{DT,TT,ΑT,FT,GT,VT}
    tableau::TableauVPRK{TT}
    Δt::TT

    sparams::SPT
    pparams::PPT

    solver::SST
    projector::STP

    scache::NonlinearFunctionCacheVPRK{DT}
    pcache::NonlinearFunctionCacheVPRKprojection{DT}

    iguess::InitialGuessIODE{DT, TT, VT, FT, IT}

    q::Vector{DT}
    p::Vector{DT}

    qₑᵣᵣ::Vector{DT}
    pₑᵣᵣ::Vector{DT}
end

function IntegratorVPRKpSymplectic(args...; kwargs...)
    IntegratorVPRKpStandard(args...; kwargs..., R=[1,1])
end

function IntegratorVPRKpStandard{DT,TT,ΑT,FT,GT,VT}(equation::IODE{DT,TT,ΑT,FT,GT,VT}, tableau::TableauVPRK{TT}, Δt::TT;
                                        R=[0,1],
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
    x̃ = zeros(DT,N)

    # create solution vector for projector
    x = zeros(DT,2*D)

    # create solution vectors
    q = zeros(DT,D)
    p = zeros(DT,D)

    # create compensated summation error vectors
    qₑᵣᵣ = zeros(DT,D)
    pₑᵣᵣ = zeros(DT,D)

    # create cache for internal stage vectors and update vectors
    scache = NonlinearFunctionCacheVPRK{DT}(D,S)
    pcache = NonlinearFunctionCacheVPRKprojection{DT}(D)

    # create solver params
    sparams = NonlinearFunctionParametersVPRK(equation.α, equation.f, Δt,
                                              tableau.q, tableau.p, d_v,
                                              q, p)

    # create rhs function for nonlinear solver
    function_stages_solver = (x,b) -> function_stages!(x, b, sparams)

    # create nonlinear solver
    solver = nonlinear_solver(x̃, function_stages_solver; nmax=nmax, atol=atol, rtol=rtol, stol=stol)

    # create projector params
    R[2] *= tableau.R∞
    pparams = NonlinearFunctionParametersVPRKpStandard(equation.α, equation.g, Δt, R, q, p)

    # create rhs function for projector
    function_stages_projector = (x,b) -> function_stages!(x, b, pparams)

    # create projector
    projector = nonlinear_solver(x, function_stages_projector; nmax=nmax, atol=atol, rtol=rtol, stol=stol)

    # create initial guess
    iguess = InitialGuessIODE(interpolation, equation, Δt)


    IntegratorVPRKpStandard{DT, TT, ΑT, FT, GT, VT, typeof(sparams), typeof(pparams), typeof(solver), typeof(projector), typeof(iguess.int)}(
                                        equation, tableau, Δt, sparams, pparams, solver, projector, scache, pcache, iguess,
                                        q, p, qₑᵣᵣ, pₑᵣᵣ)
end


"Integrate ODE with variational partitioned Runge-Kutta integrator."
function integrate!{DT,TT,ΑT,FT,GT,VT,N}(int::IntegratorVPRKpStandard{DT,TT,ΑT,FT,GT,VT}, sol::SolutionPDAE{DT,TT,N}, m1::Int, m2::Int)
    @assert m1 ≥ 1
    @assert m2 ≤ sol.ni

    # loop over initial conditions
    for m in m1:m2
        # copy initial conditions from solution
        get_initial_conditions!(sol, int.q, int.p, m)

        # initialise initial guess
        initialize!(int.iguess, sol.t[0], int.q, int.p)

        for n in 1:sol.ntime
            # set time for nonlinear and projection solver
            int.sparams.t = sol.t[0] + (n-1)*int.Δt
            int.pparams.t = sol.t[0] + (n-1)*int.Δt

            # add perturbation to solution (same vector field as previous time step)
            project_solution!(int, int.pcache, int.pparams.R1)

            # copy perturbed previous solution to initial guess
            update!(int.iguess, sol.t[0] + n*int.Δt, int.q, int.p)

            # compute initial guess
            for i in 1:int.tableau.s
                evaluate!(int.iguess, int.scache.y, int.scache.z, int.scache.v, int.tableau.q.c[i], int.tableau.p.c[i])
                for k in 1:int.equation.d
                    int.solver.x[int.equation.d*(i-1)+k] = int.scache.v[k]
                end
            end

            # call nonlinear solver
            solve!(int.solver)

            if !solverStatusOK(int.solver.status, int.solver.params)
                println(int.solver.status, ", sit=", n)
            end

            if isnan(int.solver.status.rₐ)
                break
            end

            # compute unprojected solution
            compute_stages_vprk!(int.solver.x, int.scache.Q, int.scache.V, int.scache.P, int.scache.F, int.sparams)
            update_solution!(int, int.scache)

            # set initial guess for projection
            for k in 1:int.equation.d
                int.projector.x[0*int.equation.d+k] = int.q[k]
                int.projector.x[1*int.equation.d+k] = 0
            end

            # call projection solver
            solve!(int.projector)

            if !solverStatusOK(int.projector.status, int.projector.params)
                println(int.projector.status, ", pit=", n)
            end

            if isnan(int.projector.status.rₐ)
                break
            end

            # add projection to solution
            compute_projection_vprk!(int.projector.x, int.pcache.q̅, int.pcache.p̅, int.pcache.λ, int.pcache.U, int.pcache.G, int.pparams)
            project_solution!(int, int.pcache, int.pparams.R2)

            # take care of periodic solutions
            cut_periodic_solution!(int)

            # copy to solution
            copy_solution!(sol, int.q, int.p, int.pcache.λ, n, m)
        end
    end
end
