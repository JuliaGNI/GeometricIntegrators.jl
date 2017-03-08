
"Parameters for right-hand side function of variational partitioned Runge-Kutta methods."
type NonlinearFunctionParametersVPRKpSymmetric{DT,TT,ΑT,FT,GT,D,S} <: AbstractNonlinearFunctionParametersVPRK{DT,TT,ΑT,FT,D,S}
    α::ΑT
    f::FT
    g::GT

    Δt::TT

    t_q::CoefficientsRK{TT}
    t_p::CoefficientsRK{TT}
    d_v::Vector{TT}
    R::Vector{TT}

    t::TT

    q::Vector{DT}
    p::Vector{DT}

    function NonlinearFunctionParametersVPRKpSymmetric(α, f, g, Δt, t_q, t_p, d_v, R∞, q, p)
        R = convert(Vector{TT}, [1, R∞])
        new(α, f, g, Δt, t_q, t_p, d_v, R, 0, q, p)
    end
end

@generated function compute_projection_vprk!{ST,DT,TT,ΑT,FT,GT,D,S}(x::Vector{ST}, q̅::Vector{ST}, p̅::Vector{ST}, λ::Vector{ST}, V::Matrix{ST}, U::Matrix{ST}, G::Matrix{ST}, params::NonlinearFunctionParametersVPRKpSymmetric{DT,TT,ΑT,FT,GT,D,S})
    # create temporary vectors
    tG = zeros(ST,D)

    compute_projection_vprk = quote
        local t₀::TT = params.t
        local t₁::TT = params.t + params.Δt

        # copy x to λ and q̅
        for k in 1:D
            q̅[k] = x[D*(S+0)+k]
            λ[k] = x[D*(S+1)+k]
        end

        # compute U=λ
        simd_copy_yx_first!(λ, U, 1)
        simd_copy_yx_first!(λ, U, 2)

        # compute G=g(q,λ)
        params.g(t₀, params.q, λ, $tG)
        simd_copy_yx_first!($tG, G, 1)

        params.g(t₁, q̅, λ, $tG)
        simd_copy_yx_first!($tG, G, 2)

        # compute p̅=α(q̅)p̅
        params.α(t₁, q̅, λ, p̅)
    end

    return compute_projection_vprk
end


"Compute stages of variational partitioned Runge-Kutta methods."
@generated function function_stages!{ST,DT,TT,ΑT,FT,GT,D,S}(x::Vector{ST}, b::Vector{ST}, params::NonlinearFunctionParametersVPRKpSymmetric{DT,TT,ΑT,FT,GT,D,S})
    scache = NonlinearFunctionCacheVPRK{ST}(D,S)
    pcache = NonlinearFunctionCacheVPRKprojection{ST}(D)

    function_stages = quote
        compute_stages_vprk!(x, $pcache.q̅, $pcache.p̅, $pcache.λ, $scache.Q, $scache.V, $pcache.U, $scache.P, $scache.F, $pcache.G, params)

        # compute b = - [P-AF-U]
        compute_rhs_vprk!(b, $scache.P, $scache.F, $pcache.G, params)

        # compute b = - [p-bF-G]
        compute_rhs_vprk_projection_p!(b, $pcache.p̅, $scache.F, $pcache.G, D*(S+0), params)

        # compute b = - [q-bV-U]
        compute_rhs_vprk_projection_q!(b, $pcache.q̅, $scache.V, $pcache.U, D*(S+1), params)

        compute_rhs_vprk_correction!(b, $scache.V, params)
    end

    return function_stages
end


"Variational partitioned Runge-Kutta integrator."
immutable IntegratorVPRKpSymmetric{DT,TT,ΑT,FT,GT,VT,FPT,ST,IT} <: Integrator{DT,TT}
    equation::IODE{DT,TT,ΑT,FT,GT,VT}
    tableau::TableauVPRK{TT}
    Δt::TT

    params::FPT
    solver::ST

    scache::NonlinearFunctionCacheVPRK{DT}
    pcache::NonlinearFunctionCacheVPRKprojection{DT}

    iguess::InitialGuessIODE{DT,TT,VT,FT,IT}

    q::Array{DT,1}
    p::Array{DT,1}

    qₑᵣᵣ::Vector{DT}
    pₑᵣᵣ::Vector{DT}
end

function IntegratorVPRKpSymmetric{DT,TT,ΑT,FT,GT,VT}(equation::IODE{DT,TT,ΑT,FT,GT,VT}, tableau::TableauVPRK{TT}, Δt::TT;
                                        nonlinear_solver=DEFAULT_NonlinearSolver,
                                        nmax=DEFAULT_nmax, atol=DEFAULT_atol, rtol=DEFAULT_rtol, stol=DEFAULT_stol,
                                        interpolation=HermiteInterpolation{DT})
    D = equation.d
    S = tableau.s

    N = D*(S+2)

    if isdefined(tableau, :d)
        tableau_d = tableau.d
    else
        tableau_d = DT[]
    end

    # create solution vectors
    q = zeros(DT,D)
    p = zeros(DT,D)

    # create compensated summation error vectors
    qₑᵣᵣ = zeros(DT,D)
    pₑᵣᵣ = zeros(DT,D)

    # create cache for internal stage vectors and update vectors
    scache = NonlinearFunctionCacheVPRK{DT}(D,S)
    pcache = NonlinearFunctionCacheVPRKprojection{DT}(D)

    # create params
    params = NonlinearFunctionParametersVPRKpSymmetric{DT,TT,ΑT,FT,GT,D,S}(
                                                equation.α, equation.f, equation.g, Δt,
                                                tableau.q, tableau.p, tableau_d, tableau.R∞,
                                                q, p)

    # create rhs function for nonlinear solver
    function_stages_solver = (x,b) -> function_stages!(x, b, params)

    # create solver
    solver = nonlinear_solver(zeros(DT,N), function_stages_solver; nmax=nmax, atol=atol, rtol=rtol, stol=stol)

    # create initial guess
    iguess = InitialGuessIODE(interpolation, equation, Δt)

    IntegratorVPRKpSymmetric{DT, TT, ΑT, FT, GT, VT, typeof(params), typeof(solver), typeof(iguess.int)}(
                                        equation, tableau, Δt, params, solver, scache, pcache, iguess,
                                        q, p, qₑᵣᵣ, pₑᵣᵣ)
end


"Integrate ODE with variational partitioned Runge-Kutta integrator."
function integrate!{DT,TT,ΑT,FT,GT,VT,N}(int::IntegratorVPRKpSymmetric{DT,TT,ΑT,FT,GT,VT}, sol::SolutionPDAE{DT,TT,N})
    # loop over initial conditions
    for m in 1:sol.ni
        local j::Int
        local tqᵢ::TT
        local tpᵢ::TT

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
                evaluate!(int.iguess, int.scache.y, int.scache.z, int.scache.v, int.tableau.q.c[i], int.tableau.p.c[i])
                for k in 1:int.equation.d
                    int.solver.x[int.equation.d*(i-1)+k] = int.scache.v[k]
                end
            end
            evaluate!(int.iguess, int.scache.y, int.scache.z, int.scache.v, one(TT), one(TT))
            for k in 1:int.equation.d
                int.solver.x[int.equation.d*(int.tableau.s+0)+k] = int.scache.y[k]
            end
            for k in 1:int.equation.d
                int.solver.x[int.equation.d*(int.tableau.s+1)+k] = 0
            end

            # call nonlinear solver
            solve!(int.solver)

            if !solverStatusOK(int.solver.status, int.solver.params)
                println(int.solver.status, ", it=", n)
            end

            if isnan(int.solver.status.rₐ)
                break
            end

            compute_stages_vprk!(int.solver.x,
                                 int.pcache.q̅, int.pcache.p̅, int.pcache.λ,
                                 int.scache.Q, int.scache.V, int.pcache.U,
                                 int.scache.P, int.scache.F, int.pcache.G, int.params)

            # copy solution
            simd_mult!(int.scache.y, int.scache.V, int.tableau.q.b)
            simd_mult!(int.scache.z, int.scache.F, int.tableau.p.b)
            simd_axpy!(int.Δt, int.scache.y, int.q, int.qₑᵣᵣ)
            simd_axpy!(int.Δt, int.scache.z, int.p, int.pₑᵣᵣ)

            simd_mult!(int.pcache.u, int.pcache.U, int.params.R)
            simd_mult!(int.pcache.g, int.pcache.G, int.params.R)
            simd_axpy!(int.Δt, int.pcache.u, int.q, int.qₑᵣᵣ)
            simd_axpy!(int.Δt, int.pcache.g, int.p, int.pₑᵣᵣ)

            # take care of periodic solutions
            for k in 1:int.equation.d
                if int.equation.periodicity[k] ≠ 0
                    int.q[k] = mod(int.q[k], int.equation.periodicity[k])
                end
            end

            # copy to solution
            copy_solution!(sol, int.q, int.p, int.pcache.λ, n, m)
        end
    end
end
