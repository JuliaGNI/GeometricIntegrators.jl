
"Parameters for right-hand side function of variational partitioned Runge-Kutta methods."
mutable struct NonlinearFunctionParametersVPRKpSymmetric{DT,TT,ΑT,FT,GT,D,S} <: AbstractNonlinearFunctionParametersVPRK{DT,TT,ΑT,FT,D,S}
    α::ΑT
    f::FT
    g::GT

    Δt::TT

    o::Int
    t_q::CoefficientsRK{TT}
    t_p::CoefficientsRK{TT}
    d_v::Vector{TT}
    R::Vector{TT}

    t::TT

    q::Vector{DT}
    p::Vector{DT}

    function NonlinearFunctionParametersVPRKpSymmetric{DT,TT,ΑT,FT,GT,D,S}(α, f, g, Δt, o, t_q, t_p, d_v, R∞) where {DT,TT,ΑT,FT,GT,D,S}
        R = convert(Vector{TT}, [1, R∞])
        new(α, f, g, Δt, o, t_q, t_p, d_v, R, 0, zeros(DT,D), zeros(DT,D))
    end
end


@generated function compute_projection_vprk!(x::Vector{ST}, q̅::Vector{ST}, p̅::Vector{ST}, λ::Vector{ST}, V::Matrix{ST}, U::Matrix{ST}, G::Matrix{ST}, params::NonlinearFunctionParametersVPRKpSymmetric{DT,TT,ΑT,FT,GT,D,S}) where {ST,DT,TT,ΑT,FT,GT,D,S}
    # create temporary vectors
    tG = zeros(ST,D)

    quote
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
end


"Compute stages of variational partitioned Runge-Kutta methods."
@generated function function_stages!(x::Vector{ST}, b::Vector{ST}, params::NonlinearFunctionParametersVPRKpSymmetric{DT,TT,ΑT,FT,GT,D,S}) where {ST,DT,TT,ΑT,FT,GT,D,S}
    scache = NonlinearFunctionCacheVPRK{ST}(D,S)
    pcache = NonlinearFunctionCacheVPRKprojection{ST}(D,S)

    quote
        compute_stages!(x, $pcache.q̅, $pcache.p̅, $pcache.λ, $scache.Q, $scache.V, $pcache.U, $scache.P, $scache.F, $pcache.G, params)

        # compute b = - [P-AF-U]
        compute_rhs_vprk!(b, $scache.P, $scache.F, $pcache.G, params)

        # compute b = - [p-bF-G]
        compute_rhs_vprk_projection_p!(b, $pcache.p̅, $scache.F, $pcache.G, D*(S+0), params)

        # compute b = - [q-bV-U]
        compute_rhs_vprk_projection_q!(b, $pcache.q̅, $scache.V, $pcache.U, D*(S+1), params)

        compute_rhs_vprk_correction!(b, $scache.V, params)
    end
end


"Variational partitioned Runge-Kutta integrator."
struct IntegratorVPRKpSymmetric{DT,TT,ΑT,FT,GT,VT,FPT,ST,IT} <: AbstractIntegratorVPRK{DT,TT}
    equation::IODE{DT,TT,ΑT,FT,GT,VT}
    tableau::TableauVPRK{TT}
    Δt::TT

    params::FPT
    solver::ST

    scache::NonlinearFunctionCacheVPRK{DT}
    pcache::NonlinearFunctionCacheVPRKprojection{DT}

    iguess::InitialGuessIODE{DT,TT,VT,FT,IT}

    q::Vector{Vector{Double{DT}}}
    p::Vector{Vector{Double{DT}}}
end

function IntegratorVPRKpSymmetric(equation::IODE{DT,TT,ΑT,FT,GT,VT}, tableau::TableauVPRK{TT}, Δt::TT;
                                  interpolation=HermiteInterpolation{DT}) where {DT,TT,ΑT,FT,GT,VT}
    D = equation.d
    M = equation.n
    S = tableau.s

    N = D*(S+2)

    if isdefined(tableau, :d)
        tableau_d = tableau.d
    else
        tableau_d = DT[]
    end

    # create solution vectors
    q = Array{Vector{Double{DT}}}(M)
    p = Array{Vector{Double{DT}}}(M)

    for i in 1:M
        q[i] = zeros(Double{DT},D)
        p[i] = zeros(Double{DT},D)
    end

    # create cache for internal stage vectors and update vectors
    scache = NonlinearFunctionCacheVPRK{DT}(D,S)
    pcache = NonlinearFunctionCacheVPRKprojection{DT}(D,S)

    # create params
    params = NonlinearFunctionParametersVPRKpSymmetric{DT,TT,ΑT,FT,GT,D,S}(
                                                equation.α, equation.f, equation.g, Δt,
                                                tableau.o, tableau.q, tableau.p, tableau_d, tableau.R∞)

    # create rhs function for nonlinear solver
    function_stages_solver = (x,b) -> function_stages!(x, b, params)

    # create solver
    solver = get_config(:nls_solver)(zeros(DT,N), function_stages_solver)

    # create initial guess
    iguess = InitialGuessIODE(interpolation, equation, Δt; periodicity=equation.periodicity)

    IntegratorVPRKpSymmetric{DT, TT, ΑT, FT, GT, VT, typeof(params), typeof(solver), typeof(iguess.int)}(
                                        equation, tableau, Δt, params, solver, scache, pcache, iguess, q, p)
end


function initial_guess!(int::IntegratorVPRKpSymmetric{DT,TT}, m::Int) where {DT,TT}
    for i in 1:int.tableau.s
        evaluate!(int.iguess, m, int.scache.y, int.scache.z, int.scache.v, int.tableau.q.c[i], int.tableau.p.c[i])
        for k in 1:int.equation.d
            int.solver.x[int.equation.d*(i-1)+k] = int.scache.v[k]
        end
    end
    evaluate!(int.iguess, m, int.scache.y, int.scache.z, int.scache.v, one(TT), one(TT))
    for k in 1:int.equation.d
        int.solver.x[int.equation.d*(int.tableau.s+0)+k] = int.scache.y[k]
    end
    for k in 1:int.equation.d
        int.solver.x[int.equation.d*(int.tableau.s+1)+k] = 0
    end
end

"Integrate ODE with variational partitioned Runge-Kutta integrator."
function integrate_step!(int::IntegratorVPRKpSymmetric{DT,TT,ΑT,FT,GT,VT}, sol::SolutionPDAE{DT,TT,N}, m::Int, n::Int) where {DT,TT,ΑT,FT,GT,VT,N}
    @assert m ≥ 1
    @assert m ≤ sol.ni

    @assert n ≥ 1
    @assert n ≤ sol.ntime

    # set time for nonlinear solver
    int.params.t = sol.t[0] + (n-1)*int.Δt
    int.params.q .= int.q[m]
    int.params.p .= int.p[m]

    # compute initial guess
    initial_guess!(int, m)

    # call nonlinear solver
    solve!(int.solver)

    # print solver status
    printSolverStatus(int.solver.status, int.solver.params, n)

    # check if solution contains NaNs
    checkNaN(int.solver.status, n)

    # compute vector fields at internal stages and projection vector fields
    compute_stages!(int.solver.x,
                    int.pcache.q̅, int.pcache.p̅, int.pcache.λ,
                    int.scache.Q, int.scache.V, int.pcache.U,
                    int.scache.P, int.scache.F, int.pcache.G, int.params)

    # compute unprojected solution
    update_solution!(int.q[m], int.scache.V, int.tableau.q.b, int.tableau.q.b̂, int.Δt)
    update_solution!(int.p[m], int.scache.F, int.tableau.p.b, int.tableau.p.b̂, int.Δt)

    # add projection to solution
    update_solution!(int.q[m], int.pcache.U, int.params.R, int.Δt)
    update_solution!(int.p[m], int.pcache.G, int.params.R, int.Δt)

    # copy solution to initial guess for next time step
    update!(int.iguess, m, sol.t[0] + n*int.Δt, int.q[m], int.p[m])

    # take care of periodic solutions
    cut_periodic_solution!(int.q[m], int.equation.periodicity)

    # copy to solution
    copy_solution!(sol, int.q[m], int.p[m], int.pcache.λ, n, m)
end
