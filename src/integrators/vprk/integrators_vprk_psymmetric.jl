
"Parameters for right-hand side function of variational partitioned Runge-Kutta methods."
mutable struct ParametersVPRKpSymmetric{DT, TT, ET <: IODE{DT,TT}, D, S} <: AbstractParametersVPRK{DT,TT,ET,D,S}
    equ::ET
    tab::TableauVPRK{TT}
    Δt::TT

    R::Vector{TT}

    t::TT
    q::Vector{DT}
    p::Vector{DT}
end

function ParametersVPRKpSymmetric(equ::ET, tab::TableauVPRK{TT}, Δt::TT) where {DT, TT, ET <: IODE{DT,TT}}
    R = convert(Vector{TT}, [1, tab.R∞])
    q = zeros(DT, equ.d)
    p = zeros(DT, equ.d)

    ParametersVPRKpSymmetric{DT, TT, ET, equ.d, tab.s}(equ, tab, Δt, R, 0, q, p)
end


@generated function compute_projection_vprk!(x::Vector{ST},
                q̅::Vector{ST}, p̅::Vector{ST}, λ::Vector{ST},
                V::Matrix{ST}, U::Matrix{ST}, G::Matrix{ST},
                params::ParametersVPRKpSymmetric{DT,TT,ET,D,S}
            ) where {ST,DT,TT,ET,D,S}

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
        params.equ.g(t₀, params.q, λ, $tG)
        simd_copy_yx_first!($tG, G, 1)

        params.equ.g(t₁, q̅, λ, $tG)
        simd_copy_yx_first!($tG, G, 2)

        # compute p̅=α(q̅)p̅
        params.equ.α(t₁, q̅, λ, p̅)
    end
end


"Compute stages of variational partitioned Runge-Kutta methods."
@generated function function_stages!(x::Vector{ST}, b::Vector{ST},
                params::ParametersVPRKpSymmetric{DT,TT,ET,D,S}
            ) where {ST,DT,TT,ET,D,S}

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
struct IntegratorVPRKpSymmetric{DT, TT, PT <: ParametersVPRKpSymmetric{DT,TT},
                                        ST <: NonlinearSolver{DT},
                                        IT <: InitialGuessPODE{DT,TT}} <: AbstractIntegratorVPRK{DT,TT}

    params::PT
    solver::ST
    iguess::IT

    scache::NonlinearFunctionCacheVPRK{DT}
    pcache::NonlinearFunctionCacheVPRKprojection{DT}

    q::Vector{Vector{TwicePrecision{DT}}}
    p::Vector{Vector{TwicePrecision{DT}}}
end

function IntegratorVPRKpSymmetric(equation::ET, tableau::TableauVPRK{TT}, Δt::TT) where {DT, TT, ET <: IODE{DT,TT}}
    D = equation.d
    M = equation.n
    S = tableau.s

    # create params
    params = ParametersVPRKpSymmetric(equation, tableau, Δt)

    # create solver
    solver = create_nonlinear_solver(DT, D*(S+2), params)

    # create initial guess
    iguess = InitialGuessPODE(get_config(:ig_interpolation), equation, Δt)

    # create cache for internal stage vectors and update vectors
    scache = NonlinearFunctionCacheVPRK{DT}(D,S)
    pcache = NonlinearFunctionCacheVPRKprojection{DT}(D,S)

    # create solution vectors
    q = create_solution_vector_double_double(DT, D, M)
    p = create_solution_vector_double_double(DT, D, M)

    # create integrator
    IntegratorVPRKpSymmetric{DT, TT, typeof(params), typeof(solver), typeof(iguess)}(
                params, solver, iguess, scache, pcache, q, p)
end


function initial_guess!(int::IntegratorVPRKpSymmetric{DT,TT}, m::Int) where {DT,TT}
    for i in 1:int.params.tab.s
        evaluate!(int.iguess, m, int.scache.y, int.scache.z, int.scache.v, int.params.tab.q.c[i], int.params.tab.p.c[i])
        for k in 1:int.params.equ.d
            int.solver.x[int.params.equ.d*(i-1)+k] = int.scache.v[k]
        end
    end
    evaluate!(int.iguess, m, int.scache.y, int.scache.z, int.scache.v, one(TT), one(TT))
    for k in 1:int.params.equ.d
        int.solver.x[int.params.equ.d*(int.params.tab.s+0)+k] = int.scache.y[k]
    end
    for k in 1:int.params.equ.d
        int.solver.x[int.params.equ.d*(int.params.tab.s+1)+k] = 0
    end
end

"Integrate ODE with variational partitioned Runge-Kutta integrator."
function integrate_step!(int::IntegratorVPRKpSymmetric{DT,TT}, sol::SolutionPDAE{DT,TT}, m::Int, n::Int) where {DT,TT}
    # check if m and n are compatible with solution dimensions
    check_solution_dimension_asserts(sol, m, n)

    # set time for nonlinear solver
    int.params.t = sol.t[0] + (n-1)*int.params.Δt
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

    # compute vector fields at internal stages and projection vector fields
    compute_stages!(int.solver.x,
                    int.pcache.q̅, int.pcache.p̅, int.pcache.λ,
                    int.scache.Q, int.scache.V, int.pcache.U,
                    int.scache.P, int.scache.F, int.pcache.G, int.params)

    # compute unprojected solution
    update_solution!(int.q[m], int.scache.V, int.params.tab.q.b, int.params.tab.q.b̂, int.params.Δt)
    update_solution!(int.p[m], int.scache.F, int.params.tab.p.b, int.params.tab.p.b̂, int.params.Δt)

    # add projection to solution
    update_solution!(int.q[m], int.pcache.U, int.params.R, int.params.Δt)
    update_solution!(int.p[m], int.pcache.G, int.params.R, int.params.Δt)

    # copy solution to initial guess for next time step
    update!(int.iguess, m, sol.t[0] + n*int.params.Δt, int.q[m], int.p[m])

    # take care of periodic solutions
    cut_periodic_solution!(int.q[m], int.params.equ.periodicity)

    # copy to solution
    copy_solution!(sol, int.q[m], int.p[m], int.pcache.λ, n, m)
end
