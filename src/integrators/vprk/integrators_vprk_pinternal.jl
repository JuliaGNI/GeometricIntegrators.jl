
"Parameters for right-hand side function of variational partitioned Runge-Kutta methods."
mutable struct ParametersVPRKpInternal{DT, TT, ET <: IODE{DT,TT}, D, S} <: AbstractParametersVPRK{DT,TT,ET,D,S}
    equ::ET
    tab::TableauVPRK{TT}
    Δt::TT

    t::TT
    q::Vector{DT}
    p::Vector{DT}
end

function ParametersVPRKpInternal(equ::ET, tab::TableauVPRK{TT}, Δt::TT) where {DT, TT, ET <: IODE{DT,TT}}
    D = equ.d
    S = tab.s

    q = zeros(DT,D)
    p = zeros(DT,D)

    ParametersVPRKpInternal{DT,TT,ET,D,S}(equ, tab, Δt, 0, q, p)
end


function compute_stages_vprk!(x, q̅, p̅, Q, V, Λ, P, F, R, params::ParametersVPRKpInternal)
    # copy x to V
    compute_stages_v_vprk!(x, V, params)

    # copy x to Λ
    compute_stages_λ_vprk!(x, Λ, params)

    # compute Q
    compute_stages_q_vprk!(q̅, Q, V, Λ, params)

    # compute p̅ and R
    compute_projection_vprk!(q̅, p̅, Q, V, Λ, R, params)

    # compute P and F
    compute_stages_p_vprk!(Q, V, P, F, params)
end


function compute_stages_λ_vprk!(x::Vector{ST}, Λ::Matrix{ST}, params::ParametersVPRKpInternal{DT,TT,ET,D,S}) where {ST,DT,TT,ET,D,S}
    @assert D == size(Λ,1)
    @assert S == size(Λ,2)

    # copy x to Λ
    for i in 1:S
        for k in 1:D
            Λ[k,i] = x[D*S+k]
        end
    end
end

function compute_stages_q_vprk!(q̅::Vector{ST}, Q::Matrix{ST}, V::Matrix{ST}, Λ::Matrix{ST}, params::ParametersVPRKpInternal{DT,TT,ET,D,S}) where {ST,DT,TT,ET,D,S}
    @assert D == size(Q,1) == size(V,1) == size(Λ,1) == length(q̅)
    @assert S == size(Q,2) == size(V,2) == size(Λ,2)

    local y1::ST
    local y2::ST
    local y3::ST
    local y4::ST

    # compute Q
    for i in 1:S
        for k in 1:D
            y1 = y2 = y3 = y4 = 0
            for j in 1:S
                y1 += params.tab.q.a[i,j] * V[k,j]
                y2 += params.tab.q.â[i,j] * V[k,j]
                y3 += params.tab.q.a[i,j] * Λ[k,j]
                y4 += params.tab.q.â[i,j] * Λ[k,j]
            end
            Q[k,i] = params.q[k] + params.Δt * (y1 + y2 + y3 + y4)
        end
    end

    # compute q̅
    for k in 1:D
        y1 = y2 = y3 = y4 = 0
        for j in 1:S
            y1 += params.tab.q.b[j] * V[k,j]
            y2 += params.tab.q.b̂[j] * V[k,j]
            y3 += params.tab.q.b[j] * Λ[k,j]
            y4 += params.tab.q.b̂[j] * Λ[k,j]
        end
        q̅[k] = params.q[k] + params.Δt * (y1 + y2 + y3 + y4)
    end
end


@generated function compute_projection_vprk!(q̅::Vector{ST}, p̅::Vector{ST},
                Q::Matrix{ST}, V::Matrix{ST}, Λ::Matrix{ST}, R::Matrix{ST},
                params::ParametersVPRKpInternal{DT,TT,ET,D,S}
            ) where {ST,DT,TT,ET,D,S}

    # create temporary vectors
    tQ = zeros(ST,D)
    tV = zeros(ST,D)
    tΛ = zeros(ST,D)
    tR = zeros(ST,D)

    quote
        local t₀::TT = params.t
        local t₁::TT = params.t + params.Δt
        local tᵢ::TT

        # compute p̅=ϑ(q̅)
        $tV .= 0
        params.equ.ϑ(t₁, q̅, $tV, p̅)

        for i in 1:S
            simd_copy_xy_first!($tQ, Q, i)
            simd_copy_xy_first!($tV, V, i)
            simd_copy_xy_first!($tΛ, Λ, i)

            tᵢ = t₀ + params.Δt * params.tab.p.c[i]

            params.equ.g(tᵢ, $tQ, $tΛ, $tR)

            simd_copy_yx_first!($tR, R, i)
        end
    end
end


function compute_rhs_vprk!(b::Vector{ST}, P::Matrix{ST}, F::Matrix{ST}, R::Matrix{ST},
                params::ParametersVPRKpInternal{DT,TT,ET,D,S}) where {ST,DT,TT,ET,D,S}

    local z1::ST
    local z2::ST
    local z3::ST
    local z4::ST

    # compute b = - [(P-AF-AR)]
    for i in 1:S
        for k in 1:D
            z1 = z2 = z3 = z4 = 0
            for j in 1:S
                z1 += params.tab.p.a[i,j] * F[k,j]
                z2 += params.tab.p.â[i,j] * F[k,j]
                z3 += params.tab.p.a[i,j] * R[k,j]
                z4 += params.tab.p.â[i,j] * R[k,j]
            end
            b[D*(i-1)+k] = (P[k,i] - params.p[k]) - params.Δt * (z1 + z2 + z3 + z4)
        end
    end
end

function compute_rhs_vprk_projection!(b::Vector{ST}, p̅::Vector{ST},
                F::Matrix{ST}, R::Matrix{ST}, offset::Int,
                params::ParametersVPRKpInternal{DT,TT,ET,D,S}
            ) where {ST,DT,TT,ET,D,S}

    local z1::ST
    local z2::ST
    local z3::ST
    local z4::ST

    for k in 1:D
        z1 = z2 = z3 = z4 = 0
        for j in 1:S
            z1 += params.tab.p.b[j] * F[k,j]
            z2 += params.tab.p.b̂[j] * F[k,j]
            z3 += params.tab.p.b[j] * R[k,j]
            z4 += params.tab.p.b̂[j] * R[k,j]
        end
        b[offset+k] = (p̅[k] - params.p[k]) - params.Δt * (z1 + z2 + z3 + z4)
    end
end


"Compute stages of variational partitioned Runge-Kutta methods."
@generated function function_stages!(x::Vector{ST}, b::Vector{ST},
                params::ParametersVPRKpInternal{DT,TT,ET,D,S}
            ) where {ST,DT,TT,ET,D,S}

    scache = NonlinearFunctionCacheVPRK{ST}(D,S)
    pcache = NonlinearFunctionCacheVPRKprojection{ST}(D,S)

    function_stages = quote
        compute_stages_vprk!(x, $pcache.q̅, $pcache.p̅,
                                $scache.Q, $scache.V, $pcache.Λ,
                                $scache.P, $scache.F, $pcache.R,
                                params)

        # compute b = [P-AF-AR]
        compute_rhs_vprk!(b, $scache.P, $scache.F, $pcache.R, params)

        # compute b = Φ
        compute_rhs_vprk_projection!(b, $pcache.p̅, $scache.F, $pcache.R, D*S, params)

        compute_rhs_vprk_correction!(b, $scache.V, params)
    end

    return function_stages
end


struct IntegratorVPRKpInternal{DT, TT, PT <: ParametersVPRKpInternal{DT,TT},
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

function IntegratorVPRKpInternal(equation::ET, tableau::TableauVPRK{TT}, Δt::TT) where {DT, TT, ET <: IODE{DT,TT}}
    D = equation.d
    M = equation.n
    S = tableau.s

    # create params
    params = ParametersVPRKpInternal(equation, tableau, Δt)

    # create solver
    solver = create_nonlinear_solver(DT, D*(S+1), params)

    # create initial guess
    iguess = InitialGuessPODE(get_config(:ig_interpolation), equation, Δt)

    # create cache for internal stage vectors and update vectors
    scache = NonlinearFunctionCacheVPRK{DT}(D,S)
    pcache = NonlinearFunctionCacheVPRKprojection{DT}(D,S)

    # create solution vectors
    q = create_solution_vector(DT, D, M)
    p = create_solution_vector(DT, D, M)

    # create integrator
    IntegratorVPRKpInternal{DT, TT, typeof(params), typeof(solver), typeof(iguess)}(
                params, solver, iguess, scache, pcache, q, p)
end


function initial_guess!(int::IntegratorVPRKpInternal{DT,TT}, m::Int) where {DT,TT}
    for i in 1:int.params.tab.s
        evaluate!(int.iguess, m, int.scache.y, int.scache.z, int.scache.v, int.params.tab.q.c[i], int.params.tab.p.c[i])
        for k in 1:int.params.equ.d
            int.solver.x[int.params.equ.d*(i-1)+k] = int.scache.v[k]
        end
    end
    for k in 1:int.params.equ.d
        int.solver.x[int.params.equ.d*int.params.tab.s+k] = 0
    end
end


"Integrate ODE with variational partitioned Runge-Kutta integrator."
function integrate_step!(int::IntegratorVPRKpInternal{DT,TT}, sol::SolutionPDAE{DT,TT}, m::Int, n::Int) where {DT,TT}
    # check if m and n are compatible with solution dimensions
    check_solution_dimension_asserts(sol, m, n)

    # set time and solution for nonlinear solver
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

    # compute final update
    compute_stages_vprk!(int.solver.x, int.pcache.q̅, int.pcache.p̅,
                         int.scache.Q, int.scache.V, int.pcache.Λ,
                         int.scache.P, int.scache.F, int.pcache.R,
                         int.params)

    # compute unprojected solution
    update_solution!(int.q[m], int.scache.V, int.params.tab.q.b, int.params.tab.q.b̂, int.params.Δt)
    update_solution!(int.p[m], int.scache.F, int.params.tab.p.b, int.params.tab.p.b̂, int.params.Δt)

    # add projection to solution
    update_solution!(int.q[m], int.pcache.Λ, int.params.tab.q.b, int.params.tab.q.b̂, int.params.Δt)
    update_solution!(int.p[m], int.pcache.R, int.params.tab.p.b, int.params.tab.p.b̂, int.params.Δt)

    # copy solution to initial guess for next time step
    update!(int.iguess, m, sol.t[0] + n*int.params.Δt, int.q[m], int.p[m])

    # take care of periodic solutions
    cut_periodic_solution!(int.q[m], int.params.equ.periodicity)

    # copy to solution
    copy_solution!(sol, int.q[m], int.p[m], int.pcache.λ, n, m)
end
