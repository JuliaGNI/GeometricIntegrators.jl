
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


struct IntegratorVPRKpInternal{DT, TT, PT <: ParametersVPRKpInternal{DT,TT},
                                        ST <: NonlinearSolver{DT},
                                        IT <: InitialGuessPODE{DT,TT}, D, S} <: AbstractIntegratorVPRK{DT,TT}
    params::PT
    solver::ST
    iguess::IT
    cache::IntegratorCacheVPRK{DT,D,S}
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
    cache = IntegratorCacheVPRK{DT,D,S}(true)

    # create integrator
    IntegratorVPRKpInternal{DT, TT, typeof(params), typeof(solver), typeof(iguess), D, S}(
                params, solver, iguess, cache)
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

    cache = IntegratorCacheVPRK{ST, D, S}()

    function_stages = quote
        compute_stages_vprk!(x, $cache.q̃, $cache.p̃,
                                $cache.Q, $cache.V, $cache.Λ,
                                $cache.P, $cache.F, $cache.R,
                                params)

        # compute b = [P-AF-AR]
        compute_rhs_vprk!(b, $cache.P, $cache.F, $cache.R, params)

        # compute b = Φ
        compute_rhs_vprk_projection!(b, $cache.p̃, $cache.F, $cache.R, D*S, params)

        compute_rhs_vprk_correction!(b, $cache.V, params)
    end

    return function_stages
end


function initial_guess!(int::IntegratorVPRKpInternal{DT,TT}, sol::AtomisticSolutionPODE{DT,TT}) where {DT,TT}
    for i in eachstage(int)
        evaluate!(int.iguess, sol.q, sol.p, sol.v, sol.f,
                              sol.q̅, sol.p̅, sol.v̅, sol.f̅,
                              int.cache.q̃, int.cache.ṽ,
                              tableau(int).q.c[i])
        for k in eachdim(int)
            int.solver.x[ndims(int)*(i-1)+k] = int.cache.ṽ[k]
        end
    end
    for k in eachdim(int)
        int.solver.x[ndims(int)*nstages(int)+k] = 0
    end
end


"Integrate ODE with variational partitioned Runge-Kutta integrator."
function integrate_step!(int::IntegratorVPRKpInternal{DT,TT}, sol::AtomisticSolutionPODE{DT,TT}) where {DT,TT}
    # update nonlinear solver parameters from cache
    update_params!(int.params, sol)

    # compute initial guess
    initial_guess!(int, sol)

    # reset solution
    reset!(sol, timestep(int))

    # call nonlinear solver
    solve!(int.solver)

    # print solver status
    print_solver_status(int.solver.status, int.solver.params)

    # check if solution contains NaNs or error bounds are violated
    check_solver_status(int.solver.status, int.solver.params)

    # compute final update
    compute_stages_vprk!(int.solver.x, int.cache.q̃, int.cache.p̃,
                          int.cache.Q, int.cache.V, int.cache.Λ,
                          int.cache.P, int.cache.F, int.cache.R,
                          int.params)

    # compute unprojected solution
    update_solution!(sol.q, int.cache.V, tableau(int).q.b, tableau(int).q.b̂, timestep(int))
    update_solution!(sol.p, int.cache.F, tableau(int).p.b, tableau(int).p.b̂, timestep(int))

    # add projection to solution
    update_solution!(sol.q, int.cache.Λ, tableau(int).q.b, tableau(int).q.b̂, timestep(int))
    update_solution!(sol.p, int.cache.R, tableau(int).p.b, tableau(int).p.b̂, timestep(int))

    # copy solution to initial guess
    update!(int.iguess, sol.t, sol.q, sol.p, sol.v, sol.f)
end
