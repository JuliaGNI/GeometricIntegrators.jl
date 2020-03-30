
"Parameters for right-hand side function of Variational Partitioned Runge-Kutta methods."
const ParametersVPRKpInternal = AbstractParametersVPRK{:vprk_pinternal}


struct IntegratorVPRKpInternal{DT, TT, D, S,
                PT <: ParametersVPRKpInternal{DT,TT},
                ST <: NonlinearSolver{DT},
                IT <: InitialGuessIODE{DT,TT}} <: AbstractIntegratorVPRK{DT,TT,D,S}
    params::PT
    solver::ST
    iguess::IT
    caches::CacheDict{PT}

    function IntegratorVPRKpInternal(params::ParametersVPRKpInternal{DT,TT,D,S}, solver::ST, iguess::IT, caches) where {DT,TT,D,S,ST,IT}
        new{DT, TT, D, S, typeof(params), ST, IT}(params, solver, iguess, caches)
    end

    function IntegratorVPRKpInternal{DT,D}(equations::NamedTuple, tableau::TableauVPRK{TT}, Δt::TT) where {DT,TT,D}
        # get number of stages
        S = tableau.s

        # create params
        params = ParametersVPRKpInternal{DT,D}(equations, tableau, Δt)

        # create cache dict
        caches = CacheDict(params)

        # create solver
        solver = create_nonlinear_solver(DT, D*(S+1), params, caches)

        # create initial guess
        iguess = InitialGuessIODE{DT,D}(get_config(:ig_interpolation), equations[:v], equations[:f], Δt)

        # create integrator
        IntegratorVPRKpInternal(params, solver, iguess, caches)
    end

    function IntegratorVPRKpInternal(equation::IODE{DT,TT}, tableau::TableauVPRK{TT}, Δt::TT; kwargs...) where {DT,TT}
        IntegratorVPRKpInternal{DT, ndims(equation)}(get_function_tuple(equation), tableau, Δt; kwargs...)
    end
end


function initial_guess!(int::IntegratorVPRKpInternal{DT,TT}, sol::AtomicSolutionPODE{DT,TT},
                        cache::IntegratorCacheVPRK{DT}=int.caches[DT]) where {DT,TT}
    for i in eachstage(int)
        evaluate!(int.iguess, sol.q, sol.p, sol.v, sol.f,
                              sol.q̅, sol.p̅, sol.v̅, sol.f̅,
                              cache.q̃, cache.ṽ,
                              tableau(int).q.c[i])
        for k in eachdim(int)
            int.solver.x[ndims(int)*(i-1)+k] = cache.ṽ[k]
        end
    end
    for k in eachdim(int)
        int.solver.x[ndims(int)*nstages(int)+k] = 0
    end
end


function compute_stages_vprk!(x, q, v, p, Q, V, Λ, P, F, R, params::ParametersVPRKpInternal)
    # copy x to V
    compute_stages_v_vprk!(x, V, params)

    # copy x to Λ
    compute_stages_λ_vprk!(x, Λ, params)

    # compute Q
    compute_stages_q_vprk!(q, Q, V, Λ, params)

    # compute p̅ and R
    compute_projection_vprk!(q, v, p, Q, V, Λ, R, params)

    # compute P and F
    compute_stages_p_vprk!(Q, V, P, F, params)
end


function compute_stages_λ_vprk!(x::Vector{ST}, Λ::Vector{Vector{ST}},
                params::ParametersVPRKpInternal{DT,TT,D,S}) where {ST,DT,TT,D,S}
    # @assert D == size(Λ,1)
    @assert S == length(Λ)

    # copy x to Λ
    for i in 1:S
        for k in 1:D
            Λ[i][k] = x[D*S+k]
        end
    end
end

function compute_stages_q_vprk!(q::Vector{ST}, Q::Vector{Vector{ST}}, V::Vector{Vector{ST}},
                Λ::Vector{Vector{ST}}, params::ParametersVPRKpInternal{DT,TT,D,S}) where {ST,DT,TT,D,S}
    # @assert D == size(Q,1) == size(V,1) == size(Λ,1) == length(q)
    @assert S == length(Q) == length(V) == length(Λ)

    local y1::ST
    local y2::ST
    local y3::ST
    local y4::ST

    # compute Q
    for i in 1:S
        for k in 1:D
            y1 = y2 = y3 = y4 = 0
            for j in 1:S
                y1 += params.tab.q.a[i,j] * V[j][k]
                y2 += params.tab.q.â[i,j] * V[j][k]
                y3 += params.tab.q.a[i,j] * Λ[j][k]
                y4 += params.tab.q.â[i,j] * Λ[j][k]
            end
            Q[i][k] = params.q̅[k] + params.Δt * (y1 + y2 + y3 + y4)
        end
    end

    # compute q
    for k in 1:D
        y1 = y2 = y3 = y4 = 0
        for j in 1:S
            y1 += params.tab.q.b[j] * V[j][k]
            y2 += params.tab.q.b̂[j] * V[j][k]
            # y3 += params.tab.q.b[j] * Λ[j][k]
            # y4 += params.tab.q.b̂[j] * Λ[j][k]
        end
        q[k] = params.q̅[k] + params.Δt * (y1 + y2 + y3 + y4)
    end
end


function compute_projection_vprk!(q::Vector{ST}, v::Vector{ST}, p::Vector{ST},
                Q::Vector{Vector{ST}}, V::Vector{Vector{ST}},
                Λ::Vector{Vector{ST}}, R::Vector{Vector{ST}},
                params::ParametersVPRKpInternal{DT,TT,D,S}) where {ST,DT,TT,D,S}

    local t₀::TT = params.t̅
    local t₁::TT = params.t̅ + params.Δt
    local tᵢ::TT

    # compute p=ϑ(q)
    params.equ.ϑ(t₁, q, v, p)

    for i in 1:S
        tᵢ = t₀ + params.Δt * params.tab.p.c[i]
        params.equ.g(tᵢ, Q[i], Λ[i], R[i])
    end
end


function compute_rhs_vprk!(b::Vector{ST}, P::Vector{Vector{ST}},
                F::Vector{Vector{ST}}, R::Vector{Vector{ST}},
                params::ParametersVPRKpInternal{DT,TT,D,S}) where {ST,DT,TT,D,S}

    local z1::ST
    local z2::ST
    local z3::ST
    local z4::ST

    # compute b = - [(P-AF-AR)]
    for i in 1:S
        for k in 1:D
            z1 = z2 = z3 = z4 = 0
            for j in 1:S
                z1 += params.tab.p.a[i,j] * F[j][k]
                z2 += params.tab.p.â[i,j] * F[j][k]
                z3 += params.tab.p.a[i,j] * R[j][k]
                z4 += params.tab.p.â[i,j] * R[j][k]
            end
            b[D*(i-1)+k] = (P[i][k] - params.p̅[k]) - params.Δt * (z1 + z2 + z3 + z4)
        end
    end
end

function compute_rhs_vprk_projection!(b::Vector{ST}, p::Vector{ST},
                F::Vector{Vector{ST}}, R::Vector{Vector{ST}}, offset::Int,
                params::ParametersVPRKpInternal{DT,TT,D,S}) where {ST,DT,TT,D,S}

    local z1::ST
    local z2::ST
    local z3::ST
    local z4::ST

    for k in 1:D
        z1 = z2 = z3 = z4 = 0
        for j in 1:S
            z1 += params.tab.p.b[j] * F[j][k]
            z2 += params.tab.p.b̂[j] * F[j][k]
            # z3 += params.tab.p.b[j] * R[j][k]
            # z4 += params.tab.p.b̂[j] * R[j][k]
        end
        b[offset+k] = (p[k] - params.p̅[k]) - params.Δt * (z1 + z2 + z3 + z4)
    end
end


"Compute stages of variational partitioned Runge-Kutta methods."
function Integrators.function_stages!(x::Vector{ST}, b::Vector{ST},
                params::ParametersVPRKpInternal{DT,TT,D,S},
                caches::CacheDict) where {ST,DT,TT,D,S}

    # get cache for internal stages
    cache = caches[ST]

    compute_stages_vprk!(x, cache.q̃, cache.ṽ, cache.p̃,
                            cache.Q, cache.V, cache.Λ,
                            cache.P, cache.F, cache.R,
                            params)

    # compute b = [P-AF-AR]
    compute_rhs_vprk!(b, cache.P, cache.F, cache.R, params)

    # compute b = Φ
    compute_rhs_vprk_projection!(b, cache.p̃, cache.F, cache.R, D*S, params)

    compute_rhs_vprk_correction!(b, cache.V, params)
end


"Integrate ODE with variational partitioned Runge-Kutta integrator."
function Integrators.integrate_step!(int::IntegratorVPRKpInternal{DT,TT}, sol::AtomicSolutionPODE{DT,TT},
                                     cache::IntegratorCacheVPRK{DT}=int.caches[DT]) where {DT,TT}
    # update nonlinear solver parameters from cache
    update_params!(int.params, sol)

    # compute initial guess
    initial_guess!(int, sol, cache)

    # reset solution
    reset!(sol, timestep(int))

    # call nonlinear solver
    solve!(int.solver)

    # print solver status
    print_solver_status(int.solver.status, int.solver.params)

    # check if solution contains NaNs or error bounds are violated
    check_solver_status(int.solver.status, int.solver.params)

    # compute final update
    compute_stages_vprk!(int.solver.x,
                         cache.q̃, cache.ṽ,cache.p̃,
                         cache.Q, cache.V, cache.Λ,
                         cache.P, cache.F, cache.R,
                         int.params)

    # compute unprojected solution
    update_solution!(sol.q, sol.q̃, cache.V, tableau(int).q.b, tableau(int).q.b̂, timestep(int))
    update_solution!(sol.p, sol.p̃, cache.F, tableau(int).p.b, tableau(int).p.b̂, timestep(int))

    # add projection to solution
    # update_solution!(sol.q, sol.q̃, cache.Λ, tableau(int).q.b, tableau(int).q.b̂, timestep(int), cache)
    # update_solution!(sol.p, sol.p̃, cache.R, tableau(int).p.b, tableau(int).p.b̂, timestep(int), cache)

    # copy solution to initial guess
    update_vector_fields!(int.iguess, sol.t, sol.q, sol.p, sol.v, sol.f)
end
