
# struct TableauVPRKpLegendre{T} <: AbstractTableauPRK{T}
#     @HeaderTableau
#
#     q::CoefficientsRK{T}
#     p::CoefficientsRK{T}
#
#     R∞::Int
#
#     ω::Matrix{T}
#     d::Vector{T}
#
#     function TableauVPRKpLegendre{T}(name, o, q, p, R∞, ω, d) where {T}
#         @assert q.s == p.s == length(d)
#         new(name, o, q.s, q, p, R∞, ω, d)
#     end
#
#     function TableauVPRKpLegendre{T}(name, o, q, p, R∞, ω) where {T}
#         @assert q.s == p.s
#         new(name, o, q.s, q, p, R∞, ω)
#     end
# end
#
# function TableauVPRKpLegendre(name::Symbol, order::Int, q::CoefficientsRK{T}, p::CoefficientsRK{T}, R∞::Int, ω::Matrix{T}, d::Vector{T}) where {T}
#     TableauVPRKpLegendre{T}(name, order, q, p, R∞, ω, d)
# end
#
# function TableauVPRKpLegendre(name::Symbol, order::Int, q::CoefficientsRK{T}, p::CoefficientsRK{T}, R∞::Int, ω::Matrix{T}) where {T}
#     TableauVPRKpLegendre{T}(name, order, q, p, R∞, ω)
# end
#
# function TableauVPRKpLegendre(name::Symbol, order::Int, q::CoefficientsRK{T}, R∞::Int, ω::Matrix{T}, d::Vector{T}) where {T}
#     TableauVPRKpLegendre{T}(name, order, q, get_symplectic_conjugate_coefficients(q), R∞, ω, d)
# end
#
# function TableauVPRKpLegendre(name::Symbol, order::Int, q::CoefficientsRK{T}, R∞::Int, ω::Matrix{T}) where {T}
#     TableauVPRKpLegendre{T}(name, order, q, get_symplectic_conjugate_coefficients(q), R∞, ω)
# end


"Parameters for right-hand side function of Variational Partitioned Runge-Kutta methods."
const ParametersVPRKpLegendre = AbstractParametersVPRK{:vprk_plegendre}


"Variational special partitioned additive Runge-Kutta integrator."
struct IntegratorVPRKpLegendre{DT, TT, D, S,
                PT <: ParametersVPRKpLegendre{DT,TT},
                ST <: NonlinearSolver{DT},
                IT <: InitialGuessIODE{DT,TT}} <: AbstractIntegratorVPRK{DT,TT,D,S}
    params::PT
    solver::ST
    iguess::IT
    caches::CacheDict{PT}

    function IntegratorVPRKpLegendre(params::ParametersVPRKpLegendre{DT,TT,D,S}, solver::ST, iguess::IT, caches) where {DT,TT,D,S,ST,IT}
        new{DT, TT, D, S, typeof(params), ST, IT}(params, solver, iguess, caches)
    end

    function IntegratorVPRKpLegendre{DT,D}(equations::NamedTuple, tableau::TableauVPRK{TT}, Δt::TT) where {DT,TT,D}
        # get number of stages
        S = tableau.s

        N = (3*S+2)*D

        if isdefined(tableau, :d) && length(tableau.d) > 0
            N += D
        end

        # create params
        params = ParametersVPRKpLegendre{DT,D}(equations, tableau, Δt)

        # create cache dict
        caches = CacheDict(params)

        # create solver
        solver = create_nonlinear_solver(DT, N, params, caches)

        # create initial guess
        iguess = InitialGuessIODE{DT,D}(get_config(:ig_interpolation), equations[:v̄], equations[:f̄], Δt)

        # create integrator
        IntegratorVPRKpLegendre(params, solver, iguess, caches)
    end

    function IntegratorVPRKpLegendre(equation::IODE{DT,TT}, tableau::TableauVPRK{TT}, Δt::TT; kwargs...) where {DT,TT}
        IntegratorVPRKpLegendre{DT, equation.d}(get_function_tuple(equation), tableau, Δt; kwargs...)
    end
end


function Integrators.initialize!(int::IntegratorVPRKpLegendre, sol::AtomicSolutionPODE)
    sol.t̅ = sol.t - timestep(int)

    equation(int, :v̄)(sol.t, sol.q, sol.v)
    equation(int, :f̄)(sol.t, sol.q, sol.v, sol.f)

    initialize!(int.iguess, sol.t, sol.q, sol.p, sol.v, sol.f,
                            sol.t̅, sol.q̅, sol.p̅, sol.v̅, sol.f̅)
end


function initial_guess!(int::IntegratorVPRKpLegendre{DT,TT}, sol::AtomicSolutionPODE{DT,TT},
                        cache::IntegratorCacheVPRK{DT}=int.caches[DT]) where {DT,TT}
    for i in eachstage(int)
        evaluate!(int.iguess, sol.q̅, sol.p̅, sol.v̅, sol.f̅,
                              sol.q, sol.p, sol.v, sol.f,
                              cache.q̃, cache.p̃, cache.v, cache.f,
                              tableau(int).q.c[i], tableau(int).p.c[i])

        for k in eachdim(int)
            offset = 3*(ndims(int)*(i-1)+k-1)
            int.solver.x[offset+1] = (cache.q̃[k] - sol.q[k]) / timestep(int)
            int.solver.x[offset+2] = (cache.p̃[k] - sol.p[k]) / timestep(int)
            int.solver.x[offset+3] = cache.v[k]
        end
    end

    evaluate!(int.iguess, sol.q, sol.p, sol.v, sol.f,
                          sol.q̅, sol.p̅, sol.v̅, sol.f̅,
                          cache.q̃, cache.p̃,
                          one(TT), one(TT))

    offset = (3*nstages(int)+0)*ndims(int)
    for k in eachdim(int)
        int.solver.x[offset+k] = cache.q̃[k]
    end

    offset = (3*nstages(int)+1)*ndims(int)
    for k in eachdim(int)
        int.solver.x[offset+k] = cache.p̃[k]
    end

    if isdefined(tableau(int), :d)  && length(tableau(int).d) > 0
        offset = (3*nstages(int)+2)*ndims(int)
        for k in eachdim(int)
            int.solver.x[offset+k] = 0
        end
    end
end


function compute_stages!(y::Vector{ST}, Q, V, P, F, Y, Z, Φ, q, v, p, ϕ, μ,
                params::ParametersVPRKpLegendre{DT,TT,D,S}) where {ST,DT,TT,D,S}

    local offset::Int
    local tqᵢ::TT

    for i in 1:S
        # copy y to Q and V
        for k in 1:D
            offset = 3*(D*(i-1)+k-1)
            Y[i][k] = y[offset+1]
            Z[i][k] = y[offset+2]
            V[i][k] = y[offset+3]

            Q[i][k] = params.q̅[k] + params.Δt * Y[i][k]
            P[i][k] = params.p̅[k] + params.Δt * Z[i][k]
        end

        # compute P=p(t,Q) and F=f(t,Q,V)
        tqᵢ = params.t̅ + params.Δt * params.tab.q.c[i]
        params.equ[:ϑ](tqᵢ, Q[i], V[i], Φ[i])
        params.equ[:f](tqᵢ, Q[i], V[i], F[i])
    end

    # copy y to q and p
    for k in 1:D
        q[k] = y[(3*S+0)*D+k]
        p[k] = y[(3*S+1)*D+k]
    end

    # compute p=p(t,q)
    params.equ[:ϑ](params.t̅ + params.Δt, q, v, ϕ)

    # for Lobatto-type methods, copy y to μ
    if isdefined(params.tab, :d) && length(params.tab.d) > 0
        offset = (3*S+2)*D
        for k in 1:D
            μ[k] = y[offset+k]
        end
    end
end


function compute_rhs!(b::Vector{ST},
                      Q::Vector{Vector{ST}}, V::Vector{Vector{ST}},
                      P::Vector{Vector{ST}}, F::Vector{Vector{ST}},
                      Y::Vector{Vector{ST}}, Z::Vector{Vector{ST}},
                      Φ::Vector{Vector{ST}},
                      q::Vector{ST}, p::Vector{ST}, ϕ::Vector{ST}, μ::Vector{ST},
                      params::ParametersVPRKpLegendre{DT,TT,D,S}) where {ST,DT,TT,D,S}
    local ω::ST
    local y::ST
    local z::ST

    # compute b = - (Q-q-AV)
    for i in 1:S
        offset = 0*S*D + D*(i-1)
        for k in 1:D
            y = 0
            for j in 1:S
                y += params.tab.q.a[i,j] * V[j][k]
            end
            b[offset+k] = Y[i][k] - y
        end
    end

    # compute b = - (P-p-AF)
    for i in 1:S
        offset = 1*S*D + D*(i-1)
        for k in 1:D
            z = 0
            for j in 1:S
                z += params.tab.p.a[i,j] * F[j][k]
            end
            b[offset+k] = Z[i][k] - z
        end
    end

    # compute b = - (P-Φ)
    for i in 1:S-1
        offset = 2*S*D + D*(i-1)
        for k in 1:D
            z = 0
            for j in 1:S
                ω = params.tab.p.b[j] * params.tab.p.c[j]^(i-1)
                # ω = params.tab.q.a[i+1,j]
                z += ω * (Φ[j][k] - P[j][k])
            end
            b[offset+k] = z
        end
    end

    # compute b = - (p-ϕ)
    offset = (3*S-1)*D
    for k in 1:D
        b[offset+k] = ϕ[k] - p[k]
    end

    # compute b = - (q-q-AV)
    offset = (3*S+0)*D
    for k in 1:D
        y = 0
        for j in 1:S
            y += params.tab.q.b[j] * V[j][k]
        end
        b[offset+k] = (params.q̅[k] - q[k]) + params.Δt * y
    end

    # compute b = - (p̅-p-AF)
    offset = (3*S+1)*D
    for k in 1:D
        z = 0
        for j in 1:S
            z += params.tab.p.b[j] * F[j][k]
        end
        b[offset+k] = (params.p̅[k] - p[k]) + params.Δt * z
    end

    # if isdefined(params.tab, :d) && length(params.tab.d) > 0
    #     offset = (3*S+2)*D
    #     for k in 1:D
    #         b[offset+k] = 0
    #         for i in 1:S
    #             b[offset+k] -= V[k,i] * params.tab.d[i]
    #         end
    #     end
    # end

    if isdefined(params.tab, :d) && length(params.tab.d) > 0
        sl     = div(S+1, 2)
        offset = D*(S+sl-1)

        # compute μ
        z = params.tab.p.b[sl] / params.tab.d[sl]
        for k in 1:D
            μ[k] = z * b[offset+k]
        end

        # replace equation for Pₗ with constraint on V
        for k in 1:D
            z = 0
            for i in 1:S
                z += V[k,i] * params.tab.d[i]
            end
            b[offset+k] = z
        end

        # modify P₁, ..., Pₛ except for Pₗ
        for i in 1:S
            if i ≠ sl
                offset = D*(S+i-1)
                z = params.tab.d[i] / params.tab.p.b[i]
                for k in 1:D
                    b[offset+k] -= z * μ[k]
                end
            end
        end
    end
end

"Compute stages of variational special partitioned additive Runge-Kutta methods."
function Integrators.function_stages!(y::Vector{ST}, b::Vector{ST},
                params::ParametersVPRKpLegendre{DT,TT,D,S},
                caches::CacheDict) where {ST,DT,TT,D,S}

    @assert length(y) == length(b)

    # get cache for internal stages
    cache = caches[ST]

    compute_stages!(y, cache.Q, cache.V, cache.P, cache.F, cache.Y, cache.Z, cache.Φ, cache.q̃, cache.ṽ, cache.p̃, cache.ϕ, cache.μ, params)

    compute_rhs!(b, cache.Q, cache.V, cache.P, cache.F, cache.Y, cache.Z, cache.Φ, cache.q̃, cache.p̃, cache.ϕ, cache.μ, params)

    # debug output
    # println()
    # for k in 1:D
    #     println(params.q[k], ",  ", params.p[k], ",  ", cache.Q[k,1], ",  ", cache.P[k,1], ",  ", cache.q[k], ",  ", cache.p[k])
    # end
    # println()
end


"Integrate DAE with variational special partitioned additive Runge-Kutta integrator."
function Integrators.integrate_step!(int::IntegratorVPRKpLegendre{DT,TT}, sol::AtomicSolutionPODE{DT,TT},
                                     cache::IntegratorCacheVPRK{DT}=int.caches[DT]) where {DT,TT}
    local offset::Int

    # update nonlinear solver parameters from cache
    update_params!(int.params, sol)

    # compute initial guess
    initial_guess!(int, sol, cache)

    # debug output
    # println()
    # println(int.solver.x)

    # reset solution
    reset!(sol, timestep(int))

    # call nonlinear solver
    solve!(int.solver)

    # print solver status
    print_solver_status(int.solver.status, int.solver.params)

    # check if solution contains NaNs or error bounds are violated
    check_solver_status(int.solver.status, int.solver.params)

    # compute vector fields at internal stages
    compute_stages!(int.solver.x, cache.Q, cache.V, cache.P,
                    cache.F, cache.Y, cache.Z, cache.Φ,
                    cache.q̃, cache.ṽ, cache.p̃, cache.ϕ,
                    cache.μ, int.params)

    # compute final update
    update_solution!(sol.q, cache.V, tableau(int).q.b, tableau(int).q.b̂, timestep(int))
    update_solution!(sol.p, cache.F, tableau(int).p.b, tableau(int).p.b̂, timestep(int))

    # copy solution to initial guess for next time step
    update_vector_fields!(int.iguess, sol.t, sol.q, sol.p, sol.v, sol.f)
end
