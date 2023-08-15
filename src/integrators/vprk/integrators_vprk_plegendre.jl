
# struct TableauVPRKpLegendre{T} <: AbstractPartitionedTableau{T}
#     @HeaderTableau
#
#     q::Tableau{T}
#     p::Tableau{T}
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
# function TableauVPRKpLegendre(name::Symbol, order::Int, q::Tableau{T}, p::Tableau{T}, R∞::Int, ω::Matrix{T}, d::Vector{T}) where {T}
#     TableauVPRKpLegendre{T}(name, order, q, p, R∞, ω, d)
# end
#
# function TableauVPRKpLegendre(name::Symbol, order::Int, q::Tableau{T}, p::Tableau{T}, R∞::Int, ω::Matrix{T}) where {T}
#     TableauVPRKpLegendre{T}(name, order, q, p, R∞, ω)
# end
#
# function TableauVPRKpLegendre(name::Symbol, order::Int, q::Tableau{T}, R∞::Int, ω::Matrix{T}, d::Vector{T}) where {T}
#     TableauVPRKpLegendre{T}(name, order, q, get_symplectic_conjugate_coefficients(q), R∞, ω, d)
# end
#
# function TableauVPRKpLegendre(name::Symbol, order::Int, q::Tableau{T}, R∞::Int, ω::Matrix{T}) where {T}
#     TableauVPRKpLegendre{T}(name, order, q, get_symplectic_conjugate_coefficients(q), R∞, ω)
# end


function _VPRKpLegendre_ndofs(D, tableau::PartitionedTableau)
    S = tableau.s

    N = (3*S+2)*D

    if isdefined(tableau, :d) && length(tableau.d) > 0
        N += D
    end

    return N
end

"Parameters for right-hand side function of Variational Partitioned Runge-Kutta methods."
const ParametersVPRKpLegendre = AbstractParametersVPRK{:vprk_plegendre}


function IntegratorCache(params::ParametersVPRKpLegendre{DT,TT,D,S}; kwargs...) where {DT,TT,D,S}
    IntegratorCacheVPRK{DT,D,S}(_VPRKpLegendre_ndofs(D, tableau(method)), true; kwargs...)
end

function IntegratorCache{ST}(params::ParametersVPRKpLegendre{DT,TT,D,S}; kwargs...) where {ST,DT,TT,D,S}
    IntegratorCacheVPRK{ST,D,S}(_VPRKpLegendre_ndofs(D, tableau(method)), true; kwargs...)
end


"Variational special partitioned additive Runge-Kutta integrator."
struct IntegratorVPRKpLegendre{DT, TT, D, S,
                PT <: ParametersVPRKpLegendre{DT,TT},
                ST <: NonlinearSolver,
                IT <: InitialGuessIODE{TT}} <: GeometricIntegratorVPRK{DT,TT,D,S}
    params::PT
    solver::ST
    iguess::IT
    caches::OldCacheDict{PT}

    function IntegratorVPRKpLegendre(params::ParametersVPRKpLegendre{DT,TT,D,S}, solver::ST, iguess::IT, caches) where {DT,TT,D,S,ST,IT}
        new{DT, TT, D, S, typeof(params), ST, IT}(params, solver, iguess, caches)
    end

    function IntegratorVPRKpLegendre{DT,D}(equations::NamedTuple, tableau::PartitionedTableau{TT}, nullvec, Δt::TT) where {DT,TT,D}
        # get number of variables for nonlinear solver
        N = _VPRKpLegendre_ndofs(D, tableau)

        # create params
        params = ParametersVPRKpLegendre{DT,D}(equations, tableau, nullvec, Δt)

        # create cache dict
        caches = OldCacheDict(params)

        # create solver
        solver = create_nonlinear_solver(DT, N, params, caches)

        # create initial guess
        iguess = InitialGuessIODE(get_config(:ig_extrapolation), equations[:v̄], equations[:f̄], Δt)

        # create integrator
        IntegratorVPRKpLegendre(params, solver, iguess, caches)
    end

    function IntegratorVPRKpLegendre(problem::Union{IODEProblem{DT},LODEProblem{DT}}, tableau, nullvec; kwargs...) where {DT}
        IntegratorVPRKpLegendre{DT, ndims(problem)}(functions(problem), tableau, nullvec, timestep(problem); kwargs...)
    end
end


function Integrators.initialize!(int::IntegratorVPRKpLegendre, sol::SolutionStepPODE)
    equation(int, :v̄)(sol.v, sol.t, sol.q)
    equation(int, :f̄)(sol.f, sol.t, sol.q, sol.v)

    initialize!(int.iguess, sol.t, sol.q, sol.p, sol.v, sol.f,
                            sol.t̄, sol.q̄, sol.p̄, sol.v̄, sol.f̄)
end


function initial_guess!(int::IntegratorVPRKpLegendre{DT,TT}, sol::SolutionStepPODE{DT,TT},
                        cache::IntegratorCacheVPRK{DT}=int.caches[DT]) where {DT,TT}
    for i in eachstage(int)
        evaluate!(int.iguess, sol.history.q[2], sol.history.p[2], sol.history.v[2], sol.history.f[2],
                              sol.history.q[1], sol.history.p[1], sol.history.v[1], sol.history.f[1],
                              cache.q̃, cache.p̃, cache.v, cache.f,
                              tableau(int).q.c[i], tableau(int).p.c[i])

        for k in eachdim(int)
            offset = 3*(ndims(int)*(i-1)+k-1)
            cache.x[offset+1] = (cache.q̃[k] - sol.q[k]) / timestep(int)
            cache.x[offset+2] = (cache.p̃[k] - sol.p[k]) / timestep(int)
            cache.x[offset+3] = cache.v[k]
        end
    end

    evaluate!(int.iguess, sol.history.q[2], sol.history.p[2], sol.history.v[2], sol.history.f[2],
                          sol.history.q[1], sol.history.p[1], sol.history.v[1], sol.history.f[1],
                          cache.q̃, cache.p̃,
                          one(TT), one(TT))

    offset = (3*nstages(int)+0)*ndims(int)
    for k in eachdim(int)
        cache.x[offset+k] = cache.q̃[k]
    end

    offset = (3*nstages(int)+1)*ndims(int)
    for k in eachdim(int)
        cache.x[offset+k] = cache.p̃[k]
    end

    if isdefined(tableau(int), :d)  && length(tableau(int).d) > 0
        offset = (3*nstages(int)+2)*ndims(int)
        for k in eachdim(int)
            cache.x[offset+k] = 0
        end
    end
end


function components!(y::Vector{ST}, Q, V, P, F, Y, Z, Φ, q, v, p, ϕ, μ,
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

            Q[i][k] = solstep.q̄[k] + timestep(problem) * Y[i][k]
            P[i][k] = solstep.p̄[k] + timestep(problem) * Z[i][k]
        end

        # compute P=p(t,Q) and F=f(t,Q,V)
        tqᵢ = solstep.t̄ + timestep(problem) * tableau(method).q.c[i]
        functions(problem).ϑ(Φ[i], tqᵢ, Q[i], V[i])
        functions(problem).f(F[i], tqᵢ, Q[i], V[i])
    end

    # copy y to q and p
    for k in 1:D
        q[k] = y[(3*S+0)*D+k]
        p[k] = y[(3*S+1)*D+k]
    end

    # compute p=p(t,q)
    functions(problem).ϑ(ϕ, solstep.t̄ + timestep(problem), q, v)

    # for Lobatto-type methods, copy y to μ
    if isdefined(tableau(method), :d) && length(tableau(method).d) > 0
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
                y += tableau(method).q.a[i,j] * V[j][k]
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
                z += tableau(method).p.a[i,j] * F[j][k]
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
                ω = tableau(method).p.b[j] * tableau(method).p.c[j]^(i-1)
                # ω = tableau(method).q.a[i+1,j]
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
            y += tableau(method).q.b[j] * V[j][k]
        end
        b[offset+k] = (solstep.q̄[k] - q[k]) + timestep(problem) * y
    end

    # compute b = - (p̄-p-AF)
    offset = (3*S+1)*D
    for k in 1:D
        z = 0
        for j in 1:S
            z += tableau(method).p.b[j] * F[j][k]
        end
        b[offset+k] = (solstep.p̄[k] - p[k]) + timestep(problem) * z
    end

    # if isdefined(tableau(method), :d) && length(tableau(method).d) > 0
    #     offset = (3*S+2)*D
    #     for k in 1:D
    #         b[offset+k] = 0
    #         for i in 1:S
    #             b[offset+k] -= V[k,i] * tableau(method).d[i]
    #         end
    #     end
    # end

    if isdefined(tableau(method), :d) && length(tableau(method).d) > 0
        sl     = div(S+1, 2)
        offset = D*(S+sl-1)

        # compute μ
        z = tableau(method).p.b[sl] / tableau(method).d[sl]
        for k in 1:D
            μ[k] = z * b[offset+k]
        end

        # replace equation for Pₗ with constraint on V
        for k in 1:D
            z = 0
            for i in 1:S
                z += V[k,i] * tableau(method).d[i]
            end
            b[offset+k] = z
        end

        # modify P₁, ..., Pₛ except for Pₗ
        for i in 1:S
            if i ≠ sl
                offset = D*(S+i-1)
                z = tableau(method).d[i] / tableau(method).p.b[i]
                for k in 1:D
                    b[offset+k] -= z * μ[k]
                end
            end
        end
    end
end

"Compute stages of variational special partitioned additive Runge-Kutta methods."
function Integrators.residual!(y::Vector{ST}, b::Vector{ST},
                params::ParametersVPRKpLegendre{DT,TT,D,S},
                caches::OldCacheDict) where {ST,DT,TT,D,S}

    @assert length(y) == length(b)

    # get cache for internal stages
    cache = caches[ST]

    components!(y, cache.Q, cache.V, cache.P, cache.F, cache.Y, cache.Z, cache.Φ, cache.q̃, cache.ṽ, cache.p̃, cache.ϕ, cache.μ, params)

    compute_rhs!(b, cache.Q, cache.V, cache.P, cache.F, cache.Y, cache.Z, cache.Φ, cache.q̃, cache.p̃, cache.ϕ, cache.μ, params)

    # debug output
    # println()
    # for k in 1:D
    #     println(params.q[k], ",  ", params.p[k], ",  ", cache.Q[k,1], ",  ", cache.P[k,1], ",  ", cache.q[k], ",  ", cache.p[k])
    # end
    # println()
end


function integrate_step!(int::IntegratorVPRKpLegendre{DT,TT}, sol::SolutionStepPODE{DT,TT},
                                     cache::IntegratorCacheVPRK{DT}=int.caches[DT]) where {DT,TT}
    local offset::Int

    # update nonlinear solver parameters from cache
    update_params!(int.params, sol)

    # reset solution
    reset!(sol, timestep(int))

    # compute initial guess
    initial_guess!(int, sol, cache)

    # debug output
    # println()
    # println(cache.x)

    # call nonlinear solver
    solve!(cache.x, int.solver)

    # print solver status
    # print_solver_status(int.solver.status, int.solver.params)

    # check if solution contains NaNs or error bounds are violated
    # check_solver_status(int.solver.status, int.solver.params)

    # compute vector fields at internal stages
    components!(cache.x, cache.Q, cache.V, cache.P,
                    cache.F, cache.Y, cache.Z, cache.Φ,
                    cache.q̃, cache.ṽ, cache.p̃, cache.ϕ,
                    cache.μ, int.params)

    # compute final update
    update_solution!(sol.q, cache.V, tableau(int).q.b, tableau(int).q.b̂, timestep(int))
    update_solution!(sol.p, cache.F, tableau(int).p.b, tableau(int).p.b̂, timestep(int))
end
