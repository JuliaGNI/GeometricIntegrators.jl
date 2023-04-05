
"Parameters for right-hand side function of Variational Partitioned Runge-Kutta methods."
const ParametersVPRKpStandard = AbstractParametersVPRK{:vprk_pstandard}

# function IntegratorCache(params::ParametersVPRKpStandard{DT,TT,D,S}; kwargs...) where {DT,TT,D,S}
#     IntegratorCacheVPRK{DT,D,S}(D*S, true, 2*D; kwargs...)
# end

# function IntegratorCache{ST}(params::ParametersVPRKpStandard{DT,TT,D,S}; kwargs...) where {ST,DT,TT,D,S}
#     IntegratorCacheVPRK{ST,D,S}(D*S, true, 2*D; kwargs...)
# end




# "Variational partitioned Runge-Kutta integrator with standard projection."
# const IntegratorVPRKpStandard{DT,TT} = Integrator{<:Union{IODEProblem{DT,TT},LODEProblem{DT,TT}}, <:ProjectedVPRK{<:VPRKMethod, <:VariationalProjection}}

# "Variational partitioned Runge-Kutta integrator with symplectic projection."
# const IntegratorVPRKpSymplectic{DT,TT} = Integrator{<:Union{IODEProblem{DT,TT},LODEProblem{DT,TT}}, <:ProjectedVPRK{<:VPRKMethod, <:SymplecticProjection}}

# @doc raw"""
# Variational partitioned Runge-Kutta integrator with variational projection on $(q_{n}, p_{n+1})$.
# """
# const IntegratorVPRKpVariationalQ{DT,TT} = Integrator{<:Union{IODEProblem{DT,TT},LODEProblem{DT,TT}}, <:ProjectedVPRK{<:VPRKMethod, <:VariationalProjectionOnQ}}

# @doc raw"""
# Variational partitioned Runge-Kutta integrator with variational projection on $(p_{n}, q_{n+1})$.
# """
# const IntegratorVPRKpVariationalP{DT,TT} = Integrator{<:Union{IODEProblem{DT,TT},LODEProblem{DT,TT}}, <:ProjectedVPRK{<:VPRKMethod, <:VariationalProjectionOnP}}


description(::IntegratorVPRKpStandard) = "Variational Partitioned Runge-Kutta Integrator with standard projection"
description(::IntegratorVPRKpSymplectic) = "Variational Partitioned Runge-Kutta Integrator with symplectic projection"
description(::IntegratorVPRKpVariationalQ) = @doc raw"Variational Partitioned Runge-Kutta Integrator with variational projection on $(q_{n}, p_{n+1})$"
description(::IntegratorVPRKpVariationalP) = @doc raw"Variational Partitioned Runge-Kutta Integrator with variational projection on $(p_{n}, q_{n+1})$"


# function Integrators.get_internal_variables(int::IntegratorVPRKpStandard{DT,TT,D,S}) where {DT, TT, D, S}
#     Q = create_internal_stage_vector(DT, D, S)
#     P = create_internal_stage_vector(DT, D, S)
#     V = create_internal_stage_vector(DT, D, S)
#     F = create_internal_stage_vector(DT, D, S)
#     λ = zeros(DT,D)

#     solver = get_solver_status(int.solver)

#     (Q=Q, P=P, V=V, F=F, λ=λ, solver=solver)
# end








const AbstractVPRKpStandard = Union{
    ProjectedVPRK{<:VPRKMethod, <:VariationalProjection},
    ProjectedVPRK{<:VPRKMethod, <:SymplecticProjection},
    ProjectedVPRK{<:VPRKMethod, <:VariationalProjectionOnQ},
    ProjectedVPRK{<:VPRKMethod, <:VariationalProjectionOnP}
}


solversize(problem::VPRKProblem, ::AbstractVPRKpStandard) = 2 * ndims(problem)



function Integrators.initialize!(int::IntegratorVPRKpStandard{DT}, sol::SolutionStepPODE{DT},
                                 cache::IntegratorCacheVPRK{DT}=int.caches[DT]) where {DT}
    # initialise projector
    cache.U[1] .= cache.λ
    equation(int, :g)(cache.G[1], sol.t, sol.q, sol.v, cache.λ)
end


function initial_guess_projection!(int::IntegratorVPRKpStandard{DT}, sol::SolutionStepPODE{DT},
                                   cache::IntegratorCacheVPRK{DT}=int.caches[DT]) where {DT}
    offset_q = 0
    offset_λ = ndims(int)

    for k in eachdim(int)
        cache.x̄[offset_q+k] = sol.q[k]
        cache.x̄[offset_λ+k] = 0
    end
end


function compute_projection!(x::Vector{ST}, 
                q::SolutionVector{ST}, p::SolutionVector{ST},
                v::SolutionVector{ST}, λ::SolutionVector{ST},
                U::Vector{Vector{ST}}, G::Vector{Vector{ST}},
                params::ParametersVPRKpStandard{DT,TT,D,S}) where {ST,DT,TT,D,S}

    @assert D == length(q) == length(p) == length(v) == length(λ)
    @assert D == length(U[1]) == length(U[2])
    @assert D == length(G[1]) == length(G[2])

    # copy x to q, λ
    for k in 1:D
        q[k] = x[0*D+k]
        λ[k] = x[1*D+k]
    end

    # compute u=λ and g=∇ϑ(q)⋅λ
    U[1] .= λ
    U[2] .= λ

    functions(problem).g(G[1], solstep.t̄[1], q, v, λ)
    G[2] .= G[1]

    # compute p̄=ϑ(q)
    functions(problem).ϑ(p, solstep.t̄[1], q, v)
end

"Compute stages of projected variational partitioned Runge-Kutta methods."
function Integrators.function_stages!(x::Vector{ST}, b::Vector{ST},
                params::ParametersVPRKpStandard{DT,TT,D,S},
                caches::OldCacheDict) where {ST,DT,TT,D,S}

    @assert length(x) == length(b)

    # get cache for internal stages
    cache = caches[ST]

    compute_projection!(x, cache.q̃, cache.p̃, cache.ṽ, cache.λ, cache.U, cache.G, params)

    # compute b = - [q̄-q-U]
    for k in 1:D
        b[0*D+k] = - (cache.q̃[k] - solstep.q̄[1][k]) + timestep(problem) * params.pparams[:RU][2] * cache.U[2][k]
    end

    # compute b = - [p̄-p-G]
    for k in 1:D
        b[1*D+k] = - (cache.p̃[k] - solstep.p̄[1][k]) + timestep(problem) * params.pparams[:RG][2] * cache.G[2][k]
    end
end


function integrate_step!(int::IntegratorVPRKpStandard{DT,TT}, sol::SolutionStepPODE{DT,TT},
                                     cache::IntegratorCacheVPRK{DT}=int.caches[DT]) where {DT,TT}
    # add perturbation for next time step to solution
    # (same vector field as previous time step)
    project_solution!(int, sol, int.pparams.pparams[:RU1], int.pparams.pparams[:RG1], cache)

    # update nonlinear solver parameters from cache
    update_params!(int.params, sol)

    # reset cache
    reset!(sol, timestep(int))

    # compute initial guess
    initial_guess!(int, sol, cache)

    # call nonlinear solver
    solve!(cache.x, int.solver)

    # print solver status
    # print_solver_status(int.solver.status, int.solver.params)

    # check if solution contains NaNs or error bounds are violated
    # check_solver_status(int.solver.status, int.solver.params)

    # compute vector fields at internal stages
    compute_stages!(cache.x, cache.Q, cache.V, cache.P, cache.F, int.params)

    # compute unprojected solution
    update_solution!(int, sol, cache)

    # set time and solution for projection solver
    update_params!(int.pparams, sol)

    # set initial guess for projection
    initial_guess_projection!(int, sol, cache)

    # call projection solver
    solve!(cache.x̄, int.projector)

    # print solver status
    # print_solver_status(int.projector.status, int.projector.params)

    # check if solution contains NaNs or error bounds are violated
    # check_solver_status(int.projector.status, int.projector.params)

    # compute projection vector fields
    compute_projection!(cache.x̄, cache.q̃, cache.p̃, cache.ṽ, cache.λ, cache.U, cache.G, int.pparams)

    # add projection to solution
    project_solution!(int, sol, int.pparams.pparams[:RU2], int.pparams.pparams[:RG2], cache)

    # copy solution to initial guess
    update_vector_fields!(int.iguess, sol.t, sol.q, sol.p, sol.v, sol.f)

    # copy internal stage variables
    # sol.internal.Q .= cache.Q
    # sol.internal.P .= cache.P
    # sol.internal.V .= cache.V
    # sol.internal.F .= cache.F
    # sol.internal.λ .= cache.λ

    # copy solver status
    # get_solver_status!(int.solver, sol.internal[:solver])
end
