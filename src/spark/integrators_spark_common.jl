
function internal_variables(method::AbstractSPARKMethod, problem::AbstractSPARKProblem{DT,TT}) where {DT,TT}
    S = nstages(method)
    R = pstages(method)
    D = ndims(problem)

    Qi = create_internal_stage_vector(DT, D, S)
    Pi = create_internal_stage_vector(DT, D, S)
    Vi = create_internal_stage_vector(DT, D, S)
    Φi = create_internal_stage_vector(DT, D, S)

    Qp = create_internal_stage_vector(DT, D, R)
    Pp = create_internal_stage_vector(DT, D, R)
    Λp = create_internal_stage_vector(DT, D, R)
    Φp = create_internal_stage_vector(DT, D, R)

    # solver = get_solver_status(int.solver)

    (Qi=Qi, Pi=Pi, Vi=Vi, Φi=Φi, Qp=Qp, Pp=Pp, Λp=Λp, Φp=Φp)#, solver=solver)
end


function copy_internal_variables!(solstep::SolutionStepPDAE, cache::IntegratorCacheSPARK)
    haskey(internal(solstep), :Qi) && copyto!(internal(solstep).Qi, cache.Qi)
    haskey(internal(solstep), :Pi) && copyto!(internal(solstep).Pi, cache.Pi)
    haskey(internal(solstep), :Vi) && copyto!(internal(solstep).Vi, cache.Vi)
    haskey(internal(solstep), :Φi) && copyto!(internal(solstep).Φi, cache.Φi)

    haskey(internal(solstep), :Qp) && copyto!(internal(solstep).Qp, cache.Qp)
    haskey(internal(solstep), :Pp) && copyto!(internal(solstep).Pp, cache.Pp)
    haskey(internal(solstep), :Vp) && copyto!(internal(solstep).Λp, cache.Λp)
    haskey(internal(solstep), :Φp) && copyto!(internal(solstep).Φp, cache.Φp)
end


function update!(sol, params, x::AbstractVector{DT}, int::GeometricIntegrator{<:AbstractSPARKMethod}) where {DT}
    # compute vector field at internal stages
    components!(x, sol, params, int)

    # compute final update and projection
    update!(sol.q, cache(int, DT).Vi, cache(int, DT).Up, tableau(int).q, timestep(int))
    update!(sol.p, cache(int, DT).Fi, cache(int, DT).Gp, tableau(int).p, timestep(int))
    update!(sol.λ, cache(int, DT).Λp, tableau(int).λ, timestep(int))

    return sol
end


function integrate_step!(sol, history, params, int::GeometricIntegrator{<:AbstractSPARKMethod,<:AbstractSPARKProblem})
    # call nonlinear solver
    solve!(solver(int), nlsolution(int), (sol, params, int))

    # check_jacobian(int.solver)
    # print_jacobian(int.solver)

    # print solver status
    # println(status(solver))

    # check if solution contains NaNs or error bounds are violated
    # println(meets_stopping_criteria(status(solver)))

    # compute final update
    update!(sol, params, nlsolution(int), int)
end
