
abstract type HPIMethod <: LODEMethod end

isiodemethod(::Union{HPIMethod, Type{<:HPIMethod}}) = true

default_solver(::HPIMethod) = Newton()
default_iguess(::HPIMethod) = HermiteExtrapolation()


function initial_guess!(int::GeometricIntegrator{<:HPIMethod})
    # set some local variables for convenience
    local D = ndims(int)
    local A = nparams(method(int))
    local x = nlsolution(int)

    # compute initial guess for solution q(n+1)
    initialguess!(solstep(int).t, cache(int).q, cache(int).p, solstep(int), problem(int), iguess(int))

    # copy initial guess to solution vector
    x[1:D] .= cache(int).q
    x[D+1:D+A] .= 0
end


function update!(x::AbstractVector{DT}, int::GeometricIntegrator{<:HPIMethod}) where {DT}
    # copy previous solution from solstep to cache
    reset!(cache(int, DT), current(solstep(int))...)

    # compute vector field at internal stages
    components!(x, int)

    # compute final update
    solstep(int).q .= cache(int, DT).q
    solstep(int).p .= cache(int, DT).p
end


function integrate_step!(int::GeometricIntegrator{<:HPIMethod, <:AbstractProblemIODE})
    # copy previous solution from solstep to cache
    reset!(cache(int), current(solstep(int))...)

    # call nonlinear solver
    solve!(nlsolution(int), (b,x) -> residual!(b, x, int), solver(int))

    # print solver status
    # print_solver_status(int.solver.status, int.solver.params)

    # check if solution contains NaNs or error bounds are violated
    # check_solver_status(int.solver.status, int.solver.params)

    # compute final update
    update!(nlsolution(int), int)
end
