
abstract type DVIMethod <: LODEMethod end

isexplicit(::Union{DVIMethod, Type{<:DVIMethod}}) = false
isimplicit(::Union{DVIMethod, Type{<:DVIMethod}}) = true
issymplectic(::Union{DVIMethod, Type{<:DVIMethod}}) = true

default_solver(::DVIMethod) = Newton()
default_iguess(::DVIMethod) = HermiteExtrapolation()


function integrate_step!(int::GeometricIntegrator{<:DVIMethod, <:AbstractProblemIODE})
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
