
abstract type PMVIMethod <: VIMethod end

isexplicit(method::PMVIMethod) = false
isimplicit(method::PMVIMethod) = true
issymplectic(method::PMVIMethod) = true


function initial_guess!(int::GeometricIntegrator{<:PMVIMethod})
    # compute initial guess for solution q(n+1)
    initialguess!(solstep(int).t, cache(int).q, cache(int).p, solstep(int), problem(int), iguess(int))

    # copy initial guess to solution vector
    nlsolution(int) .= cache(int).q
end


function residual!(b::Vector{ST}, x::Vector{ST}, int::GeometricIntegrator{<:PMVIMethod}) where {ST}
    # copy previous solution from solstep to cache
    reset!(cache(int, ST), current(solstep(int))...)

    # compute stages from nonlinear solver solution x
    components!(x, int)

    # compute residual vector
    residual!(b, int)
end


function update!(x::AbstractVector{DT}, int::GeometricIntegrator{<:PMVIMethod}) where {DT}
    # copy previous solution from solstep to cache
    reset!(cache(int, DT), current(solstep(int))...)

    # compute vector field at internal stages
    components!(x, int)

    # compute final update
    solstep(int).q .= cache(int).q
    solstep(int).p .= cache(int).p
end


function integrate_step!(int::GeometricIntegrator{<:PMVIMethod, <:AbstractProblemIODE})
    # call nonlinear solver
    solve!(nlsolution(int), (b,x) -> residual!(b, x, int), solver(int))

    # print solver status
    # print_solver_status(int.solver.status, int.solver.params)

    # check if solution contains NaNs or error bounds are violated
    # check_solver_status(int.solver.status, int.solver.params)

    # compute final update
    update!(nlsolution(int), int)
end
