
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


# Compute stages of Hamilton-Pontryagin integrators.
function residual!(b::AbstractVector{ST}, x::AbstractVector{ST}, sol, params, int::GeometricIntegrator{<:HPIMethod, <:AbstractProblemIODE}) where {ST}
    @assert axes(x) == axes(b)

    # copy previous solution from solstep to cache
    reset!(cache(int, ST), sol...)

    # compute stages from nonlinear solver solution x
    components!(x, sol, params, int)

    # compute residual vector
    residual!(b, int)
end


function update!(sol, params, x::AbstractVector{DT}, int::GeometricIntegrator{<:HPIMethod}) where {DT}
    # copy previous solution from solstep to cache
    reset!(cache(int, DT), sol...)

    # compute vector field at internal stages
    components!(x, sol, params, int)

    # compute final update
    sol.q .= cache(int, DT).q
    sol.p .= cache(int, DT).p
end


function integrate_step!(sol, history, params, int::GeometricIntegrator{<:HPIMethod, <:AbstractProblemIODE})
    # call nonlinear solver
    solve!(nlsolution(int), (b,x) -> residual!(b, x, sol, params, int), solver(int))

    # print solver status
    # print_solver_status(int.solver.status, int.solver.params)

    # check if solution contains NaNs or error bounds are violated
    # check_solver_status(int.solver.status, int.solver.params)

    # compute final update
    update!(sol, params, nlsolution(int), int)
end
