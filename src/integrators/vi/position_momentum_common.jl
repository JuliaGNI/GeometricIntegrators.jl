
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


function residual!(b::Vector{ST}, x::Vector{ST}, sol, params, int::GeometricIntegrator{<:PMVIMethod}) where {ST}
    # check that x and b are compatible
    @assert axes(x) == axes(b)

    # compute stages from nonlinear solver solution x
    components!(x, sol, params, int)

    # compute residual vector
    residual!(b, sol, params, int)
end


function update!(sol, params, int::GeometricIntegrator{<:PMVIMethod}, DT)
    # compute final update
    sol.q .= cache(int).q
    sol.p .= cache(int).p
end

function update!(sol, params, x::AbstractVector{DT}, int::GeometricIntegrator{<:PMVIMethod}) where {DT}
    # compute vector field at internal stages
    components!(x, sol, params, int)

    # compute final update
    update!(sol, params, int, DT)
end


function integrate_step!(sol, history, params, int::GeometricIntegrator{<:PMVIMethod, <:AbstractProblemIODE})
    # call nonlinear solver
    solve!(nlsolution(int), (b,x) -> residual!(b, x, sol, params, int), solver(int))

    # print solver status
    # print_solver_status(int.solver.status, int.solver.params)

    # check if solution contains NaNs or error bounds are violated
    # check_solver_status(int.solver.status, int.solver.params)

    # compute final update
    update!(sol, params, nlsolution(int), int)
end
