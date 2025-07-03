
abstract type PMVIMethod <: VIMethod end

isexplicit(method::PMVIMethod) = false
isimplicit(method::PMVIMethod) = true
issymplectic(method::PMVIMethod) = true


function initial_guess!(sol, history, params, int::GeometricIntegrator{<:PMVIMethod})
    # compute initial guess for solution q(n+1)
    soltmp = (
        t=sol.t,
        q=cache(int).q̃,
        p=cache(int).θ̃,
        v=cache(int).ṽ,
        f=cache(int).f̃,
    )
    solutionstep!(soltmp, history, problem(int), iguess(int))

    # copy initial guess to solution vector
    nlsolution(int) .= cache(int).q̃
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


function integrate_step!(sol, history, params, int::GeometricIntegrator{<:PMVIMethod,<:AbstractProblemIODE})
    # call nonlinear solver
    solve!(solver(int), nlsolution(int), (sol, params, int))

    # print solver status
    # print_solver_status(int.solver.status, int.solver.params)

    # check if solution contains NaNs or error bounds are violated
    # check_solver_status(int.solver.status, int.solver.params)

    # compute final update
    update!(sol, params, nlsolution(int), int)
end
