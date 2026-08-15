
abstract type HPIMethod <: LODEMethod end

isiodemethod(::Union{HPIMethod,Type{<:HPIMethod}}) = true

default_solver(::HPIMethod) = Newton()
default_iguess(::HPIMethod) = HermiteExtrapolation()


function initial_guess!(sol, history, params, int::GeometricIntegrator{<:HPIMethod})
    # set some local variables for convenience
    local D = length(cache(int).q̃)
    local A = nparams(method(int))
    local x = nlsolution(int)

    # compute initial guess for solution q(n+1)
    soltmp = (
        t=sol.t,
        q=cache(int).q̃,
        p=cache(int).θ̃,
        q̇=cache(int).ṽ,
        ṗ=cache(int).f̃,
    )
    solutionstep!(soltmp, history, problem(int), iguess(int))

    # copy initial guess to solution vector
    x[1:D] .= cache(int).q̃
    x[D+1:D+A] .= method(int).params
end


function update!(sol, params, x::AbstractVector{DT}, int::GeometricIntegrator{<:HPIMethod}) where {DT}
    # compute vector field at internal stages
    components!(x, sol, params, int)

    # compute final update
    sol.q .= cache(int, DT).q
    sol.p .= cache(int, DT).p
end


function integrate_step!(sol, history, params, int::GeometricIntegrator{<:HPIMethod,<:AbstractProblemIODE})
    # call nonlinear solver and act on the outcome it reports
    solverstatus = solve_with_status!(nlsolution(int), solver(int), solverstate(int), (sol, params, int))
    check_solver_status(solverstatus, int)

    # compute final update
    update!(sol, params, nlsolution(int), int)
end
