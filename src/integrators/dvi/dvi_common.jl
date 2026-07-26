
abstract type DVIMethod <: LODEMethod end

isexplicit(::Union{DVIMethod,Type{<:DVIMethod}}) = false
isimplicit(::Union{DVIMethod,Type{<:DVIMethod}}) = true
issymplectic(::Union{DVIMethod,Type{<:DVIMethod}}) = true

default_solver(::DVIMethod) = Newton()
default_iguess(::DVIMethod) = HermiteExtrapolation()


"""
    check_dvi_dimension(D)

Throw an `ArgumentError` unless the configuration-space dimension `D` is even.

All degenerate variational integrators split the coordinates into two halves of
equal size: the first `D÷2` components carry the update rule, the second `D÷2` are
determined by the constraint `p = ϑ(q)`. For odd `D` there is no such splitting;
the `div(D, 2)` bounds in the residuals then drop equations — for `D = 1` the
update rules are lost entirely — and the nonlinear system degenerates. Without
this check the methods fail with a `SingularException` instead of a useful message.
"""
function check_dvi_dimension(D)
    iseven(D) || throw(ArgumentError(
        "Degenerate variational integrators require an even-dimensional " *
        "configuration space, got D = $(D)."))
    return nothing
end


function residual!(b::AbstractVector{ST}, x::AbstractVector{ST}, sol, params, int::GeometricIntegrator{<:DVIMethod,<:AbstractProblemIODE}) where {ST}
    # check that x and b are compatible
    @assert axes(x) == axes(b)

    # compute stages from nonlinear solver solution x
    components!(x, sol, params, int)

    # compute residual vector
    residual!(b, sol, params, int)
end


function update!(sol, params, x::AbstractVector{DT}, int::GeometricIntegrator{<:DVIMethod,<:AbstractProblemIODE}) where {DT}
    # compute vector field at internal stages
    components!(x, sol, params, int)

    # compute final update
    update!(sol, params, int, DT)
end


function integrate_step!(sol, history, params, int::GeometricIntegrator{<:DVIMethod,<:AbstractProblemIODE})
    # call nonlinear solver
    solve!(nlsolution(int), solver(int), solverstate(int), (sol, params, int))

    # print solver status
    # print_solver_status(int.solver.status, int.solver.params)

    # check if solution contains NaNs or error bounds are violated
    # check_solver_status(int.solver.status, int.solver.params)

    # compute final update
    update!(sol, params, nlsolution(int), int)
end
