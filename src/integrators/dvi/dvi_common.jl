
abstract type DVIMethod <: LODEMethod end

isexplicit(::Union{DVIMethod, Type{<:DVIMethod}}) = false
isimplicit(::Union{DVIMethod, Type{<:DVIMethod}}) = true
issymplectic(::Union{DVIMethod, Type{<:DVIMethod}}) = true

default_solver(::DVIMethod) = Newton()
default_iguess(::DVIMethod) = HermiteExtrapolation()


function residual!(b::AbstractVector{ST}, x::AbstractVector{ST}, sol, params, int::GeometricIntegrator{<:DVIMethod, <:AbstractProblemIODE}) where {ST}
    # check that x and b are compatible
    @assert axes(x) == axes(b)

    # compute stages from nonlinear solver solution x
    components!(x, sol, params, int)

    # compute residual vector
    residual!(b, sol, params, int)
end


function update!(sol, params, x::AbstractVector{DT}, int::GeometricIntegrator{<:DVIMethod, <:AbstractProblemIODE}) where {DT}
    # compute vector field at internal stages
    components!(x, sol, params, int)

    # compute final update
    update!(sol, params, int, DT)
end


function integrate_step!(sol, history, params, int::GeometricIntegrator{<:DVIMethod, <:AbstractProblemIODE})
    # call nonlinear solver
    solve!(nlsolution(int), (b,x) -> residual!(b, x, sol, params, int), solver(int))

    # print solver status
    # print_solver_status(int.solver.status, int.solver.params)

    # check if solution contains NaNs or error bounds are violated
    # check_solver_status(int.solver.status, int.solver.params)

    # compute final update
    update!(sol, params, nlsolution(int), int)
end
