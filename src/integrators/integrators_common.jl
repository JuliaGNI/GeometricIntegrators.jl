
# general helper functions for integrators.

_default_solver() = SimpleSolvers.NewtonSolver

"""
Create nonlinear solver object for a system of `N` equations with data type `DT`.
The function ``f(x)=0`` to be solved for is determined by a julia function
`function_stages!(x, b, params)`, where `x` is the current solution and `b` is
the output vector, s.th. ``b = f(x)``. `params` are a set of parameters depending
on the equation and integrator that is used.
The solver type is obtained from the config dictionary (`:nls_solver`).
"""
function create_nonlinear_solver(DT, N, params, caches, solver::Type{<:NonlinearSolver}=_default_solver(), F::Function=function_stages!)
    # create solution vector for nonlinear solver
    x = zeros(DT, N)

    # create wrapper function f!(x,b) that calls `F(x, b, params)`
    # with the appropriate `params`
    f! = (b,x) -> F(x, b, params, caches)

    # create nonlinear solver with solver type obtained from config dictionary
    solver(x, zero(x), f!; linesearch = Backtracking(), config = Options(min_iterations = 1, x_abstol = 8eps(), f_abstol = 8eps()))
end

function create_nonlinear_solver(DT, N, params, caches, i::Int, solver::Type{<:NonlinearSolver}=_default_solver(), F::Function=function_stages!)
    # create solution vector for nonlinear solver
    x = zeros(DT, N)

    # create wrapper function f(x,b) that calls `F(x, b, params)`
    # with the appropriate `params`
    f! = (b,x) -> F(x, b, params, caches, i)

    # create nonlinear solver with solver type obtained from config dictionary
    solver(x, zero(x), f!; linesearch = Backtracking(), config = Options(min_iterations = 1, x_abstol = 8eps(), f_abstol = 8eps()))
end

function create_nonlinear_solver_with_jacobian(DT, N, params, caches, solver::Type{<:NonlinearSolver}=_default_solver(), F=function_stages!, J=jacobian!)
    # create solution vector for nonlinear solver
    x = zeros(DT, N)

    # create wrapper function f!(x,b) that calls `F(x, b, params)`
    # with the appropriate `params`
    f! = (b,x) -> F(x, b, params, caches)

    # create wrapper function j!(x,df) that calls `J(x, df, params)`
    # with the appropriate `params`
    j! = (df,x) -> J(x, df, caches[DT], params)

    # create nonlinear solver with solver type obtained from config dictionary
    solver(x, zero(x), f!; J! = j!, linesearch = Backtracking(), config = Options(min_iterations = 1, x_abstol = 8eps(), f_abstol = 8eps()))
end


function cut_periodic_solution!(x::SolutionVector{T}, periodicity::AbstractVector{T}) where {T}
    @assert length(x) == length(periodicity)

    for k in eachindex(x, periodicity)
        if periodicity[k] ≠ 0
            while x[k] < 0
                x[k] += periodicity[k]
            end
            while x[k] ≥ periodicity[k]
                x[k] -= periodicity[k]
            end
        end
    end
end

function cut_periodic_solution!(x::SolutionVector{T}, periodicity::AbstractVector{T}, shift::AbstractVector{T}) where {T}
    @assert length(x) == length(periodicity)
    shift .= 0
    for k in eachindex(x, periodicity, shift)
        if periodicity[k] ≠ 0
            while x[k] + shift[k] < 0
                shift[k] += periodicity[k]
            end
            while x[k] + shift[k] ≥ periodicity[k]
                shift[k] -= periodicity[k]
            end
        end
    end
end
