

"""
Create nonlinear solver object for a system of `N` equations with data type `DT`.
The function ``f(x)=0`` to be solved for is determined by a julia function
`function_stages!(x, b, params)`, where `x` is the current solution and `b` is
the output vector, s.th. ``b = f(x)``. `params` are a set of parameters depending
on the equation and integrator that is used.
The solver type is obtained from the config dictionary (`:nls_solver`).
"""
function create_nonlinear_solver(DT, N, params)
    # create solution vector for nonlinear solver
    x = zeros(DT, N)

    # create wrapper function f(x,b) that calls `function_stages!(x, b, params)`
    # with the appropriate `params`
    f = (x,b) -> function_stages!(x, b, params)

    # create nonlinear solver with solver type obtained from config dictionary
    s = get_config(:nls_solver)(x, f)
end


"""
Create a solution vector of type `Double{DT}` for a problem with `D` dimensions
and `M` independent initial conditions.
"""
function create_solution_vector_double_double(DT, D, M)
    x = Array{Vector{Double{DT}}}(M)

    for i in 1:M
        x[i] = zeros(Double{DT}, D)
    end

    return x
end


function check_solution_dimension_asserts(sol::Solution, m::Int, n::Int)
    @assert m ≥ 1
    @assert m ≤ sol.ni

    @assert n ≥ 1
    @assert n ≤ sol.ntime
end
