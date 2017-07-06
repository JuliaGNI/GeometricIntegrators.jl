
"Parameters for right-hand side function of fully implicit Runge-Kutta methods."
mutable struct NonlinearFunctionParametersFIRK{DT,TT,VT,D,S} <: NonlinearFunctionParameters{DT}
    v::VT
    Δt::TT

    a::Matrix{TT}
    c::Vector{TT}

    t::TT

    q::Vector{DT}
end

function NonlinearFunctionParametersFIRK(v::VT, Δt::TT, tab, q::Vector{DT}) where {DT,TT,VT}
    NonlinearFunctionParametersFIRK{DT,TT,VT,length(q),tab.s}(v, Δt, tab.a, tab.c, 0, q)
end

struct NonlinearFunctionCacheFIRK{DT}
    Q::Matrix{DT}
    V::Matrix{DT}
    Y::Matrix{DT}

    tQ::Vector{DT}
    tV::Vector{DT}

    function NonlinearFunctionCacheFIRK{DT}(d, s) where {DT}

        # create internal stage vectors
        Q = zeros(DT,d,s)
        V = zeros(DT,d,s)
        Y = zeros(DT,d,s)

        # create temporary vectors
        tQ = zeros(DT,d)
        tV = zeros(DT,d)

        new(Q, V, Y, tQ, tV)
    end
end

function compute_stages_firk!(x::Vector{DT}, Q::Matrix{DT}, V::Matrix{DT}, Y::Matrix{DT},
                              q::Vector, a::Matrix{TT}, c::Vector{TT}, Δt::TT, t::TT, v::VT,
                              tQ::Vector{DT}, tV::Vector{DT}) where {DT,TT,VT}

    local d::Int = length(q)
    local s::Int = length(c)
    local tᵢ::TT

    @assert d == length(tQ) == length(tV)
    @assert s == size(a,1) == size(a,2)
    @assert d == size(Q,1) == size(V,1) == size(Y,1)
    @assert s == size(Q,2) == size(V,2) == size(Y,2)

    # copy x to Y and compute Q = q + Δt Y
    for i in 1:size(Y,2)
        for k in 1:size(Y,1)
            Y[k,i] = x[d*(i-1)+k]
            Q[k,i] = q[k] + Δt * Y[k,i]
        end
    end

    # compute V = v(Q)
    for i in 1:s
        tᵢ = t + Δt * c[i]
        simd_copy_xy_first!(tQ, Q, i)
        v(tᵢ, tQ, tV)
        simd_copy_yx_first!(tV, V, i)
    end

    nothing
end

"Compute stages of fully implicit Runge-Kutta methods."
@generated function function_stages!(x::Vector{ST}, b::Vector{ST}, params::NonlinearFunctionParametersFIRK{DT,TT,VT,D,S}) where {ST,DT,TT,VT,D,S}

    cache = NonlinearFunctionCacheFIRK{ST}(D, S)

    function_stages = quote

        compute_stages_firk!(x, $cache.Q, $cache.V, $cache.Y, params.q, params.a, params.c, params.Δt, params.t, params.v, $cache.tQ, $cache.tV)

        # compute b = - (Y-AV)
        for i in 1:S
            for k in 1:D
                b[D*(i-1)+k] = - $cache.Y[k,i]
                for j in 1:S
                    b[D*(i-1)+k] += params.a[i,j] * $cache.V[k,j]
                end
            end
        end
    end

    return function_stages
end


"Fully implicit Runge-Kutta integrator."
struct IntegratorFIRK{DT, TT, FT, SPT, ST, IT} <: Integrator{DT,TT}
    equation::ODE{DT,TT,FT}
    tableau::TableauFIRK{TT}
    Δt::TT

    Sparams::SPT
    solver::ST
    iguess::InitialGuessODE{DT, TT, FT, IT}

    q::Vector{DT}
    v::Vector{DT}
    y::Vector{DT}

    Q::Matrix{DT}
    V::Matrix{DT}
    Y::Matrix{DT}

    tQ::Vector{DT}
    tV::Vector{DT}
end

function IntegratorFIRK(equation::ODE{DT,TT,FT}, tableau::TableauFIRK{TT}, Δt::TT;
                        interpolation=HermiteInterpolation{DT}) where {DT,TT,FT}
    D = equation.d
    S = tableau.q.s

    # create solution vector for internal stages / nonlinear solver
    x = zeros(DT, D*S)

    # create solution vectors
    q = zeros(DT,D)
    v = zeros(DT,D)

    # create compensated summation error vectors
    qₑᵣᵣ = zeros(DT,D)

    # create update vector
    y = zeros(DT,D)

    # create internal stage vectors
    Q = zeros(DT,D,S)
    V = zeros(DT,D,S)
    Y = zeros(DT,D,S)

    # create temporary vectors
    tQ = zeros(DT,D)
    tV = zeros(DT,D)

    # create params
    params = NonlinearFunctionParametersFIRK(equation.v, Δt, tableau.q, q)

    # create rhs function for nonlinear solver
    function_stages = (x,b) -> function_stages!(x, b, params)

    # create solver
    solver = get_config(:nls_solver)(x, function_stages)

    # create initial guess
    iguess = InitialGuessODE(interpolation, equation, Δt)

    # create integrator
    IntegratorFIRK{DT, TT, FT, typeof(params), typeof(solver), typeof(iguess.int)}(
                                        equation, tableau, Δt, params, solver, iguess,
                                        q, qₑᵣᵣ, v, y, Q, V, Y, tQ, tV)
end


function initialize!(int::IntegratorFIRK{DT,TT}, sol::SolutionODE, m::Int) where {DT,TT}
    @assert m ≥ 1
    @assert m ≤ sol.ni

    # copy initial conditions from solution
    get_initial_conditions!(sol, int.q, m)

    # initialise initial guess
    initialize!(int.iguess, sol.t[0], int.q)

    # reset compensated summation error
    int.qₑᵣᵣ .= 0
end

"Integrate ODE with fully implicit Runge-Kutta integrator."
function integrate_step!(int::IntegratorFIRK{DT, TT, FT, SPT, ST, IT}, sol::SolutionODE{DT,TT,N}, m::Int, n::Int) where {DT,TT,FT,SPT,ST,IT,N}
    # set time for nonlinear solver
    int.Sparams.t = sol.t[0] + (n-1)*int.Δt

    # copy previous solution to initial guess
    update!(int.iguess, sol.t[0] + (n-1)*int.Δt, int.q)

    # compute initial guess for internal stages
    for i in 1:int.tableau.q.s
        evaluate!(int.iguess, int.y, int.v, int.tableau.q.c[i])
        for k in 1:int.equation.d
            int.V[k,i] = int.v[k]
        end
    end
    for i in 1:int.tableau.q.s
        for k in 1:int.equation.d
            int.solver.x[int.equation.d*(i-1)+k] = 0
            for j in 1:int.tableau.q.s
                int.solver.x[int.equation.d*(i-1)+k] += int.tableau.q.a[i,j] * int.V[k,j]
            end
        end
    end

    # call nonlinear solver
    solve!(int.solver)

    if !solverStatusOK(int.solver.status, int.solver.params)
        println(int.solver.status)
    end

    compute_stages_firk!(int.solver.x, int.Q, int.V, int.Y, int.q,
                         int.tableau.q.a, int.tableau.q.c, int.Δt, sol.t[0] + (n-1)*int.Δt,
                         int.equation.v, int.tQ, int.tV)

    # compute final update
    for k in 1:size(int.V, 1)
        for i in 1:size(int.V, 2)
            int.q[k], int.qₑᵣᵣ[k] = compensated_summation(int.Δt * int.tableau.q.b[i] * int.V[k,i], int.q[k], int.qₑᵣᵣ[k])
            int.q[k], int.qₑᵣᵣ[k] = compensated_summation(int.Δt * int.tableau.q.b̂[i] * int.V[k,i], int.q[k], int.qₑᵣᵣ[k])
        end
    end

    # take care of periodic solutions
    for k in 1:int.equation.d
        if int.equation.periodicity[k] ≠ 0
            while int.q[k] < 0
                (int.q[k], int.qₑᵣᵣ[k]) = compensated_summation(+int.equation.periodicity[k], int.q[k], int.qₑᵣᵣ[k])
            end
            while int.q[k] ≥ int.equation.periodicity[k]
                (int.q[k], int.qₑᵣᵣ[k]) = compensated_summation(-int.equation.periodicity[k], int.q[k], int.qₑᵣᵣ[k])
            end
        end
    end

    # copy to solution
    copy_solution!(sol, int.q, n, m)
end

"Integrate partitioned ODE with fully implicit Runge-Kutta integrator."
function integrate!(int::IntegratorFIRK, s::SolutionPODE)
    # TODO
end
