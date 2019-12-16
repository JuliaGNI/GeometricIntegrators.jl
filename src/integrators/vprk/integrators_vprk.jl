@doc raw"""
`TableauVPRK`: Tableau of a Variational Partitioned Runge-Kutta method
```math
\begin{align*}
P_{n,i} &= \dfrac{\partial L}{\partial v} (Q_{n,i}, V_{n,i}) , &
Q_{n,i} &= q_{n} + h \sum \limits_{j=1}^{s} a_{ij} \, V_{n,j} , &
q_{n+1} &= q_{n} + h \sum \limits_{i=1}^{s} b_{i} \, V_{n,i} , \\
F_{n,i} &= \dfrac{\partial L}{\partial q} (Q_{n,i}, V_{n,i}) , &
P_{n,i} &= p_{n} + h \sum \limits_{i=1}^{s} \bar{a}_{ij} \, F_{n,j} - d_i \lambda , &
p_{n+1} &= p_{n} + h \sum \limits_{i=1}^{s} \bar{b}_{i} \, F_{n,i} , \\
&&
0 &= \sum \limits_{i=1}^{s} d_i V_i , &&
\end{align*}
```
satisfying the symplecticity conditions
```math
\begin{align*}
b_{i} \bar{a}_{ij} + b_{j} a_{ji} &= b_{i} b_{j} , &
\bar{b}_i &= b_i .
\end{align*}
```
"""
struct TableauVPRK{T} <: AbstractTableauPRK{T}
    @HeaderTableau

    q::CoefficientsRK{T}
    p::CoefficientsRK{T}

    R∞::Int

    d::Vector{T}

    function TableauVPRK{T}(name, o, q, p, R∞, d) where {T}
        @assert q.s == p.s == length(d)
        new(name, o, q.s, q, p, R∞, d)
    end

    function TableauVPRK{T}(name, o, q, p, R∞) where {T}
        @assert q.s == p.s
        new(name, o, q.s, q, p, R∞)
    end
end

function TableauVPRK(name::Symbol, order::Int, q::CoefficientsRK{T}, p::CoefficientsRK{T}, R∞::Int, d::Vector{T}) where {T}
    TableauVPRK{T}(name, order, q, p, R∞, d)
end

function TableauVPRK(name::Symbol, order::Int, q::CoefficientsRK{T}, p::CoefficientsRK{T}, R∞::Int) where {T}
    TableauVPRK{T}(name, order, q, p, R∞)
end

function TableauVPRK(name::Symbol, order::Int, q::CoefficientsRK{T}, R∞::Int, d::Vector{T}) where {T}
    TableauVPRK{T}(name, order, q, get_symplectic_conjugate_coefficients(q), R∞, d)
end

function TableauVPRK(name::Symbol, order::Int, q::CoefficientsRK{T}, R∞::Int) where {T}
    TableauVPRK{T}(name, order, q, get_symplectic_conjugate_coefficients(q), R∞)
end


"Parameters for right-hand side function of variational partitioned Runge-Kutta methods."
mutable struct ParametersVPRK{DT, TT, ET <: IODE{DT,TT}, D, S} <: AbstractParametersVPRK{DT,TT,ET,D,S}
    equ::ET
    tab::TableauVPRK{TT}
    Δt::TT

    t̅::TT
    q̅::Vector{DT}
    p̅::Vector{DT}
end

function ParametersVPRK(equ::ET, tab::TableauVPRK{TT}, Δt::TT) where {DT, TT, ET <: IODE{DT,TT}}
    q = zeros(DT, equ.d)
    p = zeros(DT, equ.d)

    ParametersVPRK{DT, TT, ET, equ.d, tab.s}(equ, tab, Δt, 0, q, p)
end


"Compute stages of variational partitioned Runge-Kutta methods."
@generated function function_stages!(x::Vector{ST}, b::Vector{ST},
                params::ParametersVPRK{DT,TT,ET,D,S}
            ) where {ST,DT,TT,ET,D,S}

    cache = IntegratorCacheVPRK{ST, D, S}()

    quote
        @assert length(x) == length(b)

        compute_stages!(x, $cache.Q, $cache.V, $cache.P, $cache.F, params)
        compute_rhs_vprk!(b, $cache.P, $cache.F, params)
        compute_rhs_vprk_correction!(b, $cache.V, params)
    end
end


"Variational partitioned Runge-Kutta integrator."
struct IntegratorVPRK{DT, TT, PT <: ParametersVPRK{DT,TT},
                              ST <: NonlinearSolver{DT},
                              IT <: InitialGuessPODE{DT,TT}, D, S} <: AbstractIntegratorVPRK{DT,TT}
    params::PT
    solver::ST
    iguess::IT
    cache::IntegratorCacheVPRK{DT,D,S}
end

function IntegratorVPRK(equation::ET, tableau::TableauVPRK{TT}, Δt::TT) where {DT, TT, ET <: IODE{DT,TT}}
    D = equation.d
    S = tableau.s

    # create params
    params = ParametersVPRK(equation, tableau, Δt)

    # create nonlinear solver
    solver = create_nonlinear_solver(DT, D*S, params)

    # create initial guess
    iguess = InitialGuessPODE(get_config(:ig_interpolation), equation, Δt)

    # create cache
    cache = IntegratorCacheVPRK{DT,D,S}()

    # create integrator
    IntegratorVPRK{DT, TT, typeof(params), typeof(solver), typeof(iguess), D, S}(params, solver, iguess, cache)
end


IntegratorVPRKpNone = IntegratorVPRK


equation(int::IntegratorVPRK) = int.params.equ
timestep(int::IntegratorVPRK) = int.params.Δt
tableau(integrator::IntegratorVPRK) = integrator.params.tab
nstages(integrator::IntegratorVPRK) = integrator.params.tab.s


function initial_guess!(int::IntegratorVPRK, sol::AtomisticSolutionPODE)
    for i in eachstage(int)
        evaluate!(int.iguess, sol.q, sol.p, sol.v, sol.f,
                              sol.q̅, sol.p̅, sol.v̅, sol.f̅,
                              int.cache.q̃, int.cache.ṽ,
                              tableau(int).q.c[i])

        for k in eachdim(int)
            int.solver.x[ndims(int)*(i-1)+k] = int.cache.ṽ[k]
        end
    end
end


"Integrate ODE with variational partitioned Runge-Kutta integrator."
function integrate_step!(int::IntegratorVPRK{DT,TT}, sol::AtomisticSolutionPODE{DT,TT}) where {DT,TT}
    # update nonlinear solver parameters from cache
    update_params!(int.params, sol)

    # compute initial guess
    initial_guess!(int, sol)

    # reset solution
    reset!(sol, timestep(int))

    # call nonlinear solver
    solve!(int.solver)

    # print solver status
    print_solver_status(int.solver.status, int.solver.params)

    # check if solution contains NaNs or error bounds are violated
    check_solver_status(int.solver.status, int.solver.params)

    # compute vector fields at internal stages
    compute_stages!(int.solver.x, int.cache.Q, int.cache.V, int.cache.P, int.cache.F, int.params)

    # compute final update
    update_solution!(int, sol)

    # copy solution to initial guess
    update!(int.iguess, sol.t, sol.q, sol.p, sol.v, sol.f)
end
