
using ForwardDiff

"Parameters for right-hand side function of fully implicit Runge-Kutta methods."
mutable struct ParametersFIRK{DT, TT, D, S, ET <: NamedTuple, FT, JT} <: Parameters{DT,TT}
    equs::ET
    tab::Tableau{TT}
    Δt::TT

    F::FT
    Jconfig::JT

    t::TT
    q::Vector{DT}

    function ParametersFIRK{DT,D}(equs::ET, tab::Tableau{TT}, Δt::TT) where {DT, TT, D, ET <: NamedTuple}
        F = (v,q) -> equs[:v](zero(TT), q, v)
        tq = zeros(DT,D)
        tv = zeros(DT,D)
        Jconfig = ForwardDiff.JacobianConfig(F, tv, tq)

        new{DT, TT, D, tab.s, ET, typeof(F), typeof(Jconfig)}(equs, tab, Δt, F, Jconfig, zero(TT), zeros(DT,D))
    end
end


@doc raw"""
Fully implicit Runge-Kutta integrator cache.

### Fields

* `q̃`: initial guess of solution
* `ṽ`: initial guess of vector field
* `s̃`: holds shift due to periodicity of solution
* `Q`: internal stages of solution
* `V`: internal stages of vector field
* `Y`: vector field of internal stages
"""
struct IntegratorCacheFIRK{DT,D,S} <: ODEIntegratorCache{DT,D}
    q̃::Vector{DT}
    ṽ::Vector{DT}
    s̃::Vector{DT}

    Q::Vector{Vector{DT}}
    V::Vector{Vector{DT}}
    Y::Vector{Vector{DT}}
    J::Vector{Matrix{DT}}

    function IntegratorCacheFIRK{DT,D,S}() where {DT,D,S}
        Q = create_internal_stage_vector(DT, D, S)
        V = create_internal_stage_vector(DT, D, S)
        Y = create_internal_stage_vector(DT, D, S)
        J = [zeros(DT,D,D) for i in 1:S]
        new(zeros(DT,D), zeros(DT,D), zeros(DT,D), Q, V, Y, J)
    end
end

function IntegratorCache{ST}(params::ParametersFIRK{DT,TT,D,S}; kwargs...) where {ST,DT,TT,D,S}
    IntegratorCacheFIRK{ST,D,S}(; kwargs...)
end

@inline CacheType(ST, params::ParametersFIRK{DT,TT,D,S}) where {DT,TT,D,S} = IntegratorCacheFIRK{ST,D,S}


@doc raw"""
Fully implicit Runge-Kutta integrator solving the system
```math
\begin{aligned}
V_{n,i} &= v (Q_{n,i}, P_{n,i}) , &
Q_{n,i} &= q_{n} + h \sum \limits_{j=1}^{s} a_{ij} \, V_{n,j} , &
q_{n+1} &= q_{n} + h \sum \limits_{i=1}^{s} b_{i} \, V_{n,i} .
\end{aligned}
```
"""
struct IntegratorFIRK{DT, TT, D, S, PT <: ParametersFIRK{DT,TT},
                                    ST <: NonlinearSolver{DT},
                                    IT <: InitialGuessODE{TT}} <: AbstractIntegratorRK{DT,TT}
    params::PT
    solver::ST
    iguess::IT
    caches::CacheDict{PT}

    function IntegratorFIRK(params::ParametersFIRK{DT,TT,D,S}, solver::ST, iguess::IT, caches) where {DT,TT,D,S,ST,IT}
        new{DT, TT, D, S, typeof(params), ST, IT}(params, solver, iguess, caches)
    end

    function IntegratorFIRK{DT,D}(equations::NamedTuple, tableau::Tableau{TT}, Δt::TT; exact_jacobian=true) where {DT,TT,D}
        # get number of stages
        S = tableau.s

        # check if tableau is fully implicit
        if get_config(:verbosity) ≥ 1
            if isexplicit(tableau)
                @warn "Initializing IntegratorFIRK with explicit tableau $(tableau.name).\nYou might want to use IntegratorERK instead."
            elseif isdiagnonallyimplicit(tableau)
                @warn "Initializing IntegratorFIRK with diagonally implicit tableau $(tableau.name).\nYou might want to use IntegratorDIRK instead."
            end
        end

        # create params
        params = ParametersFIRK{DT,D}(equations, tableau, Δt)

        # create cache dict
        caches = CacheDict(params)

        # create solver
        if exact_jacobian
            solver = create_nonlinear_solver_with_jacobian(DT, D*S, params, caches)
        else
            solver = create_nonlinear_solver(DT, D*S, params, caches)
        end

        # create initial guess
        iguess = InitialGuessODE(get_config(:ig_extrapolation), equations[:v], Δt)

        # create integrator
        IntegratorFIRK(params, solver, iguess, caches)
    end

    function IntegratorFIRK{DT,D}(v::Function, tableau::Tableau{TT}, Δt::TT; kwargs...) where {DT,TT,D}
        IntegratorFIRK{DT,D}(NamedTuple{(:v,)}((v,)), tableau, Δt; kwargs...)
    end

    function IntegratorFIRK{DT,D}(v::Function, h::Function, tableau::Tableau{TT}, Δt::TT; kwargs...) where {DT,TT,D}
        IntegratorFIRK{DT,D}(NamedTuple{(:v,:h)}((v,h)), tableau, Δt; kwargs...)
    end

    function IntegratorFIRK(equation::ODE{DT}, tableau::Tableau{TT}, Δt::TT; kwargs...) where {DT,TT}
        IntegratorFIRK{DT, ndims(equation)}(get_functions(equation), tableau, Δt; kwargs...)
    end
end


@inline Base.ndims(::IntegratorFIRK{DT,TT,D,S}) where {DT,TT,D,S} = D


function initialize!(int::IntegratorFIRK, sol::AtomicSolutionODE)
    sol.t̄ = sol.t - timestep(int)

    equations(int)[:v](sol.t, sol.q, sol.v)

    initialize!(int.iguess, sol.t, sol.q, sol.v,
                            sol.t̄, sol.q̄, sol.v̄)
end


function update_params!(int::IntegratorFIRK, sol::AtomicSolutionODE)
    # set time for nonlinear solver and copy previous solution
    int.params.t  = sol.t
    int.params.q .= sol.q
end


function initial_guess!(int::IntegratorFIRK{DT}, sol::AtomicSolutionODE{DT},
                        cache::IntegratorCacheFIRK{DT}=int.caches[DT]) where {DT}

    local offset::Int

    # compute initial guess for internal stages
    for i in eachstage(int)
        evaluate!(int.iguess, sol.q̄, sol.v̄, sol.q, sol.v, cache.Q[i], cache.V[i], tableau(int).c[i])
    end
    for i in eachstage(int)
        offset = ndims(int)*(i-1)
        for k in eachdim(int)
            int.solver.x[offset+k] = 0
            for j in eachstage(int)
                int.solver.x[offset+k] += timestep(int) * tableau(int).a[i,j] * cache.V[j][k]
            end
        end
    end
end


function compute_stages!(x::Vector{ST}, Q::Vector{Vector{ST}}, V::Vector{Vector{ST}}, Y::Vector{Vector{ST}},
                         params::ParametersFIRK{DT,TT,D}) where {ST,DT,TT,D}

    local tᵢ::TT

    # copy x to Y and compute Q = q + Δt Y
    for i in eachindex(Q,Y)
        for k in eachindex(Q[i],Y[i])
            Y[i][k] = x[D*(i-1)+k]
            Q[i][k] = params.q[k] + Y[i][k]
        end
    end

    # compute V = v(Q)
    for i in eachindex(Q,V)
        tᵢ = params.t + params.Δt * params.tab.c[i]
        params.equs[:v](tᵢ, Q[i], V[i])
    end
end

# Compute stages of fully implicit Runge-Kutta methods.
function function_stages!(x::Vector{ST}, b::Vector{ST}, params::ParametersFIRK{DT,TT,D},
                          caches::CacheDict) where {ST,DT,TT,D}
    # temporary variables
    local y1::ST
    local y2::ST

    # get cache for internal stages
    cache = caches[ST]

    # compute stages from nonlinear solver solution x
    compute_stages!(x, cache.Q, cache.V, cache.Y, params)

    # compute b = - (Y-AV)
    for i in eachindex(cache.Y)
        for k in eachindex(cache.Y[i])
            y1 = 0
            y2 = 0
            for j in eachindex(cache.V)
                y1 += params.tab.a[i,j] * cache.V[j][k]
                y2 += params.tab.â[i,j] * cache.V[j][k]
            end
            b[D*(i-1)+k] = - cache.Y[i][k] + params.Δt * (y1 + y2)
        end
    end
end

# Compute Jacobian of fully implicit Runge-Kutta methods.
function jacobian!(x::Vector{DT}, jac::Matrix{DT}, cache::IntegratorCacheFIRK{DT,D,S},
                   params::ParametersFIRK{DT,TT,D,S}) where {DT,TT,D,S}
    local tᵢ::TT

    for i in 1:S
        for k in 1:D
            cache.Y[i][k] = x[D*(i-1)+k]
            cache.Q[i][k] = params.q[k] + cache.Y[i][k]
        end
        tᵢ = params.t + params.Δt * params.tab.c[i]
        F = (v,q) -> params.equs[:v](tᵢ, q, v)
        ForwardDiff.jacobian!(cache.J[i], params.F, cache.ṽ, cache.Q[i], params.Jconfig)
    end

    jac .= 0

    @inbounds for j in 1:S
        for i in 1:S
            for l in 1:D
                for k in 1:D
                    m = D*(i-1)+k
                    n = D*(j-1)+l

                    if i == j && k == l
                        jac[m,n] += -1
                    end

                    jac[m,n] += params.Δt * params.tab.a[i,j] * cache.J[j][k,l]
                    jac[m,n] += params.Δt * params.tab.â[i,j] * cache.J[j][k,l]
                end
            end
        end
    end
end


function integrate_step!(int::IntegratorFIRK{DT,TT}, sol::AtomicSolutionODE{DT,TT},
                         cache::IntegratorCacheFIRK{DT}=int.caches[DT]) where {DT,TT}

    # update nonlinear solver parameters from atomic solution
    update_params!(int, sol)

    # compute initial guess
    initial_guess!(int, sol, cache)

    # reset atomic solution
    reset!(sol, timestep(int))

    # call nonlinear solver
    solve!(int.solver)

    # print solver status
    print_solver_status(int.solver.status, int.solver.params)

    # check if solution contains NaNs or error bounds are violated
    check_solver_status(int.solver.status, int.solver.params)

    # compute vector field at internal stages
    compute_stages!(int.solver.x, cache.Q, cache.V, cache.Y, int.params)

    # compute final update
    update_solution!(sol.q, sol.q̃, cache.V, tableau(int).b, tableau(int).b̂, timestep(int))

    # copy solution to initial guess
    update_vector_fields!(int.iguess, sol.t, sol.q, sol.v)
end
