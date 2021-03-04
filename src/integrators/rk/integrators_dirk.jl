
"Parameters for right-hand side function of diagonally implicit Runge-Kutta methods."
mutable struct ParametersDIRK{DT, TT, D, S, ET <: NamedTuple} <: Parameters{DT,TT}
    equs::ET
    tab::Tableau{TT}
    Δt::TT

    t::TT
    q::Vector{DT}
    V::Vector{Vector{DT}}

    function ParametersDIRK{DT,D}(equs::ET, tab::Tableau{TT}, Δt::TT) where {DT, TT, D, ET <: NamedTuple}
        new{DT, TT, D, tab.s, ET}(equs, tab, Δt, zero(TT), zeros(DT,D), create_internal_stage_vector(DT, D, tab.s))
    end
end


"""
Diagonally implicit Runge-Kutta integrator cache.
"""
struct IntegratorCacheDIRK{DT,D,S} <: ODEIntegratorCache{DT,D}
    q̃::Vector{DT}
    ṽ::Vector{DT}
    s̃::Vector{DT}
    y::Vector{DT}

    Q::Vector{Vector{DT}}
    V::Vector{Vector{DT}}
    Y::Vector{Vector{DT}}

    function IntegratorCacheDIRK{DT,D,S}() where {DT,D,S}
        Q = create_internal_stage_vector(DT, D, S)
        V = create_internal_stage_vector(DT, D, S)
        Y = create_internal_stage_vector(DT, D, S)
        new(zeros(DT,D), zeros(DT,D), zeros(DT,D), zeros(DT,D), Q, V, Y)
    end
end

function IntegratorCache{ST}(params::ParametersDIRK{DT,TT,D,S}; kwargs...) where {ST,DT,TT,D,S}
    IntegratorCacheDIRK{ST,D,S}(; kwargs...)
end

@inline CacheType(ST, params::ParametersDIRK{DT,TT,D,S}) where {DT,TT,D,S} = IntegratorCacheDIRK{ST,D,S}


"Diagonally implicit Runge-Kutta integrator."
struct IntegratorDIRK{DT, TT, D, S, PT <: ParametersDIRK{DT,TT},
                                    ST,# <: NonlinearSolver{DT},
                                    IT <: InitialGuessODE{TT}} <: AbstractIntegratorRK{DT,TT}
    params::PT
    solver::ST
    iguess::IT
    caches::CacheDict{PT}

    function IntegratorDIRK(params::ParametersDIRK{DT,TT,D,S}, solver::ST, iguess::IT, caches) where {DT,TT,D,S,ST,IT}
        new{DT, TT, D, S, typeof(params), ST, IT}(params, solver, iguess, caches)
    end

    function IntegratorDIRK{DT,D}(equations::NamedTuple, tableau::Tableau{TT}, Δt::TT) where {DT, TT, D}
        # get number of stages
        S = tableau.s

        # check if tableau is diagonally implicit
        @assert !isfullyimplicit(tableau)

        if get_config(:verbosity) ≥ 1
            if isexplicit(tableau)
                @warn "Initializing IntegratorDIRK with explicit tableau $(tableau.name).\nYou might want to use IntegratorERK instead."
            end
        end

        # create params
        params = ParametersDIRK{DT,D}(equations, tableau, Δt)

        # create cache dict
        caches = CacheDict(params)

        # create solvers
        solvers = [create_nonlinear_solver(DT, D, params, caches, i) for i in 1:S]

        # create initial guess
        iguess = InitialGuessODE(get_config(:ig_interpolation), equations[:v], Δt)

        # create integrator
        IntegratorDIRK(params, solvers, iguess, caches)
    end

    function IntegratorDIRK{DT,D}(v::Function, tableau::Tableau{TT}, Δt::TT; kwargs...) where {DT,TT,D}
        IntegratorDIRK{DT,D}(NamedTuple{(:v,)}((v,)), tableau, Δt; kwargs...)
    end

    function IntegratorDIRK{DT,D}(v::Function, h::Function, tableau::Tableau{TT}, Δt::TT; kwargs...) where {DT,TT,D}
        IntegratorDIRK{DT,D}(NamedTuple{(:v,:h)}((v,h)), tableau, Δt; kwargs...)
    end

    function IntegratorDIRK(equation::ODE{DT}, tableau::Tableau{TT}, Δt::TT; kwargs...) where {DT,TT}
        IntegratorDIRK{DT, ndims(equation)}(get_function_tuple(equation), tableau, Δt; kwargs...)
    end
end


@inline Base.ndims(::IntegratorDIRK{DT,TT,D,S}) where {DT,TT,D,S} = D
@inline has_initial_guess(int::IntegratorDIRK) = true


"Initialise initial guess"
function initialize!(int::IntegratorDIRK, cache::IntegratorCacheDIRK)
    # initialise initial guess
    cache.t̄ = cache.t - timestep(int)

    int.params.equs[:v](cache.t, cache.q, cache.v)

    initialize!(int.iguess, cache.t, cache.q, cache.v, cache.t̄, cache.q̄, cache.v̄)
end


function update_params!(int::IntegratorDIRK, sol::AtomicSolutionODE)
    # set time for nonlinear solver and copy previous solution
    int.params.t  = sol.t
    int.params.q .= sol.q
end


"Compute initial guess for internal stages."
function initial_guess!(int::IntegratorDIRK{DT}, sol::AtomicSolutionODE{DT},
                        cache::IntegratorCacheDIRK{DT}=int.caches[DT]) where {DT}

    for i in eachstage(int)
        evaluate!(int.iguess, sol.q̄, sol.v̄, sol.q, sol.v, cache.q̃, cache.ṽ, int.params.tab.c[i])
        for k in eachindex(cache.V[i], cache.ṽ)
            cache.V[i][k] = cache.ṽ[k]
        end
    end
    for i in eachstage(int)
        for k in eachdim(int)
            int.solver[i].x[k] = 0
            for j in eachstage(int)
                int.solver[i].x[k] += int.params.tab.a[i,j] * cache.V[j][k]
            end
        end
    end
end


function compute_stages!(x::Vector{ST}, Q::Vector{ST}, V::Vector{ST}, Y::Vector{ST},
                         params::ParametersDIRK{DT,TT,D,S}, i::Int) where {ST,DT,TT,D,S}

    local tᵢ::TT

    @assert D == length(Q) == length(V) == length(Y)

    # copy x to Y and compute Q = q + Δt Y
    for k in 1:D
        Y[k] = x[k]
        Q[k] = params.q[k] + params.Δt * Y[k]
    end

    # compute V = v(Q)
    tᵢ = params.t + params.Δt * params.tab.c[i]
    params.equs[:v](tᵢ, Q, V)
end


"Compute stages of fully implicit Runge-Kutta methods."
function function_stages!(x::Vector{ST}, b::Vector{ST}, params::ParametersDIRK{DT,TT,D,S},
                          caches::CacheDict, i::Int) where {ST,DT,TT,D,S}

    # temporary variables
    local y1::ST
    local y2::ST

    # get cache for internal stages
    cache = caches[ST]

    # compute stages from nonlinear solver solution x
    compute_stages!(x, cache.Q[i], cache.V[i], cache.Y[i], params, i)

    # compute b = - (Y-AV)
    for k in 1:D
        y1 = params.tab.a[i,i] * cache.V[i][k]
        y2 = params.tab.â[i,i] * cache.V[i][k]
        for j in 1:i-1
            y1 += params.tab.a[i,j] * params.V[j][k]
            y2 += params.tab.â[i,j] * params.V[j][k]
        end
        b[k] = - cache.Y[i][k] + (y1 + y2)
    end
end


function integrate_step!(int::IntegratorDIRK{DT,TT}, sol::AtomicSolutionODE{DT,TT},
                         cache::IntegratorCacheDIRK{DT}=int.caches[DT]) where {DT,TT}

     # update nonlinear solver parameters from atomic solution
     update_params!(int, sol)

    # compute initial guess
    initial_guess!(int, sol, cache)

    # reset atomic solution
    reset!(sol, timestep(int))

    # consecutively solve for all stages
    for i in eachstage(int)
        # call nonlinear solver
        solve!(int.solver[i])

        # print solver status
        print_solver_status(int.solver[i].status, int.solver[i].params)

        # check if solution contains NaNs or error bounds are violated
        check_solver_status(int.solver[i].status, int.solver[i].params)

        # compute vector field at internal stages
        compute_stages!(int.solver[i].x, cache.Q[i], cache.V[i], cache.Y[i], int.params, i)

        # copy velocity field to parameters
        int.params.V[i] .= cache.V[i]
    end

    # compute final update
    update_solution!(sol.q, sol.q̃, cache.V, tableau(int).b, tableau(int).b̂, timestep(int))

    # update vector field for initial guess
    update_vector_fields!(int.iguess, sol.t, sol.q, sol.v)
end
