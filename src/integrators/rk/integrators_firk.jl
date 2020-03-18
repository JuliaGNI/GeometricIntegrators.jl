
using ForwardDiff

"Holds the tableau of a fully implicit Runge-Kutta method."
struct TableauFIRK{T} <: AbstractTableauIRK{T}
    @HeaderTableau

    q::CoefficientsRK{T}

    function TableauFIRK{T}(q) where {T}
        if (q.s > 1 && istrilstrict(q.a)) || (q.s==1 && q.a[1,1] == 0)
            @warn "Initializing TableauFIRK with explicit tableau $(q.name).\nYou might want to use TableauERK instead."
        elseif q.s > 1 && istril(q.a)
            @warn "Initializing TableauFIRK with diagonally implicit tableau $(q.name).\nYou might want to use TableauDIRK instead."
        end

        new(q.name, q.o, q.s, q)
    end
end

function TableauFIRK(q::CoefficientsRK{T}) where {T}
    TableauFIRK{T}(q)
end

function TableauFIRK(name::Symbol, order::Int, a::Matrix{T}, b::Vector{T}, c::Vector{T}) where {T}
    TableauFIRK{T}(CoefficientsRK(name, order, a, b, c))
end

# TODO function readTableauFIRKFromFile(dir::AbstractString, name::AbstractString)


"Parameters for right-hand side function of fully implicit Runge-Kutta methods."
mutable struct ParametersFIRK{DT, TT, D, S, ET <: NamedTuple, FT, JT} <: Parameters{DT,TT}
    equs::ET
    tab::TableauFIRK{TT}
    Δt::TT

    F::FT
    Jconfig::JT

    t::TT
    q::Vector{DT}

    function ParametersFIRK{DT,D}(equs::ET, tab::TableauFIRK{TT}, Δt::TT) where {DT, TT, D, ET <: NamedTuple}
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


function IntegratorCache(params::ParametersFIRK{DT, TT, D, S}; kwargs...) where {DT, TT, D, S}
    IntegratorCacheFIRK{DT,D,S}(; kwargs...)
end


"Fully implicit Runge-Kutta integrator."
struct IntegratorFIRK{DT, TT, D, S, PT <: ParametersFIRK{DT,TT},
                                    ST <: NonlinearSolver{DT},
                                    IT <: InitialGuessODE{DT,TT}} <: IntegratorRK{DT,TT}
    params::PT
    solver::ST
    iguess::IT
    cache::IntegratorCacheFIRK{DT,D,S}

    function IntegratorFIRK(params::ParametersFIRK{DT,TT,D,S}, solver::ST, iguess::IT, cache) where {DT,TT,D,S,ST,IT}
        new{DT, TT, D, S, typeof(params), ST, IT}(params, solver, iguess, cache)
    end

    function IntegratorFIRK{DT,D}(equations::NamedTuple, tableau::TableauFIRK{TT}, Δt::TT; exact_jacobian=false) where {DT,TT,D}
        # get number of stages
        S = tableau.s

        # create params
        params = ParametersFIRK{DT,D}(equations, tableau, Δt)

        # create solver
        if exact_jacobian
            solver = create_nonlinear_solver_with_jacobian(DT, D*S, params)
        else
            solver = create_nonlinear_solver(DT, D*S, params)
        end

        # create initial guess
        iguess = InitialGuessODE{DT,D}(get_config(:ig_interpolation), equations[:v], Δt)

        # create cache
        cache = IntegratorCacheFIRK{DT,D,S}()

        # create integrator
        IntegratorFIRK(params, solver, iguess, cache)
    end

    function IntegratorFIRK{DT,D}(v::Function, tableau::TableauFIRK{TT}, Δt::TT; kwargs...) where {DT,TT,D}
        IntegratorFIRK{DT,D}(NamedTuple{(:v,)}((v,)), tableau, Δt; kwargs...)
    end

    function IntegratorFIRK{DT,D}(v::Function, h::Function, tableau::TableauFIRK{TT}, Δt::TT; kwargs...) where {DT,TT,D}
        IntegratorFIRK{DT,D}(NamedTuple{(:v,:h)}((v,h)), tableau, Δt; kwargs...)
    end

    function IntegratorFIRK(equation::ODE{DT,TT}, tableau::TableauFIRK{TT}, Δt::TT; kwargs...) where {DT,TT}
        IntegratorFIRK{DT, ndims(equation)}(get_function_tuple(equation), tableau, Δt; kwargs...)
    end
end


@inline Base.ndims(int::IntegratorFIRK{DT,TT,D,S}) where {DT,TT,D,S} = D


function initialize!(int::IntegratorFIRK, sol::AtomicSolutionODE)
    sol.t̅ = sol.t - timestep(int)

    equations(int)[:v](sol.t, sol.q, sol.v)

    initialize!(int.iguess, sol.t, sol.q, sol.v,
                            sol.t̅, sol.q̅, sol.v̅)
end


function update_params!(int::IntegratorFIRK, sol::AtomicSolutionODE)
    # set time for nonlinear solver and copy previous solution
    int.params.t  = sol.t
    int.params.q .= sol.q
end


function initial_guess!(int::IntegratorFIRK, sol::AtomicSolutionODE)
    local offset::Int

    # compute initial guess for internal stages
    for i in eachstage(int)
        evaluate!(int.iguess, sol.q, sol.v, sol.q̅, sol.v̅, int.cache.Q[i], int.cache.V[i], tableau(int).q.c[i])
    end
    for i in eachstage(int)
        offset = ndims(int)*(i-1)
        for k in eachdim(int)
            int.solver.x[offset+k] = 0
            for j in eachstage(int)
                int.solver.x[offset+k] += timestep(int) * tableau(int).q.a[i,j] * int.cache.V[j][k]
            end
        end
    end
end


function compute_stages!(x::Vector{ST}, Q::Vector{Vector{ST}}, V::Vector{Vector{ST}}, Y::Vector{Vector{ST}},
                         params::ParametersFIRK{DT,TT,D,S}) where {ST,DT,TT,D,S}

    local tᵢ::TT

    @assert S == length(Q) == length(V) == length(Y)

    # copy x to Y and compute Q = q + Δt Y
    for i in 1:S
        @assert D == length(Q[i]) == length(V[i]) == length(Y[i])
        for k in 1:D
            Y[i][k] = x[D*(i-1)+k]
            Q[i][k] = params.q[k] + Y[i][k]
        end
    end

    # compute V = v(Q)
    for i in 1:S
        tᵢ = params.t + params.Δt * params.tab.q.c[i]
        params.equs[:v](tᵢ, Q[i], V[i])
    end
end

"Compute stages of fully implicit Runge-Kutta methods."
function function_stages!(x::Vector{ST}, b::Vector{ST}, params::ParametersFIRK{DT,TT,D,S}) where {ST,DT,TT,D,S}
    # temporary variables
    local y1::ST
    local y2::ST

    # create caches for internal stages
    Q = create_internal_stage_vector(ST, D, S)
    V = create_internal_stage_vector(ST, D, S)
    Y = create_internal_stage_vector(ST, D, S)

    # compute stages from nonlinear solver solution x
    compute_stages!(x, Q, V, Y, params)

    # compute b = - (Y-AV)
    for i in 1:S
        for k in 1:D
            y1 = 0
            y2 = 0
            for j in 1:S
                y1 += params.tab.q.a[i,j] * V[j][k]
                y2 += params.tab.q.â[i,j] * V[j][k]
            end
            b[D*(i-1)+k] = - Y[i][k] + params.Δt * (y1 + y2)
        end
    end
end

"Compute stages of fully implicit Runge-Kutta methods."
function jacobian!(x::Vector{DT}, jac::Matrix{DT}, cache::IntegratorCacheFIRK{DT,D,S}, params::ParametersFIRK{DT,TT,D,S}) where {DT,TT,D,S}
    local tᵢ::TT

    for i in 1:S
        for k in 1:D
            cache.Y[i][k] = x[D*(i-1)+k]
            cache.Q[i][k] = params.q[k] + cache.Y[i][k]
        end
        tᵢ = params.t + params.Δt * params.tab.q.c[i]
        F = (v,q) -> params.equs[:v](tᵢ, q, v)
        ForwardDiff.jacobian!(cache.J[i], params.F, cache.ṽ, cache.Q[i], params.Jconfig)
    end

    jac .= 0

    for i in 1:S
        for j in 1:S
            for k in 1:D
                for l in 1:D
                    m = D*(i-1)+k
                    n = D*(j-1)+l

                    if i == j && k == l
                        jac[m,n] += -1
                    end

                    jac[m,n] += params.Δt * params.tab.q.a[i,j] * cache.J[j][k,l]
                    jac[m,n] += params.Δt * params.tab.q.â[i,j] * cache.J[j][k,l]
                end
            end
        end
    end
end


"Integrate ODE with fully implicit Runge-Kutta integrator."
function integrate_step!(int::IntegratorFIRK{DT,TT}, sol::AtomicSolutionODE{DT,TT}) where {DT,TT}
    # update nonlinear solver parameters from cache
    update_params!(int, sol)

    # compute initial guess
    initial_guess!(int, sol)

    # reset cache
    reset!(sol, timestep(int))

    # call nonlinear solver
    solve!(int.solver)

    # print solver status
    print_solver_status(int.solver.status, int.solver.params)

    # check if solution contains NaNs or error bounds are violated
    check_solver_status(int.solver.status, int.solver.params)

    # compute vector field at internal stages
    compute_stages!(int.solver.x, int.cache.Q, int.cache.V, int.cache.Y, int.params)

    # compute final update
    update_solution!(sol.q, sol.q̃, int.cache.V, int.params.tab.q.b, int.params.tab.q.b̂, int.params.Δt)

    # copy solution to initial guess
    update_vector_fields!(int.iguess, sol.t, sol.q, sol.v)
end
