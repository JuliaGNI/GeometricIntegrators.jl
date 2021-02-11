
using ForwardDiff

"Parameters for right-hand side function of formal Lagrangian Runge-Kutta methods."
mutable struct ParametersFLRK{DT, TT, D, S, ET <: NamedTuple} <: Parameters{DT,TT}
    equs::ET
    tab::Tableau{TT}
    Δt::TT

    t::TT
    q::Vector{DT}
    p::Vector{DT}

    function ParametersFLRK{DT,D}(equs::ET, tab::Tableau{TT}, Δt::TT) where {DT, TT, D, ET <: NamedTuple}
        new{DT,TT,D,tab.s,ET}(equs, tab, Δt, zero(TT), zeros(DT,D), zeros(DT,D))
    end
end

struct IntegratorCacheFLRK{DT,D,S} <: ODEIntegratorCache{DT,D}
    q̃::Vector{DT}
    ṽ::Vector{DT}
    s̃::Vector{DT}

    Q::Vector{Vector{DT}}
    V::Vector{Vector{DT}}
    Y::Vector{Vector{DT}}

    function IntegratorCacheFLRK{DT,D,S}() where {DT,D,S}
        # create internal stage vectors
        Q = create_internal_stage_vector(DT, D, S)
        V = create_internal_stage_vector(DT, D, S)
        Y = create_internal_stage_vector(DT, D, S)

        new(zeros(DT,D), zeros(DT,D), zeros(DT,D), Q, V, Y)
    end
end

function IntegratorCache{ST}(params::ParametersFLRK{DT,TT,D,S}; kwargs...) where {ST,DT,TT,D,S}
    IntegratorCacheFLRK{ST,D,S}(; kwargs...)
end

@inline CacheType(ST, params::ParametersFLRK{DT,TT,D,S}) where {DT,TT,D,S} = IntegratorCacheFLRK{ST,D,S}


"Formal Lagrangian Runge-Kutta integrator."
struct IntegratorFLRK{DT, TT, D, S, PT <: ParametersFLRK{DT,TT},
                                    ST <: NonlinearSolver{DT},
                                    IT <: InitialGuessODE{TT}} <: AbstractIntegratorPRK{DT,TT}
    params::PT
    solver::ST
    iguess::IT
    caches::CacheDict{PT}

    ϑ::Vector{Vector{DT}}
    P::Vector{Vector{DT}}
    F::Vector{Vector{DT}}
    G::Vector{Vector{DT}}
    Z::Vector{Vector{DT}}
    J::Vector{Matrix{DT}}
    A::Array{DT,2}

    function IntegratorFLRK(params::ParametersFLRK{DT,TT,D,S}, solver::ST, iguess::IT, caches) where {DT,TT,D,S,ST,IT}
        # create internal stage vectors
        ϑ = create_internal_stage_vector(DT, D, S)
        P = create_internal_stage_vector(DT, D, S)
        F = create_internal_stage_vector(DT, D, S)
        G = create_internal_stage_vector(DT, D, S)
        Z = create_internal_stage_vector(DT, D, S)

        J = [zeros(DT,D,D) for i in 1:S]
        A = zeros(DT,D*S,D*S)

        new{DT, TT, D, S, typeof(params), ST, IT}(params, solver, iguess, caches, ϑ, P, F, G, Z, J, A)
    end

    function IntegratorFLRK{DT,D}(equations::NamedTuple, tableau::Tableau{TT}, Δt::TT) where {DT,TT,D}
        # get number of stages
        S = tableau.s

        # check if tableau is fully implicit
        if get_config(:verbosity) ≥ 1
            if isexplicit(tableau)
                @warn "Initializing IntegratorFIRK with explicit tableau $(q.name).\nYou might want to use IntegratorERK instead."
            elseif isdiagnonallyimplicit(tableau)
                @warn "Initializing IntegratorFIRK with diagonally implicit tableau $(q.name).\nYou might want to use IntegratorDIRK instead."
            end
        end

        # create params
        params = ParametersFLRK{DT,D}(equations, tableau, Δt)

        # create cache dict
        caches = CacheDict(params)

        # create solver
        solver = create_nonlinear_solver(DT, D*S, params, caches)

        # create initial guess
        iguess = InitialGuessODE(get_config(:ig_interpolation), equations[:v̄], Δt)

        # create integrator
        IntegratorFLRK(params, solver, iguess, caches)
    end

    function IntegratorFLRK(equation::LODE{DT,TT}, tableau::Tableau{TT}, Δt::TT; kwargs...) where {DT,TT}
        IntegratorFLRK{DT, equation.d}(get_function_tuple(equation), tableau, Δt; kwargs...)
    end

end


@inline Base.ndims(::IntegratorFLRK{DT,TT,D,S}) where {DT,TT,D,S} = D


function initialize!(int::IntegratorFLRK, sol::AtomicSolutionODE)
    sol.t̄ = sol.t - timestep(int)

    equations(int)[:v̄](sol.t, sol.q, sol.v)

    initialize!(int.iguess, sol.t, sol.q, sol.v,
                            sol.t̄, sol.q̄, sol.v̄)
end


function update_params!(int::IntegratorFLRK, sol::AtomicSolutionODE)
    # set time for nonlinear solver and copy previous solution
    int.params.t  = sol.t
    int.params.q .= sol.q
end

function update_params!(int::IntegratorFLRK, sol::AtomicSolutionPODE)
    # set time for nonlinear solver and copy previous solution
    int.params.t  = sol.t
    int.params.q .= sol.q
    int.params.p .= sol.p
end


function initial_guess!(int::IntegratorFLRK{DT}, sol::Union{AtomicSolutionODE{DT},AtomicSolutionPODE{DT}},
                        cache::IntegratorCacheFLRK{DT}=int.caches[DT]) where {DT}
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
                         params::ParametersFLRK{DT,TT,D,S}) where {ST,DT,TT,D,S}

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
        tᵢ = params.t + params.Δt * params.tab.c[i]
        params.equs[:v̄](tᵢ, Q[i], V[i])
    end
end

"Compute stages of formal Lagrangian Runge-Kutta methods."
function function_stages!(x::Vector{ST}, b::Vector{ST}, params::ParametersFLRK{DT,TT,D,S},
                          caches::CacheDict) where {ST,DT,TT,D,S}
    # temporary variables
    local y1::ST
    local y2::ST

    # get cache for internal stages
    cache = caches[ST]

    # compute stages from nonlinear solver solution x
    compute_stages!(x, cache.Q, cache.V, cache.Y, params)

    # compute b = - (Y-AV)
    for i in 1:S
        for k in 1:D
            y1 = 0
            y2 = 0
            for j in 1:S
                y1 += params.tab.a[i,j] * cache.V[j][k]
                y2 += params.tab.â[i,j] * cache.V[j][k]
            end
            b[D*(i-1)+k] = - cache.Y[i][k] + params.Δt * (y1 + y2)
        end
    end
end


function integrate_step!(int::IntegratorFLRK, sol::AtomicSolutionODE)
    integrate_step_flrk!(int, sol)
end

function integrate_step!(int::IntegratorFLRK, sol::AtomicSolutionPODE)
    integrate_step_flrk!(int, sol)
    integrate_diag_flrk!(int, sol)
end

function integrate_step_flrk!(int::IntegratorFLRK{DT,TT}, sol::AtomicSolutionPODE{DT,TT},
                              cache::IntegratorCacheFLRK{DT}=int.caches[DT]) where {DT,TT}
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

function integrate_diag_flrk!(int::IntegratorFLRK{DT,TT,D,S}, sol::AtomicSolutionPODE{DT,TT},
                              cache::IntegratorCacheFLRK{DT}=int.caches[DT]) where {DT,TT,D,S}
    # create temporary arrays
    δP = zeros(DT, D*S)

    # compute ϑ = ϑ(t,Q), V(Q) = v(t, Q, V)
    # and f_0(Q, V(Q)) = f(t, Q, V, F)
    for i in 1:S
        tᵢ = int.params.t + timestep(int) * tableau(int).c[i]
        int.params.equs[:ϑ](tᵢ, cache.Q[i], cache.V[i], int.ϑ[i])
        int.params.equs[:v̄](tᵢ, cache.Q[i], cache.V[i])
        int.params.equs[:g](tᵢ, cache.Q[i], cache.V[i], int.F[i])
    end

    # compute Jacobian of v via ForwardDiff
    for i in 1:S
        tᵢ = int.params.t + timestep(int) * tableau(int).c[i]
        v_rev! = (v,q) -> int.params.equs[:v̄](tᵢ,q,v)
        ForwardDiff.jacobian!(int.J[i], v_rev!, cache.ṽ, cache.Q[i])
    end

    # contract J with ϑ and add to G
    for l in 1:S
        for i in 1:D
            int.G[l][i] = 0
            for j in 1:D
                int.G[l][i] += int.ϑ[l][j] * int.J[l][j,i]
            end
        end
    end

    # solve linear system AP=δP for P

    # compute δP
    for l in 1:S
        for i in 1:D
            # set δP = p
            δP[(l-1)*D+i] = int.params.p[i]
            # add A(F+G) to δP
            for k in 1:S
                δP[(l-1)*D+i] += timestep(int) * tableau(int).a[l,k] * (int.F[k][i] + int.G[k][i])
            end
        end
    end

    # construct A = identity(sd×sd) + A ⊗ J
    for k in 1:S
        for l in 1:S
            for i in 1:D
                for j in 1:D
                    int.A[(k-1)*D+i, (l-1)*D+j] = timestep(int) * tableau(int).a[k,l] * int.J[l][j,i]
                end
            end
        end
    end
    for i in 1:D*S
        int.A[i,i] += one(DT)
    end

    # solve AP = δP and copy result to int.P
    lu = LUSolver(int.A, δP)
    factorize!(lu)
    solve!(lu)
    tP = reshape(lu.x, (D,S))
    for l in 1:S
        int.P[l] .= tP[:,l]
    end

    # contract J with P and subtract from G, so that G = (ϑ-P)J
    for l in 1:S
        for i in 1:D
            int.G[l][i] = 0
            for j in 1:D
                int.G[l][i] += (int.ϑ[l][j] - int.P[l][j]) * int.J[l][j,i]
            end
        end
    end

    # println(int.ϑ)
    # println(int.P)
    # println(int.ϑ .- int.P)
    # println()

    # compute final update for p
    update_solution!(sol.p, sol.p̃, int.F, tableau(int).b, tableau(int).b̂, timestep(int))
    update_solution!(sol.p, sol.p̃, int.G, tableau(int).b, tableau(int).b̂, timestep(int))
end
