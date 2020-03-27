"""
Holds the tableau of a weak implicit Runge-Kutta method.

    Reference: Wang, Hong, Xu, "Construction of Symplectic Runge-Kutta Methods for Stochastic Hamiltonian Systems",
    Commun. Comput. Phys. 21(1), 2017.

Order of the tableau is not included, because unlike in the deterministic
setting, it depends on the properties of the noise (e.g., the dimension of
the Wiener process and the commutativity properties of the diffusion matrix)

Orders stored in qdrift and qdiff are understood as the classical orders of these methods.

`qdrift0`, `qdrift1` correspond to A0, A1 in the paper
`qdiff0`, `qdiff1`, `qdiff3` correspond to B0, B1, B3
`qdrift0.b = alpha`
`qdiff0.b  = beta`
`qdrift0.c = c0`
`qdrift1.c = c1`
"""
struct TableauWIRK{T} <: AbstractTableauIRK{T}
    name::Symbol
    s::Int
    qdrift0::CoefficientsRK{T}
    qdrift1::CoefficientsRK{T}
    qdiff0::CoefficientsRK{T}
    qdiff1::CoefficientsRK{T}
    qdiff3::CoefficientsRK{T}

    function TableauWIRK{T}(name, s, qdrift0, qdrift1, qdiff0, qdiff1, qdiff3) where {T}
        @assert s == qdrift0.s == qdrift1.s == qdiff0.s == qdiff1.s == qdiff3.s
        new(name, s, qdrift0, qdrift1, qdiff0, qdiff1, qdiff3)
    end
end

function TableauWIRK(name::Symbol, qdrift0::CoefficientsRK{T}, qdrift1::CoefficientsRK{T}, qdiff0::CoefficientsRK{T}, qdiff1::CoefficientsRK{T}, qdiff3::CoefficientsRK{T}) where {T}
    TableauWIRK{T}(name, qdrift0.s, qdrift0, qdrift1, qdiff0, qdiff1, qdiff3)
end

function TableauWIRK(name::Symbol, A0::Matrix{T}, A1::Matrix{T},
                                   B0::Matrix{T}, B1::Matrix{T}, B3::Matrix{T},
                                    α::Vector{T}, β1::Vector{T},
                                   c0::Vector{T}, c1::Vector{T} ) where {T}
    TableauWIRK(name, CoefficientsRK(name, 0, A0, α, c0),
                      CoefficientsRK(name, 0, A1, α, c1),
                      CoefficientsRK(name, 0, B0, β1, c0),
                      CoefficientsRK(name, 0, B1, β1, c1),
                      CoefficientsRK(name, 0, B3, β1, c1))
end



"""
Parameters for right-hand side function of weak implicit Runge-Kutta methods.
"""
mutable struct ParametersWIRK{DT, TT, D, M, S, ET <: NamedTuple} <: Parameters{DT,TT}
    equ::ET
    tab::TableauWIRK{TT}
    Δt::TT
    ΔW::Vector{DT}
    ΔZ::Vector{DT}

    t::TT
    q::Vector{DT}

    function ParametersWIRK{DT,D,M}(equ::ET, tab::TableauWIRK{TT}, Δt::TT, ΔW::Vector{DT}, ΔZ::Vector{DT}) where {DT, TT, D, M, ET <: NamedTuple}
        @assert M == length(ΔW) == length(ΔZ)
        new{DT, TT, D, M, tab.s, ET}(equ, tab, Δt, ΔW, ΔZ, zero(TT), zeros(DT,D))
    end
end


"""
Structure for holding the internal stages `Q0`, and `Q1` the values of the drift vector
and the diffusion matrix evaluated at the internal stages `V=v(Q0)`, `B=B(Q1)`,
and the increments `Y = Δt*a_drift*v(Q) + a_diff*B(Q)*ΔW`.
"""
struct IntegratorCacheWIRK{DT,D,M,S} <: SDEIntegratorCache{DT,D,M}
    Q::OffsetArray{Array{DT,1},2,Array{Array{DT,1},2}}     # Q[l,i][k] - the k-th component of the internal stage Q^(l)_i, 0 ≤ l ≤ M
    Y::OffsetArray{Array{DT,1},2,Array{Array{DT,1},2}}     # Y[l,i][k] - the increment of the internal stage Q^(l)_i, 0 ≤ l ≤ M
    V::Vector{Vector{DT}}     # V [i][k]   - the k-th component of the drift vector v(Q0[i][:])
    B::Vector{Matrix{DT}}     # B [i][:,:] - the diffusion matrix at the stage i such that the l-th column B[i][:,l] is evaluated at Q^(l)_i

    v::Vector{DT}
    b::Matrix{DT}
    y::Vector{DT}

    Δy::Vector{DT}

    function IntegratorCacheWIRK{DT,D,M,S}() where {DT,D,M,S}
        # create internal stage vectors
        Q  = create_internal_stage_vector_with_zero(DT, D, M, S)
        Y  = create_internal_stage_vector_with_zero(DT, D, M, S)
        V  = create_internal_stage_vector(DT, D, S)
        B  = create_internal_stage_matrix(DT, D, M, S)

        # create velocity and update vector
        v = zeros(DT,D)
        b = zeros(DT,D,M)
        y = zeros(DT,D)

        # create temporary vectors
        Δy  = zeros(DT,M)

        new(Q, Y, V, B, v, b, y, Δy)
    end
end

function Integrators.IntegratorCache{ST}(params::ParametersWIRK{DT,TT,D,M,S}; kwargs...) where {ST,DT,TT,D,M,S}
    IntegratorCacheWIRK{ST,D,M,S}(; kwargs...)
end

@inline Integrators.CacheType(ST, params::ParametersWIRK{DT,TT,D,M,S}) where {DT,TT,D,M,S} = IntegratorCacheWIRK{ST,D,M,S}


"""
Stochastic implicit Runge-Kutta integrator.
"""
struct IntegratorWIRK{DT, TT, D, M, S,
                PT <: ParametersWIRK{DT,TT},
                ST <: NonlinearSolver{DT}} <: StochasticIntegratorRK{DT,TT,D,M,S}
    params::PT
    solver::ST
    caches::CacheDict{PT}

    function IntegratorWIRK(params::ParametersWIRK{DT,TT,D,M,S}, solver::ST, caches) where {DT,TT,D,M,S,ST}
        new{DT, TT, D, M, S, typeof(params), ST}(params, solver, caches)
    end

    function IntegratorWIRK{DT,D,M}(equations::NamedTuple, tableau::TableauWIRK{TT}, Δt::TT) where {DT,TT,D,M}
        # get number of stages
        S = tableau.s

        # create params
        params = ParametersWIRK{DT,D,M}(equations, tableau, Δt, zeros(DT,M), zeros(DT,M))

        # create cache dict
        caches = CacheDict(params)

        # create solver
        solver = create_nonlinear_solver(DT, D*S*(M+1), params, caches)

        # create integrator
        IntegratorWIRK(params, solver, caches)
    end

    function IntegratorWIRK(equation::SDE{DT,TT}, tableau::TableauWIRK{TT}, Δt::TT; kwargs...) where {DT,TT}
        IntegratorWIRK{DT, ndims(equation), equation.m}(get_function_tuple(equation), tableau, Δt; kwargs...)
    end
end


"""
This function computes initial guesses for `Y` and assigns them to int.solver.x
For WIRK we are NOT IMPLEMENTING an InitialGuess.

Using an explicit integrator to predict the next step's value (like in SIRK)
does not seem to be a good idea here, because the integrators are convergent
in the weak sense only, and there is no guarantee that the explicit integrator
will produce anything close to the desired solution...

The simplest initial guess for `Y` is 0.
"""
function initial_guess!(int::IntegratorWIRK{DT,TT}, sol::AtomicSolutionSDE{DT,TT},
                        cache::IntegratorCacheWIRK{DT}=int.caches[DT]) where {DT,TT}
    int.solver.x .= 0
end


"""
Unpacks the data stored in
`x = (Y0[1][1], Y0[1][2]], ... Y0[1][D], ... Y0[S][D], Y1[1][1,1], Y1[1][2,1], ... Y1[1][D,1], Y1[1][1,2], Y1[1][2,2], ... Y1[1][D,2], ... Y1[S][D,M] )`
into `Y0` and `Y1`, calculates the internal stages `Q0` and `Q1`, the values of the RHS
of the SDE ( `v(Q0)` and `B(Q1)` ), and assigns them to `V` and `B`.
Unlike for FIRK, here `Y = Δt a v(Q) + â B(Q) ΔW`.
"""
function compute_stages!(x::Vector{ST}, Q::AbstractMatrix{Vector{ST}}, V::Vector{Vector{ST}},
                                        B::Vector{Matrix{ST}}, Y::AbstractMatrix{Vector{ST}},
                                        params::ParametersWIRK{DT,TT,D,M,S}) where {ST,DT,TT,D,M,S}

    @assert S == size(Q,2) == size(Y,2) == length(V) == length(B)

    local tᵢ::TT

    # copy x to Y and calculate Q
    for i in axes(Q,2)
        for k in eachindex(Y[0,i])
            Y[0,i][k] = x[D*(i-1)+k]
        end
        for l in 1:M
            for k in eachindex(Y[l,i])
                Y[l,i][k] = x[ D*S + D*M*(i-1) + D*(l-1) + k ]
            end
        end
        for l in axes(Q,1)
            for k in eachindex(Q[l,i])
                Q[l,i][k] = params.q[k] + Y[l,i][k]
            end
        end
    end

    # compute V = v(Q0) and B=B(Q1)
    for i in 1:S
        # time point with the coefficient c0
        tᵢ = params.t + params.Δt * params.tab.qdrift0.c[i]
        # calculates v(t,tQ) and assigns to the i-th column of V
        params.equ[:v](tᵢ, Q[0,i], V[i])

        # time point with the coefficient c1
        tᵢ = params.t + params.Δt * params.tab.qdrift1.c[i]

        for l=1:M
            # calculates the l-th column of B(t,Q) and assigns to the l-th column of the Matrix B
            params.equ[:B](tᵢ, Q[l,i], B[i], l)
        end
    end
end

"""
Compute stages of weak implicit Runge-Kutta methods.
"""
function Integrators.function_stages!(x::Vector{ST}, b::Vector{ST},
                params::ParametersWIRK{DT,TT,D,M,S},
                caches::CacheDict) where {ST,DT,TT,D,M,S}

    # get cache for internal stages
    cache = caches[ST]

    compute_stages!(x, cache.Q, cache.V, cache.B, cache.Y, params)

    local y1::ST
    local y2::ST

    # compute b = - (Y-AV)
    # the terms corresponding to Y0
    for i in 1:S
        for k in 1:D
            y1 = 0
            y2 = 0
            for j in 1:S
                y1 += params.tab.qdrift0.a[i,j] * cache.V[j][k] * params.Δt
                y2 += params.tab.qdrift0.â[i,j] * cache.V[j][k] * params.Δt
                for l in 1:M
                    y1 += params.tab.qdiff0.a[i,j] * cache.B[j][k,l] * params.ΔW[l]
                    y2 += params.tab.qdiff0.â[i,j] * cache.B[j][k,l] * params.ΔW[l]
                end
            end
            b[D*(i-1)+k] = - cache.Y[0,i][k] + (y1 + y2)
        end
    end

    # compute b = - (Y-AV)
    # the terms corresponding to Y1
    for i in 1:S
        for l in 1:M
            for k in 1:D
                y1 = 0
                y2 = 0
                for j in 1:S
                    # The drift terms
                    y1 += params.tab.qdrift1.a[i,j] * cache.V[j][k] * params.Δt
                    y2 += params.tab.qdrift1.â[i,j] * cache.V[j][k] * params.Δt

                    # The noise terms, calculated either with the B1 or B3 tableau
                    for noise_idx in 1:M
                        if noise_idx==l
                            y1 += params.tab.qdiff1.a[i,j] * cache.B[j][k,noise_idx] * params.ΔW[noise_idx]
                            y2 += params.tab.qdiff1.â[i,j] * cache.B[j][k,noise_idx] * params.ΔW[noise_idx]
                        else
                            y1 += params.tab.qdiff3.a[i,j] * cache.B[j][k,noise_idx] * params.ΔW[noise_idx]
                            y2 += params.tab.qdiff3.â[i,j] * cache.B[j][k,noise_idx] * params.ΔW[noise_idx]
                        end
                    end
                end
                b[D*S + D*M*(i-1) + D*(l-1) + k] = - cache.Y[l,i][k] + (y1 + y2)
            end
        end
    end
end


"""
Integrate SDE with a stochastic implicit Runge-Kutta integrator.
"""
function Integrators.integrate_step!(int::IntegratorWIRK{DT,TT}, sol::AtomicSolutionSDE{DT,TT},
                                     cache::IntegratorCacheWIRK{DT}=int.caches[DT]) where {DT,TT}
    # update nonlinear solver parameters from cache
    update_params!(int, sol)

    # compute initial guess and assign to int.solver.x
    initial_guess!(int, sol, cache)

    # reset cache
    reset!(sol, timestep(int))

    # call nonlinear solver
    solve!(int.solver)

    # print solver status
    print_solver_status(int.solver.status, int.solver.params)

    # check if solution contains NaNs or error bounds are violated
    check_solver_status(int.solver.status, int.solver.params)

    # compute the drift vector field and the diffusion matrix at internal stages
    compute_stages!(int.solver.x, cache.Q, cache.V, cache.B, cache.Y, int.params)

    # compute final update (same update function as for SIRK)
    update_solution!(sol, cache.V, cache.B, tableau(int).qdrift0.b, tableau(int).qdiff0.b, timestep(int), int.params.ΔW, cache.Δy)
    update_solution!(sol, cache.V, cache.B, tableau(int).qdrift0.b̂, tableau(int).qdiff0.b̂, timestep(int), int.params.ΔW, cache.Δy)
end
