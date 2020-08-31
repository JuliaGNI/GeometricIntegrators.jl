@doc raw"""
Holds the tableau of a stochastic implicit partitioned Runge-Kutta method.
qdrift, pdrift hold the RK coefficients for the drift part,
and qdiff, pdiff hold the RK coefficients for the diffusion part of the SDE.

Order of the tableau is not included, because unlike in the deterministic
setting, it depends on the properties of the noise (e.g., the dimension of
the Wiener process and the commutativity properties of the diffusion matrix).

Orders stored in qdrift and qdiff are understood as the classical orders of these methods.
"""
struct TableauSIPRK{T} <: AbstractTableauIRK{T}
    name::Symbol
    s::Int
    qdrift::CoefficientsRK{T}
    qdiff::CoefficientsRK{T}
    pdrift::CoefficientsRK{T}
    pdiff::CoefficientsRK{T}

    function TableauSIPRK{T}(name, s, qdrift, qdiff, pdrift, pdiff) where {T}
        @assert s == qdrift.s == qdiff.s == pdrift.s == pdiff.s
        new(name, s, qdrift, qdiff, pdrift, pdiff)
    end
end

function TableauSIPRK(name::Symbol, qdrift::CoefficientsRK{T}, qdiff::CoefficientsRK{T}, pdrift::CoefficientsRK{T}, pdiff::CoefficientsRK{T}) where {T}
    TableauSIPRK{T}(name, qdrift.s, qdrift, qdiff, pdrift, pdiff)
end

function TableauSIPRK(name::Symbol, qorder_drift::Int, qa_drift::Matrix{T}, qb_drift::Vector{T}, qc_drift::Vector{T}, qorder_diff::Int, qa_diff::Matrix{T}, qb_diff::Vector{T}, qc_diff::Vector{T},
                                    porder_drift::Int, pa_drift::Matrix{T}, pb_drift::Vector{T}, pc_drift::Vector{T}, porder_diff::Int, pa_diff::Matrix{T}, pb_diff::Vector{T}, pc_diff::Vector{T}) where {T}
    TableauSIPRK{T}(name, length(qc_drift), CoefficientsRK(name, qorder_drift, qa_drift, qb_drift, qc_drift), CoefficientsRK(name, qorder_diff, qa_diff, qb_diff, qc_diff), CoefficientsRK(name, porder_drift, pa_drift, pb_drift, pc_drift), CoefficientsRK(name, porder_diff, pa_diff, pb_diff, pc_diff))
end

# TODO function readTableauSFIRKFromFile(dir::AbstractString, name::AbstractString)


@doc raw"""
Parameters for right-hand side function of implicit Runge-Kutta methods.

`A` - if positive, the upper bound of the Wiener process increments; if A=0.0, no truncation
"""
mutable struct ParametersSIPRK{DT, TT, D, M, S, ET <: NamedTuple} <: Parameters{DT,TT}
    equ::ET
    tab::TableauSIPRK{TT}
    Δt::TT
    ΔW::Vector{DT}
    ΔZ::Vector{DT}
    A::DT

    t::TT
    q::Vector{DT}
    p::Vector{DT}

    function ParametersSIPRK{DT,D,M}(equ::ET, tab::TableauSIPRK{TT}, Δt::TT, ΔW::Vector{DT}, ΔZ::Vector{DT}, A::DT) where {DT, TT, D, M, ET <: NamedTuple}
        @assert M == length(ΔW) == length(ΔZ)
        new{DT, TT, D, M, tab.s, ET}(equ, tab, Δt, ΔW, ΔZ, A, zero(TT), zeros(DT,D), zeros(DT,D))
    end
end


@doc raw"""
Structure for holding the internal stages Q, the values of the drift vector
and the diffusion matrix evaluated at the internal stages `V=v(Q)`, `B=B(Q)`,
and the increments `Y = Δt*a_drift*v(Q) + a_diff*B(Q)*ΔW`.
"""
struct IntegratorCacheSIPRK{DT,D,M,S} <: PSDEIntegratorCache{DT,D,M}
    Q::Vector{Vector{DT}}
    P::Vector{Vector{DT}}
    V::Vector{Vector{DT}}
    F::Vector{Vector{DT}}
    B::Vector{Matrix{DT}}
    G::Vector{Matrix{DT}}
    Y::Vector{Vector{DT}}
    Z::Vector{Vector{DT}}

    v::Vector{DT}
    f::Vector{DT}
    b::Matrix{DT}
    g::Matrix{DT}
    y::Vector{DT}
    z::Vector{DT}

    V1::Vector{DT}
    V2::Vector{DT}
    F1::Vector{DT}
    F2::Vector{DT}
    B1::Matrix{DT}
    B2::Matrix{DT}
    G1::Matrix{DT}
    G2::Matrix{DT}
    Δw::Vector{DT}
    ΔQ::Vector{DT}
    Y1::Vector{DT}
    Y2::Vector{DT}
    ΔP::Vector{DT}
    Z1::Vector{DT}
    Z2::Vector{DT}
    Δy::Vector{DT}
    Δz::Vector{DT}

    function IntegratorCacheSIPRK{DT,D,M,S}() where {DT,D,M,S}

        # create internal stage vectors
        Q = create_internal_stage_vector(DT, D, S)
        P = create_internal_stage_vector(DT, D, S)
        V = create_internal_stage_vector(DT, D, S)
        F = create_internal_stage_vector(DT, D, S)
        B = create_internal_stage_matrix(DT, D, M, S)
        G = create_internal_stage_matrix(DT, D, M, S)
        Y = create_internal_stage_vector(DT, D, S)
        Z = create_internal_stage_vector(DT, D, S)

        # create velocity and update vector
        v = zeros(DT,D)
        f = zeros(DT,D)
        b = zeros(DT,D,M)
        g = zeros(DT,D,M)
        y = zeros(DT,D)
        z = zeros(DT,D)

        # create temporary arrays
        V1 = zeros(DT,D)
        V2 = zeros(DT,D)
        F1 = zeros(DT,D)
        F2 = zeros(DT,D)
        B1 = zeros(DT,D,M)
        B2 = zeros(DT,D,M)
        G1 = zeros(DT,D,M)
        G2 = zeros(DT,D,M)
        ΔQ = zeros(DT,D)
        Y1 = zeros(DT,D)
        Y2 = zeros(DT,D)
        ΔP = zeros(DT,D)
        Z1 = zeros(DT,D)
        Z2 = zeros(DT,D)
        Δw = zeros(DT,M)
        Δy = zeros(DT,M)
        Δz = zeros(DT,M)

        new(Q, P, V, F, B, G, Y, Z, v, f, b, g, y, z, V1, V2, F1, F2, B1, B2, G1, G2, Δw, ΔQ, Y1, Y2, ΔP, Z1, Z2, Δy, Δz)
    end
end

function Integrators.IntegratorCache{ST}(params::ParametersSIPRK{DT,TT,D,M,S}; kwargs...) where {ST,DT,TT,D,M,S}
    IntegratorCacheSIPRK{ST,D,M,S}(; kwargs...)
end

@inline Integrators.CacheType(ST, params::ParametersSIPRK{DT,TT,D,M,S}) where {DT,TT,D,M,S} = IntegratorCacheSIPRK{ST,D,M,S}


"Stochastic implicit partitioned Runge-Kutta integrator."
struct IntegratorSIPRK{DT, TT, D, M, S,
                PT <: ParametersSIPRK{DT,TT},
                ST <: NonlinearSolver{DT}} <: StochasticIntegratorPRK{DT,TT,D,M,S}
    params::PT
    solver::ST
    caches::CacheDict{PT}

    function IntegratorSIPRK(params::ParametersSIPRK{DT,TT,D,M,S}, solver::ST, caches) where {DT,TT,D,M,S,ST}
        new{DT, TT, D, M, S, typeof(params), ST}(params, solver, caches)
    end

    function IntegratorSIPRK{DT,D,M}(equations::NamedTuple, tableau::TableauSIPRK{TT}, Δt::TT; K::Int=0) where {DT,TT,D,M}
        # get number of stages
        S = tableau.s

        # K - the integer in the bound A = √(2 K Δt |log Δt|) due to Milstein & Tretyakov; K=0 no truncation
        K==0 ? A = 0.0 : A = sqrt( 2*K*Δt*abs(log(Δt)) )

        # create params
        params = ParametersSIPRK{DT,D,M}(equations, tableau, Δt, zeros(DT,M), zeros(DT,M), A)

        # create cache dict
        caches = CacheDict(params)

        # create solver
        solver = create_nonlinear_solver(DT, 2*D*S, params, caches)

        # create integrator
        IntegratorSIPRK(params, solver, caches)
    end

    function IntegratorSIPRK(equation::PSDE{DT,TT}, tableau::TableauSIPRK{TT}, Δt::TT; kwargs...) where {DT,TT}
        IntegratorSIPRK{DT, ndims(equation), equation.m}(get_function_tuple(equation), tableau, Δt; kwargs...)
    end
end


@doc raw"""
This function computes initial guesses for `Y`, `Z` and assigns them to int.solver.x
The prediction is calculated using an explicit integrator.

SIMPLE SOLUTION
The simplest initial guess for `Y`, `Z` is 0:
`int.solver.x .= zeros(eltype(int), 2*tableau(int).s*ndims(int))`

USING AN EXPLICIT INTEGRATOR TO COMPUTE AN INITIAL GUESS
Below we use the R2 method of Burrage & Burrage to calculate
the internal stages at the times c[1]...c[s].
This approach seems to give very good approximations if the time step
and magnitude of noise are not too large. If the noise intensity is too big,
one may have to perform a few iterations of the explicit method with a smaller
time step, use a higher-order explicit method (e.g. CL or G5), or use
the simple solution above.

When calling this function, int.params should contain the data:
`int.params.q` - the q solution at the previous time step
`int.params.p` - the p solution at the previous time step
`int.params.t` - the time of the previous step
`int.params.ΔW`- the increment of the Brownian motion for the current step
"""
function initial_guess!(int::IntegratorSIPRK{DT,TT}, sol::AtomicSolutionPSDE{DT,TT},
                        cache::IntegratorCacheSIPRK{DT}=int.caches[DT]) where {DT,TT}

    local t2::TT
    local Δt_local::TT

    # Evaluating the functions v and B at t,q - same for all stages
    int.params.equ[:v](int.params.t, int.params.q, int.params.p, cache.V1)
    int.params.equ[:B](int.params.t, int.params.q, int.params.p, cache.B1)
    int.params.equ[:f](int.params.t, int.params.q, int.params.p, cache.F1)
    int.params.equ[:G](int.params.t, int.params.q, int.params.p, cache.G1)

    # Calculating the positions q at the points qdrift.c[i]
    # if qdrift.c==pdrift.c, then also calculating the momenta
    for i in eachstage(int)
        # Taking the c[i] from the qdrift tableau.
        Δt_local = tableau(int).qdrift.c[i]  * timestep(int)
        cache.Δw .= tableau(int).qdrift.c[i] .* int.params.ΔW

        mul!(cache.ΔQ, cache.B1, cache.Δw)
        mul!(cache.ΔP, cache.G1, cache.Δw)
        @. cache.Q[i] = int.params.q + 2. / 3. * Δt_local * cache.V1 + 2. / 3. * cache.ΔQ
        @. cache.P[i] = int.params.p + 2. / 3. * Δt_local * cache.F1 + 2. / 3. * cache.ΔP

        t2 = int.params.t + 2. / 3. * Δt_local

        int.params.equ.v(t2, cache.Q[i], cache.P[i], cache.V2)
        int.params.equ.B(t2, cache.Q[i], cache.P[i], cache.B2)
        int.params.equ.f(t2, cache.Q[i], cache.P[i], cache.F2)
        int.params.equ.G(t2, cache.Q[i], cache.P[i], cache.G2)

        #Calculating the Y's and assigning them to the array int.solver.x as initial guesses
        mul!(cache.Y1, cache.B1, cache.Δw)
        mul!(cache.Y2, cache.B2, cache.Δw)
        for j in eachdim(int)
            int.solver.x[(i-1)*ndims(int)+j] =  Δt_local*(1. / 4. * cache.V1[j] + 3. / 4. * cache.V2[j]) + 1. / 4. * cache.Y1[j] + 3. / 4. * cache.Y2[j]
        end

        # if the collocation points are the same for both q and p parts
        mul!(cache.Z1, cache.G1, cache.Δw)
        mul!(cache.Z2, cache.G2, cache.Δw)
        if tableau(int).qdrift.c==tableau(int).pdrift.c
            for j in 1:ndims(int)
                int.solver.x[(tableau(int).s+i-1)*ndims(int)+j] =  Δt_local*(1. / 4. * cache.F1[j] + 3. / 4. * cache.F2[j]) + 1. / 4. * cache.Z1[j] + 3. / 4. * cache.Z2[j]
            end
        end
    end

    # If qdrift.c != pdrift.c, then calculating the momenta p at the points pdrift.c[i]
    if tableau(int).qdrift.c != tableau(int).pdrift.c
        for i in eachstage(int)
            # Taking the c[i] from the pdrift tableau.
            Δt_local = tableau(int).pdrift.c[i]  * timestep(int)
            cache.Δw .= tableau(int).pdrift.c[i] .* int.params.ΔW

            mul!(cache.ΔQ, cache.B1, cache.Δw)
            mul!(cache.ΔP, cache.G1, cache.Δw)
            @. cache.Q[i] = int.params.q + 2. / 3. * Δt_local * cache.V1 + 2. / 3. * cache.ΔQ
            @. cache.P[i] = int.params.p + 2. / 3. * Δt_local * cache.F1 + 2. / 3. * cache.ΔP

            t2 = int.params.t + 2. / 3. * Δt_local

            int.params.equ[:v](t2, cache.Q[i], cache.P[i], cache.V2)
            int.params.equ[:B](t2, cache.Q[i], cache.P[i], cache.B2)
            int.params.equ[:f](t2, cache.Q[i], cache.P[i], cache.F2)
            int.params.equ[:G](t2, cache.Q[i], cache.P[i], cache.G2)

            # Calculating the Z's and assigning them to the array int.solver.x as initial guesses
            # The guesses for the Y's have already been written to x above
            mul!(cache.Z1, cache.G1, cache.Δw)
            mul!(cache.Z2, cache.G2, cache.Δw)
            for j in eachdim(int)
                int.solver.x[(tableau(int).s+i-1)*ndims(int)+j] =  Δt_local*(1. / 4. * cache.F1[j] + 3. / 4. * cache.F2[j]) + 1. / 4. * cache.Z1[j] + 3. / 4. * cache.Z2[j]
            end
        end
    end
end


@doc raw"""
Unpacks the data stored in `x = (Y[1][1], Y[1][2], ... Y[1][D], Y[2][1], ..., Z[1][1], Z[1][2], ... Z[1][D], Z[2][1], ...)`
into `Y`, `Z`, calculates the internal stages `Q`, `P`, the values of the RHS
of the SDE ( `v(Q,P)`, `f(Q,P)`, `B(Q,P)` and `G(Q,P)` ), and assigns them to `V`, `F`, `B` and `G`.
Unlike for FIRK, here
`Y = Δt a_drift v(Q,P) + a_diff B(Q,P) ΔW`,
`Z = Δt â_drift v(Q,P) + â_diff B(Q,P) ΔW`.
"""
function compute_stages!(x::Vector{ST}, Q::Vector{Vector{ST}}, P::Vector{Vector{ST}},
                                        V::Vector{Vector{ST}}, F::Vector{Vector{ST}},
                                        B::Vector{Matrix{ST}}, G::Vector{Matrix{ST}},
                                        Y::Vector{Vector{ST}}, Z::Vector{Vector{ST}},
                                        params::ParametersSIPRK{DT,TT,D,M,S}) where {ST,DT,TT,D,M,S}

    local tqᵢ::TT       #times for the q internal stages
    local tpᵢ::TT       #times for the p internal stages


    @assert S == length(Q) == length(V) == length(B) == length(P) == length(F) == length(G)

    # copy x to Y, Z and calculate Q, P
    for i in 1:S
        @assert D == length(Q[i]) == length(V[i]) == size(B[i],1) == length(P[i]) == length(F[i]) == size(G[i],1)
        @assert M == size(B[i],2) == size(G[i],2)

        for k in 1:D
            Y[i][k] = x[D*(  i-1)+k]
            Z[i][k] = x[D*(S+i-1)+k]
        end

        Q[i] .= params.q .+ Y[i]
        P[i] .= params.p .+ Z[i]
    end

    # compute V=v(t,Q,P), F=f(t,Q,P), B=B(t,Q,P) and G=G(t,Q,P)
    for i in 1:S
        tqᵢ = params.t + params.Δt * params.tab.qdrift.c[i]
        tpᵢ = params.t + params.Δt * params.tab.pdrift.c[i]
        # calculates v(t,Q,P) and f(t,Q,P) and assigns to the i-th column of V and F
        params.equ[:v](tqᵢ, Q[i], P[i], V[i])
        params.equ[:f](tpᵢ, Q[i], P[i], F[i])
        # calculates B(t,Q,P) and assigns to the matrices B[i] and G[i]
        params.equ[:B](tqᵢ, Q[i], P[i], B[i])
        params.equ[:G](tpᵢ, Q[i], P[i], G[i])
    end
end

"Compute stages of implicit Runge-Kutta methods."
function Integrators.function_stages!(x::Vector{ST}, b::Vector{ST},
                params::ParametersSIPRK{DT,TT,D,M,S},
                caches::CacheDict) where {ST,DT,TT,D,M,S}

    # get cache for internal stages
    cache = caches[ST]

    compute_stages!(x, cache.Q, cache.P, cache.V, cache.F, cache.B, cache.G, cache.Y, cache.Z, params)

    local y1::ST
    local y2::ST
    local z1::ST
    local z2::ST

    # compute b = - (Y-AV)
    for i in 1:S
        for k in 1:D
            y1 = 0
            y2 = 0
            z1 = 0
            z2 = 0
            for j in 1:S
                y1 += params.tab.qdrift.a[i,j] * cache.V[j][k] * params.Δt
                y2 += params.tab.qdrift.â[i,j] * cache.V[j][k] * params.Δt
                z1 += params.tab.pdrift.a[i,j] * cache.F[j][k] * params.Δt
                z2 += params.tab.pdrift.â[i,j] * cache.F[j][k] * params.Δt
                for l in 1:M
                    y1 += params.tab.qdiff.a[i,j] * cache.B[j][k,l] * params.ΔW[l]
                    y2 += params.tab.qdiff.â[i,j] * cache.B[j][k,l] * params.ΔW[l]
                    z1 += params.tab.pdiff.a[i,j] * cache.G[j][k,l] * params.ΔW[l]
                    z2 += params.tab.pdiff.â[i,j] * cache.G[j][k,l] * params.ΔW[l]
                end
            end
            b[D*(  i-1)+k] = - cache.Y[i][k] + (y1 + y2)
            b[D*(S+i-1)+k] = - cache.Z[i][k] + (z1 + z2)
        end
    end
end


function update_solution!(sol::AtomicSolutionPSDE{T}, tab::TableauSIPRK{T}, Δt::T, ΔW::Vector{T},
                          cache::IntegratorCacheSIPRK{T}) where {T}

    update_solution!(sol, cache.V, cache.F, cache.B, cache.G,
                     tab.qdrift.b, tab.qdiff.b,
                     tab.pdrift.b, tab.pdiff.b,
                     Δt, ΔW, cache.Δy, cache.Δz)

    update_solution!(sol, cache.V, cache.F, cache.B, cache.G,
                     tab.qdrift.b̂, tab.qdiff.b̂,
                     tab.pdrift.b̂, tab.pdiff.b̂,
                     Δt, ΔW, cache.Δy, cache.Δz)
end


"Integrate PSDE with a stochastic implicit partitioned Runge-Kutta integrator."
function Integrators.integrate_step!(int::IntegratorSIPRK{DT,TT}, sol::AtomicSolutionPSDE{DT,TT},
                                     cache::IntegratorCacheSIPRK{DT}=int.caches[DT]) where {DT,TT}
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
    compute_stages!(int.solver.x, cache.Q, cache.P, cache.V, cache.F, cache.B, cache.G, cache.Y, cache.Z, int.params)

    # compute final update
    update_solution!(sol, tableau(int), timestep(int), int.params.ΔW, cache)
end
