@doc raw"""
Holds the tableau of a stochastic implicit Runge-Kutta method.

`qdrift` holds the RK coefficients for the drift part,
`qdiff` holds the RK coefficients for the diffusion part of the SDE.

Order of the tableau is not included, because unlike in the deterministic
setting, it depends on the properties of the noise (e.g., the dimension of
the Wiener process and the commutativity properties of the diffusion matrix).

Orders stored in `qdrift` and `qdiff` are understood as the classical orders of these methods.
"""
struct TableauSIRK{T} <: AbstractTableauIRK{T}
    name::Symbol
    s::Int
    qdrift::CoefficientsRK{T}
    qdiff::CoefficientsRK{T}

    function TableauSIRK{T}(name, s, qdrift, qdiff) where {T}
        @assert s == qdrift.s == qdiff.s
        new(name, s, qdrift, qdiff)
    end
end

function TableauSIRK(name::Symbol, qdrift::CoefficientsRK{T}, qdiff::CoefficientsRK{T}) where {T}
    TableauSIRK{T}(name, qdrift.s, qdrift, qdiff)
end

function TableauSIRK(name::Symbol, order_drift::Int, a_drift::Matrix{T}, b_drift::Vector{T}, c_drift::Vector{T}, order_diff::Int, a_diff::Matrix{T}, b_diff::Vector{T}, c_diff::Vector{T}) where {T}
    TableauSIRK{T}(name, length(c_drift), CoefficientsRK(name, order_drift, a_drift, b_drift, c_drift), CoefficientsRK(name, order_diff, a_diff, b_diff, c_diff))
end

# TODO function readTableauSIRKFromFile(dir::AbstractString, name::AbstractString)


@doc raw"""
Parameters for right-hand side function of implicit Runge-Kutta methods.

`A` - if positive, the upper bound of the Wiener process increments; if `A=0.0`, no truncation
"""
mutable struct ParametersSIRK{DT, TT, D, M, S, ET <: NamedTuple} <: Parameters{DT,TT}
    equ::ET
    tab::TableauSIRK{TT}
    Δt::TT

    ΔW::Vector{DT}
    ΔZ::Vector{DT}
    A::DT

    t::TT
    q::Vector{DT}

    function ParametersSIRK{DT,D,M}(equ::ET, tab::TableauSIRK{TT}, Δt::TT, ΔW::Vector{DT}, ΔZ::Vector{DT}, A::DT) where {DT, TT, D, M, ET <: NamedTuple}
        @assert M == length(ΔW) == length(ΔZ)
        new{DT, TT, D, M, tab.s, ET}(equ, tab, Δt, ΔW, ΔZ, A, zero(TT), zeros(DT,D))
    end
end


@doc raw"""
Structure for holding the internal stages `Q`, the values of the drift vector
and the diffusion matrix evaluated at the internal stages `V=v(Q)`, `B=B(Q)`,
and the increments `Y = a_drift*v(Q)*Δt + a_diff*B(Q)*ΔW`.
"""
struct IntegratorCacheSIRK{DT,D,M,S} <: SDEIntegratorCache{DT,D,M}
    Q::Vector{Vector{DT}}
    V::Vector{Vector{DT}}
    B::Vector{Matrix{DT}}
    Y::Vector{Vector{DT}}

    v::Vector{DT}
    b::Matrix{DT}
    y::Vector{DT}

    ΔQ::Vector{DT}
    Y1::Vector{DT}
    Y2::Vector{DT}
    V1::Vector{DT}
    V2::Vector{DT}
    B1::Matrix{DT}
    B2::Matrix{DT}
    Δw::Vector{DT}
    Δy::Vector{DT}

    function IntegratorCacheSIRK{DT,D,M,S}() where {DT,D,M,S}
        # create internal stage vectors
        Q = create_internal_stage_vector(DT, D, S)
        V = create_internal_stage_vector(DT, D, S)
        B = create_internal_stage_matrix(DT, D, M, S)
        Y = create_internal_stage_vector(DT, D, S)

        # create velocity and update vector
        v = zeros(DT,D)
        b = zeros(DT,D,M)
        y = zeros(DT,D)

        # create temporary vectors
        ΔQ = zeros(DT,D)
        Y1 = zeros(DT,D)
        Y2 = zeros(DT,D)
        V1 = zeros(DT,D)
        V2 = zeros(DT,D)
        B1 = zeros(DT,D,M)
        B2 = zeros(DT,D,M)
        Δw = zeros(DT,M)
        Δy = zeros(DT,M)

        new(Q, V, B, Y, v, b, y, ΔQ, Y1, Y2, V1, V2, B1, B2, Δw, Δy)
    end
end

function Integrators.IntegratorCache{ST}(params::ParametersSIRK{DT,TT,D,M,S}; kwargs...) where {ST,DT,TT,D,M,S}
    IntegratorCacheSIRK{ST,D,M,S}(; kwargs...)
end

@inline Integrators.CacheType(ST, params::ParametersSIRK{DT,TT,D,M,S}) where {DT,TT,D,M,S} = IntegratorCacheSIRK{ST,D,M,S}


"Stochastic implicit Runge-Kutta integrator."
struct IntegratorSIRK{DT, TT, D, M, S,
                PT <: ParametersSIRK{DT,TT},
                ST <: NonlinearSolver{DT}} <: StochasticIntegratorRK{DT,TT,D,M,S}
    params::PT
    solver::ST
    caches::CacheDict{PT}

    function IntegratorSIRK(params::ParametersSIRK{DT,TT,D,M,S}, solver::ST, caches) where {DT,TT,D,M,S,ST}
        new{DT, TT, D, M, S, typeof(params), ST}(params, solver, caches)
    end

    function IntegratorSIRK{DT,D,M}(equations::NamedTuple, tableau::TableauSIRK{TT}, Δt::TT; K::Int=0) where {DT,TT,D,M}
        # get number of stages
        S = tableau.s

        # K - the integer in the bound A = √(2 K Δt |log Δt|) due to Milstein & Tretyakov; K=0 no truncation
        K==0 ? A = 0.0 : A = sqrt( 2*K*Δt*abs(log(Δt)) )

        # create params
        params = ParametersSIRK{DT,D,M}(equations, tableau, Δt, zeros(DT,M), zeros(DT,M), A)

        # create cache dict
        caches = CacheDict(params)

        # create solver
        solver = create_nonlinear_solver(DT, D*S, params, caches)

        # create integrator
        IntegratorSIRK(params, solver, caches)
    end

    function IntegratorSIRK(equation::SDE{DT,TT}, tableau::TableauSIRK{TT}, Δt::TT; kwargs...) where {DT,TT}
        IntegratorSIRK{DT, ndims(equation), equation.m}(get_function_tuple(equation), tableau, Δt; kwargs...)
    end
end


@doc raw"""
This function computes initial guesses for Y and assigns them to `int.solver.x`.
The prediction is calculated using an explicit integrator.

SIMPLE SOLUTION
The simplest initial guess for `Y` is 0:
`int.solver.x .= zeros(eltype(int), tableau(int).s*ndims(int))`

USING AN EXPLICIT INTEGRATOR TO COMPUTE AN INITIAL GUESS
Below we use the R2 method of Burrage & Burrage to calculate
the internal stages at the times `c[1]...c[s]`.
This approach seems to give very good approximations if the time step
and magnitude of noise are not too large. If the noise intensity is too big,
one may have to perform a few iterations of the explicit method with a smaller
time step, use a higher-order explicit method (e.g. CL or G5), or use
the simple solution above.

When calling this function, int.params should contain the data:
`int.params.q`  - the solution at the previous time step
`int.params.t`  - the time of the previous step
`int.params.ΔW` - the increment of the Brownian motion for the current step
"""
function initial_guess!(int::IntegratorSIRK{DT,TT}, sol::AtomicSolutionSDE{DT,TT},
                        cache::IntegratorCacheSIRK{DT}=int.caches[DT]) where {DT,TT}
    local t2::TT
    local Δt_local::TT

    # Evaluating the functions v and B at t,q - same for all stages
    int.params.equ[:v](int.params.t, int.params.q, cache.V1)
    int.params.equ[:B](int.params.t, int.params.q, cache.B1)

    for i in eachstage(int)
        Δt_local = tableau(int).qdrift.c[i]  * int.params.Δt
        cache.Δw .= tableau(int).qdrift.c[i] .* int.params.ΔW

        mul!(cache.ΔQ, cache.B1, cache.Δw)
        @. cache.Q[i] = int.params.q + 2. / 3. * Δt_local * cache.V1 + 2. / 3. * cache.ΔQ

        t2 = int.params.t + 2. / 3. * Δt_local

        int.params.equ[:v](t2, cache.Q[i], cache.V2)
        int.params.equ[:B](t2, cache.Q[i], cache.B2)

        # Calculating the Y's and assigning them to the array int.solver.x as initial guesses
        mul!(cache.Y1, cache.B1, cache.Δw)
        mul!(cache.Y2, cache.B2, cache.Δw)

        for j in eachdim(int)
            int.solver.x[(i-1)*ndims(int)+j] = Δt_local*(1. / 4. * cache.V1[j] + 3. / 4. * cache.V2[j]) + 1. / 4. * cache.Y1[j] + 3. / 4. * cache.Y2[j]
        end
    end
end


@doc raw"""
Unpacks the data stored in `x = (Y[1][1], Y[1][2], ... Y[1][D], Y[2][1], ...)`
into `Y`, calculates the internal stages `Q`, the values of the RHS
of the SDE ( v(Q) and B(Q) ), and assigns them to `V` and `B`.
Unlike for FIRK, here `Y = a v(Q) Δt + â B(Q) ΔW`
"""
function compute_stages!(x::Vector{ST}, Q::Vector{Vector{ST}}, V::Vector{Vector{ST}},
                         B::Vector{Matrix{ST}}, Y::Vector{Vector{ST}},
                         params::ParametersSIRK{DT,TT,D}) where {ST,DT,TT,D}

    local tᵢ::TT

    @assert length(Y) == length(Q) == length(V) == length(B)

    # copy x to Y and calculate Q
    for i in eachindex(Q,Y)
        @assert size(B[i],1) == length(Q[i]) == length(V[i])
        for k in eachindex(Q[i])
            Y[i][k] = x[D*(i-1)+k]
        end
        Q[i] .= params.q .+ Y[i]
    end

    # compute V = v(Q) and B=B(Q)
    for i in eachindex(Q,V,B)
        tᵢ = params.t + params.Δt * params.tab.qdrift.c[i]
        # calculates v(t,Q[i]) and assigns to the i-th column of V
        params.equ[:v](tᵢ, Q[i], V[i])
        # calculates B(t,Q[i]) and assigns to the matrix B[i]
        params.equ[:B](tᵢ, Q[i], B[i])
    end
end

"Compute stages of stochastic implicit Runge-Kutta methods."
function Integrators.function_stages!(x::Vector{ST}, b::Vector{ST},
                params::ParametersSIRK{DT,TT,D,M,S},
                caches::CacheDict) where {ST,DT,TT,D,M,S}

    # get cache for internal stages
    cache = caches[ST]

    compute_stages!(x, cache.Q, cache.V, cache.B, cache.Y, params)

    local y1::ST
    local y2::ST
    local y3::ST
    local y4::ST

    # compute b = - (Y-AV)
    for i in 1:S
        for k in 1:D
            y1 = 0
            y2 = 0
            y3 = 0
            y4 = 0
            for j in 1:S
                y1 += params.tab.qdrift.a[i,j] * cache.V[j][k] * params.Δt
                y2 += params.tab.qdrift.â[i,j] * cache.V[j][k] * params.Δt
                for l in 1:M
                    y3 += params.tab.qdiff.a[i,j] * cache.B[j][k,l] * params.ΔW[l]
                    y4 += params.tab.qdiff.â[i,j] * cache.B[j][k,l] * params.ΔW[l]
                end
            end
            b[D*(i-1)+k] = - cache.Y[i][k] + (y1 + y2) + (y3 + y4)
        end
    end
end


"Integrate SDE with a stochastic implicit Runge-Kutta integrator."
function Integrators.integrate_step!(int::IntegratorSIRK{DT,TT}, sol::AtomicSolutionSDE{DT,TT},
                                     cache::IntegratorCacheSIRK{DT}=int.caches[DT]) where {DT,TT}
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

    # compute final update
    update_solution!(sol, cache.V, cache.B, tableau(int).qdrift.b, tableau(int).qdiff.b, int.params.Δt, int.params.ΔW, cache.Δy)
    update_solution!(sol, cache.V, cache.B, tableau(int).qdrift.b̂, tableau(int).qdiff.b̂, int.params.Δt, int.params.ΔW, cache.Δy)
end
