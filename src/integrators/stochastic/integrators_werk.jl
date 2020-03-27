@doc raw"""
Holds the tableau of a weak explicit Runge-Kutta method.

    Reference: Andreas Rossler, "Second order Runge-Kutta methods for Stratonovich stochastic differential equations",
    BIT Numerical Mathematics (2007) 47, equation (5.1).

Order of the tableau is not included, because unlike in the deterministic
setting, it depends on the properties of the noise (e.g., the dimension of
the Wiener process and the commutativity properties of the diffusion matrix)

Orders stored in `qdrift`, `qdiff` are irrelevant and set to 0.

`qdrift0`, `qdrift1`, `qdrift2` correspond to A0, A1, A2 in the paper
`qdiff0`, `qdiff1`, `qdiff2`, `qdiff3` correspond to B0, B1, B2, B3
`qdrift0.b = alpha`
`qdiff0.b  = beta1`
`qdiff3.b  = beta2`
`qdrift0.c = c0`
`qdrift1.c = c1`
`qdrift2.c = c2`
"""
struct TableauWERK{T} <: AbstractTableauERK{T}
    name::Symbol
    s::Int

    qdrift0::CoefficientsRK{T}
    qdrift1::CoefficientsRK{T}
    qdrift2::CoefficientsRK{T}

    qdiff0::CoefficientsRK{T}
    qdiff1::CoefficientsRK{T}
    qdiff2::CoefficientsRK{T}
    qdiff3::CoefficientsRK{T}

    function TableauWERK{T}(name, qdrift0, qdrift1, qdrift2, qdiff0, qdiff1, qdiff2, qdiff3) where {T}
        @assert qdrift0.s == qdrift1.s == qdrift2.s == qdiff0.s == qdiff1.s == qdiff2.s == qdiff3.s
        @assert qdrift0.c[1] == qdrift1.c[1] == qdrift2.c[1] == qdiff0.c[1] == qdiff1.c[1] == qdiff2.c[1] == qdiff3.c[1] == 0.
        @assert istrilstrict(qdrift0.a)
        @assert istrilstrict(qdrift1.a)
        @assert istrilstrict(qdrift2.a)
        @assert istrilstrict(qdiff0.a)
        @assert istrilstrict(qdiff1.a)
        @assert istrilstrict(qdiff2.a)
        @assert istrilstrict(qdiff3.a)
        @assert !(qdrift0.s==1. && qdrift0.a[1,1] ≠ 0.)
        @assert !(qdrift1.s==1. && qdrift1.a[1,1] ≠ 0.)
        @assert !(qdrift2.s==1. && qdrift2.a[1,1] ≠ 0.)
        @assert !(qdiff0.s==1. && qdiff0.a[1,1] ≠ 0.)
        @assert !(qdiff1.s==1. && qdiff1.a[1,1] ≠ 0.)
        @assert !(qdiff2.s==1. && qdiff2.a[1,1] ≠ 0.)
        @assert !(qdiff3.s==1. && qdiff3.a[1,1] ≠ 0.)
        new(name, qdrift0.s, qdrift0, qdrift1, qdrift2, qdiff0, qdiff1, qdiff2, qdiff3)
    end
end

function TableauWERK(name, qdrift0::CoefficientsRK{T}, qdrift1::CoefficientsRK{T}, qdrift2::CoefficientsRK{T},
                            qdiff0::CoefficientsRK{T}, qdiff1::CoefficientsRK{T}, qdiff2::CoefficientsRK{T}, qdiff3::CoefficientsRK{T}) where {T}
    TableauWERK{T}(name, qdrift0, qdrift1, qdrift2, qdiff0, qdiff1, qdiff2, qdiff3)
end

function TableauWERK(name::Symbol, A0::Matrix{T}, A1::Matrix{T}, A2::Matrix{T},
                                   B0::Matrix{T}, B1::Matrix{T}, B2::Matrix{T}, B3::Matrix{T},
                                    α::Vector{T}, β1::Vector{T}, β2::Vector{T},
                                   c0::Vector{T}, c1::Vector{T}, c2::Vector{T} ) where {T}
    TableauWERK{T}(name, CoefficientsRK(name, 0, A0, α, c0), CoefficientsRK(name, 0, A1, α, c1), CoefficientsRK(name, 0, A2, α, c2),
                         CoefficientsRK(name, 0, B0, β1, c0), CoefficientsRK(name, 0, B1, β1, c1), CoefficientsRK(name, 0, B2, β1, c2), CoefficientsRK(name, 0, B3, β2, c1))
end


"Parameters for weak stochastic explicit Runge-Kutta methods."
struct ParametersWERK{DT, TT, D, M, S, ET <: NamedTuple} <: Parameters{DT,TT}
    equs::ET
    tab::TableauWERK{TT}
    Δt::TT

    function ParametersWERK{DT,D,M}(equs::ET, tab::TableauWERK{TT}, Δt::TT) where {DT, TT, D, M, ET <: NamedTuple}
        new{DT, TT, D, M, tab.s, ET}(equs, tab, Δt)
    end
end


"Weak Stochastic Explicit Runge-Kutta integrator cache."
struct IntegratorCacheWERK{DT,D,M,S} <: SDEIntegratorCache{DT,D,M}
    Δy::Vector{DT}

    Q0::Vector{DT}           # Q0[1..n]       - the internal stage H^(0)_i for a given i
    Q1::Vector{Vector{DT}}   # Q1[1..m][1..n] - the internal stages H^(l)_i for a given i
    Q2::Vector{Vector{DT}}   # Q2[1..m][1..n] - the internal stages \hat H^(l)_i for a given i
    V ::Vector{Vector{DT}}
    B1::Vector{Matrix{DT}}   # B1[i] holds the values of the diffusion matrix such that B1[i][:,l] is evaluated at H^(l)_i
    B2::Vector{Matrix{DT}}   # B2[i] holds the values of the diffusion matrix such that B2[i][:,l] is evaluated at \hat H^(l)_i

    tB::Vector{DT}           # the value of the l-th column of the diffusion term evaluated at the internal stage H^(l)_i for given i and l

    function IntegratorCacheWERK{DT,D,M,S}() where {DT,D,M,S}
        # create internal stage vectors
        Q0 = zeros(DT, D)
        Q1 = create_internal_stage_vector(DT, D, M)
        Q2 = create_internal_stage_vector(DT, D, M)
        V  = create_internal_stage_vector(DT, D, S)
        B1 = create_internal_stage_matrix(DT, D, M, S)
        B2 = create_internal_stage_matrix(DT, D, M, S)

        new(zeros(DT,M), Q0, Q1, Q2, V, B1, B2, zeros(DT,D))
    end
end

function Integrators.IntegratorCache{ST}(params::ParametersWERK{DT,TT,D,M,S}; kwargs...) where {ST,DT,TT,D,M,S}
    IntegratorCacheWERK{ST,D,M,S}(; kwargs...)
end

@inline Integrators.CacheType(ST, params::ParametersWERK{DT,TT,D,M,S}) where {DT,TT,D,M,S} = IntegratorCacheWERK{ST,D,M,S}


"Weak Stochastic Explicit Runge-Kutta integrator."
struct IntegratorWERK{DT, TT, D, M, S, ET} <: StochasticIntegratorRK{DT,TT,D,M,S}
    params::ParametersWERK{DT,TT,D,M,S,ET}
    caches::CacheDict{ParametersWERK{DT,TT,D,M,S,ET}}


    function IntegratorWERK(params::ParametersWERK{DT,TT,D,M,S,ET}, caches) where {DT,TT,D,M,S,ET}
        new{DT, TT, D, M, S, ET}(params, caches)
    end

    function IntegratorWERK{DT,D,M}(equations::ET, tableau::TableauWERK{TT}, Δt::TT) where {DT, TT, D, M, ET <: NamedTuple}
        # get number of stages
        S = tableau.s

        # create params
        params = ParametersWERK{DT,D,M}(equations, tableau, Δt)

        # create cache dict
        caches = CacheDict(params)

        # create integrator
        IntegratorWERK(params, caches)
    end

    function IntegratorWERK(equation::SDE{DT,TT}, tableau::TableauWERK{TT}, Δt::TT; kwargs...) where {DT,TT}
        IntegratorWERK{DT, ndims(equation), equation.m}(get_function_tuple(equation), tableau, Δt; kwargs...)
    end
end


"Integrate SDE with explicit Runge-Kutta integrator."
function Integrators.integrate_step!(int::IntegratorWERK{DT,TT}, sol::AtomicSolutionSDE{DT,TT},
                                     cache::IntegratorCacheWERK{DT}=int.caches[DT]) where {DT,TT}
    local tᵢ::TT
    local ydrift::DT
    local ydiff2::DT

    # reset cache
    reset!(sol, timestep(int))

    # THE FIRST INTERNAL STAGES ARE ALL EQUAL TO THE PREVIOUS STEP SOLUTION
    # calculates v(t,tQ0) and assigns to the 1st column of V
    equation(int, :v)(sol.t̅, sol.q̅, cache.V[1])
    # calculates B(t,Q) and assigns to the matrix B1[1]
    # no need to calculate the columns of B separately, because all internal stages
    # for i=1 are equal to int.q[r,m]
    equation(int, :B)(sol.t̅, sol.q̅, cache.B1[1])
    # copy B1[i] to B2[i] (they're equal for i=1)
    cache.B2[1] .= cache.B1[1]

    for i in 2:tableau(int).s
        # Calculating the internal stage H^(0)_i
        for k in eachindex(cache.Q0)
            # contribution from the drift part
            ydrift = 0
            for j = 1:i-1
                ydrift += tableau(int).qdrift0.a[i,j] * cache.V[j][k]
            end

            # ΔW contribution from the diffusion part
            cache.Δy .= 0
            for j = 1:i-1
                for l in eachindex(cache.Δy)
                    cache.Δy[l] += tableau(int).qdiff0.a[i,j] * cache.B1[j][k,l]
                end
            end

            cache.Q0[k] = sol.q̅[k] + timestep(int) * ydrift + dot(cache.Δy, sol.ΔW)
        end

        # Calculating the internal stages H^(l)_i for l=1..sol.nm
        for k in eachindex(cache.Q1[1])
            # contribution from the drift part (same for all noises)
            ydrift = 0
            for j = 1:i-1
                ydrift += tableau(int).qdrift1.a[i,j] * cache.V[j][k]
            end

            # contribution from the diffusion part is different for each noise
            for noise_idx in eachindex(cache.Q1)

                # ΔW contribution from the diffusion part
                cache.Δy .= 0
                for j = 1:i-1
                    for l in eachindex(cache.Δy)
                        if l==noise_idx
                            cache.Δy[l] += tableau(int).qdiff1.a[i,j] * cache.B1[j][k,l]
                        else
                            cache.Δy[l] += tableau(int).qdiff3.a[i,j] * cache.B1[j][k,l]
                        end
                    end
                end

                cache.Q1[noise_idx][k] = sol.q̅[k] + timestep(int) * ydrift + dot(cache.Δy, sol.ΔW)
            end
        end

        # Calculating the internal stages \hat H^(l)_i for l=1..sol.nm
        for k in eachindex(cache.Q2[1])
            # contribution from the drift part (same for all noises)
            ydrift = 0
            for j = 1:i-1
                ydrift += tableau(int).qdrift2.a[i,j] * cache.V[j][k]
            end

            # contribution from the diffusion part is different for each noise
            for noise_idx in eachindex(cache.Q2)

                # ΔW contribution from the diffusion part
                ydiff2 = 0

                for j = 1:i-1
                    for l in eachindex(sol.ΔW, sol.ΔZ)
                        # calculating the terms with the random variables \hat I_{noise_idx,l}
                        if l<noise_idx
                            ydiff2 += tableau(int).qdiff2.a[i,j] * cache.B1[j][k,l] * sol.ΔW[noise_idx] * sol.ΔZ[l]
                        elseif l>noise_idx
                            ydiff2 += -tableau(int).qdiff2.a[i,j] * cache.B1[j][k,l] * sol.ΔW[l] * sol.ΔZ[noise_idx]
                        end
                    end
                end

                cache.Q2[noise_idx][k] = sol.q̅[k] + timestep(int) * ydrift + ydiff2/sqrt(timestep(int))
            end
        end

        # CALCULATING THE NEW VALUES OF V
        tᵢ = sol.t̅ + timestep(int) * tableau(int).qdrift0.c[i]
        equation(int, :v)(tᵢ, cache.Q0, cache.V[i])

        # CALCULATING THE NEW VALUES OF B1
        # each column of B evaluated at a different internal stage
        tᵢ = sol.t̅ + timestep(int) * tableau(int).qdrift1.c[i]

        for l in eachnoise(int)
            equation(int, :B)(tᵢ, cache.Q1[l], cache.B1[i], l)
        end

        # CALCULATING THE NEW VALUES OF B2
        # each column of B evaluated at a different internal stage
        tᵢ = sol.t̅ + timestep(int) * tableau(int).qdrift2.c[i]

        for l in eachnoise(int)
            equation(int, :B)(tᵢ, cache.Q2[l], cache.B2[i], l)
        end
    end

    # compute final update
    update_solution!(sol, cache.V, cache.B1, cache.B2, tableau(int).qdrift0.b, tableau(int).qdiff0.b, tableau(int).qdiff3.b, timestep(int), sol.ΔW, cache.Δy)
    # update_solution!(sol, cache.V, cache.B1, cache.B2, tableau(int).qdrift0.b̂, tableau(int).qdiff0.b̂, tableau(int).qdiff3.b̂, timestep(int), sol.ΔW, cache.Δy)
end
