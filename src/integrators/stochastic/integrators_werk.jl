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



"Stochastic Explicit Runge-Kutta integrator."
struct IntegratorWERK{DT, TT, ET <: SDE{DT,TT}} <: StochasticIntegrator{DT,TT}
    equation::ET
    tableau::TableauWERK{TT}
    Δt::TT

    Δy::Vector{DT}

    Q0::Vector{DT}           # Q0[1..n]      - the internal stage H^(0)_i for a given i
    Q1::Vector{Vector{DT}}   # Q1[1..m][1..n] - the internal stages H^(l)_i for a given i
    Q2::Vector{Vector{DT}}   # Q2[1..m][1..n] - the internal stages \hat H^(l)_i for a given i
    V ::Vector{Vector{DT}}
    B1::Vector{Matrix{DT}}   # B1[i] holds the values of the diffusion matrix such that B1[i][:,l] is evaluated at H^(l)_i
    B2::Vector{Matrix{DT}}   # B2[i] holds the values of the diffusion matrix such that B2[i][:,l] is evaluated at \hat H^(l)_i
    tB::Vector{DT}           # the value of the l-th column of the diffusion term evaluated at the internal stage H^(l)_i for given i and l

    function IntegratorWERK{DT,TT}(equation::ET, tableau, Δt::TT) where {DT, TT, ET <: SDE{DT,TT}}
        D = equation.d
        M = equation.m
        NS= max(equation.ns,equation.ni)
        S = tableau.s

        # create internal stage vectors
        Q0 = zeros(DT, D)
        Q1 = create_internal_stage_vector(DT, D, M)
        Q2 = create_internal_stage_vector(DT, D, M)
        V  = create_internal_stage_vector(DT, D, S)
        B1 = create_internal_stage_vector(DT, D, M, S)
        B2 = create_internal_stage_vector(DT, D, M, S)

        new{DT,TT,ET}(equation, tableau, Δt, zeros(DT,M), Q0, Q1, Q2, V, B1, B2, zeros(DT,D))
    end
end

function IntegratorWERK(equation::SDE{DT,TT}, tableau::TableauWERK{TT}, Δt::TT) where {DT,TT}
    IntegratorWERK{DT,TT}(equation, tableau, Δt)
end

@inline equation(integrator::IntegratorWERK) = integrator.equation
@inline timestep(integrator::IntegratorWERK) = integrator.Δt
@inline tableau(integrator::IntegratorWERK)  = integrator.tableau
@inline nstages(integrator::IntegratorWERK)  = nstages(tableau(integrator))
@inline eachstage(integrator::IntegratorWERK) = 1:nstages(integrator)
@inline Base.eltype(integrator::IntegratorWERK{DT}) where {DT} = DT


"Integrate SDE with explicit Runge-Kutta integrator."
function Integrators.integrate_step!(int::IntegratorWERK{DT,TT}, sol::AtomicSolutionSDE{DT,TT}) where {DT,TT}
    local tᵢ::TT
    local ydrift::DT
    local ydiff2::DT

    # reset cache
    reset!(sol, timestep(int))

    # THE FIRST INTERNAL STAGES ARE ALL EQUAL TO THE PREVIOUS STEP SOLUTION
    # calculates v(t,tQ0) and assigns to the 1st column of V
    int.equation.v(sol.t̅, sol.q̅, int.V[1])
    # calculates B(t,Q) and assigns to the matrix B1[1]
    # no need to calculate the columns of B separately, because all internal stages
    # for i=1 are equal to int.q[r,m]
    int.equation.B(sol.t̅, sol.q̅, int.B1[1])
    # copy B1[i] to B2[i] (they're equal for i=1)
    int.B2[1] .= int.B1[1]

    for i in 2:int.tableau.s

        # Calculating the internal stage H^(0)_i
        @inbounds for k in eachindex(int.Q0)
            # contribution from the drift part
            ydrift = 0
            for j = 1:i-1
                ydrift += int.tableau.qdrift0.a[i,j] * int.V[j][k]
            end

            # ΔW contribution from the diffusion part
            int.Δy .= 0
            for j = 1:i-1
                for l in eachindex(int.Δy)
                    int.Δy[l] += int.tableau.qdiff0.a[i,j] * int.B1[j][k,l]
                end
            end

            int.Q0[k] = sol.q̅[k] + int.Δt * ydrift + dot(int.Δy,sol.ΔW)
        end

        # Calculating the internal stages H^(l)_i for l=1..sol.nm
        @inbounds for k in eachindex(int.Q1[1])
            # contribution from the drift part (same for all noises)
            ydrift = 0
            for j = 1:i-1
                ydrift += int.tableau.qdrift1.a[i,j] * int.V[j][k]
            end

            # contribution from the diffusion part is different for each noise
            @inbounds for noise_idx in eachindex(int.Q1)

                # ΔW contribution from the diffusion part
                int.Δy .= 0
                for j = 1:i-1
                    for l in eachindex(int.Δy)
                        if l==noise_idx
                            int.Δy[l] += int.tableau.qdiff1.a[i,j] * int.B1[j][k,l]
                        else
                            int.Δy[l] += int.tableau.qdiff3.a[i,j] * int.B1[j][k,l]
                        end
                    end
                end

                int.Q1[noise_idx][k] = sol.q̅[k] + int.Δt * ydrift + dot(int.Δy,sol.ΔW)
            end
        end

        # Calculating the internal stages \hat H^(l)_i for l=1..sol.nm
        @inbounds for k in eachindex(int.Q2[1])
            # contribution from the drift part (same for all noises)
            ydrift = 0
            for j = 1:i-1
                ydrift += int.tableau.qdrift2.a[i,j] * int.V[j][k]
            end

            # contribution from the diffusion part is different for each noise
            @inbounds for noise_idx in eachindex(int.Q2)

                # ΔW contribution from the diffusion part
                ydiff2 = 0

                for j = 1:i-1
                    for l in eachindex(sol.ΔW, sol.ΔZ)
                        # calculating the terms with the random variables \hat I_{noise_idx,l}
                        if l<noise_idx
                            ydiff2 += int.tableau.qdiff2.a[i,j] * int.B1[j][k,l] * sol.ΔW[noise_idx] * sol.ΔZ[l]
                        elseif l>noise_idx
                            ydiff2 += -int.tableau.qdiff2.a[i,j] * int.B1[j][k,l] * sol.ΔW[l] * sol.ΔZ[noise_idx]
                        end
                    end
                end

                int.Q2[noise_idx][k] = sol.q̅[k] + int.Δt * ydrift + ydiff2/sqrt(int.Δt)
            end
        end

        # CALCULATING THE NEW VALUES OF V
        tᵢ = sol.t̅ + int.Δt * int.tableau.qdrift0.c[i]
        int.equation.v(tᵢ, int.Q0, int.V[i])

        # CALCULATING THE NEW VALUES OF B1
        # each column of B evaluated at a different internal stage
        tᵢ = sol.t̅ + int.Δt * int.tableau.qdrift1.c[i]

        for l = 1:equation(int).m
            int.equation.B(tᵢ, int.Q1[l], int.B1[i], l)
        end

        # CALCULATING THE NEW VALUES OF B2
        # each column of B evaluated at a different internal stage
        tᵢ = sol.t̅ + int.Δt * int.tableau.qdrift2.c[i]

        for l = 1:equation(int).m
            int.equation.B(tᵢ, int.Q2[l], int.B2[i], l)
        end

    end

    # compute final update
    update_solution!(sol, int.V, int.B1, int.B2, int.tableau.qdrift0.b, int.tableau.qdiff0.b, int.tableau.qdiff3.b, int.Δt, sol.ΔW, int.Δy)
    # update_solution!(sol, int.V, int.B1, int.B2, int.tableau.qdrift0.b̂, int.tableau.qdiff0.b̂, int.tableau.qdiff3.b̂, int.Δt, sol.ΔW, int.Δy)
end
