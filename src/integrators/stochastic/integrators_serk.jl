
"Holds the tableau of a stochastic explicit Runge-Kutta method."
struct TableauSERK{T} <: AbstractTableauERK{T}
    name::Symbol
    s::Int

    qdrift::CoefficientsRK{T}
    qdiff::CoefficientsRK{T}
    qdiff2::CoefficientsRK{T}

    # Order of the tableau is not included, because unlike in the deterministic
    # setting, it depends on the properties of the noise (e.g., the dimension of
    # the Wiener process and the commutativity properties of the diffusion matrix)
    #
    # Orders stored in qdrift, qdiff and qdiff2 are understood as the classical orders of these methods.

    function TableauSERK{T}(name, qdrift, qdiff, qdiff2) where {T}
        @assert qdrift.s == qdiff.s == qdiff2.s
        @assert qdrift.c[1] == qdiff.c[1] == qdiff2.c[1] == 0.
        @assert istrilstrict(qdrift.a)
        @assert istrilstrict(qdiff.a)
        @assert istrilstrict(qdiff2.a)
        @assert !(qdrift.s==1. && qdrift.a[1,1] ≠ 0.)
        @assert !(qdiff.s==1. && qdiff.a[1,1] ≠ 0.)
        @assert !(qdiff2.s==1. && qdiff2.a[1,1] ≠ 0.)
        new(name, qdrift.s, qdrift, qdiff, qdiff2)
    end
end

function TableauSERK(name, qdrift::CoefficientsRK{T}, qdiff::CoefficientsRK{T}, qdiff2::CoefficientsRK{T}) where {T}
    TableauSERK{T}(name, qdrift, qdiff, qdiff2)
end

function TableauSERK(name, qdrift::CoefficientsRK{T}, qdiff::CoefficientsRK{T}) where {T}
    TableauSERK{T}(name, qdrift, qdiff, CoefficientsRK(T, :NULL, 0, zero(qdrift.a), zero(qdrift.b), zero(qdrift.c)) )
end

function TableauSERK(name::Symbol, order_drift::Int, a_drift::Matrix{T}, b_drift::Vector{T}, c_drift::Vector{T},
                                    order_diff::Int, a_diff::Matrix{T}, b_diff::Vector{T}, c_diff::Vector{T},
                                    order_diff2::Int, a_diff2::Matrix{T}, b_diff2::Vector{T}, c_diff2::Vector{T}) where {T}
    TableauSERK{T}(name, CoefficientsRK(name, order_drift, a_drift, b_drift, c_drift), CoefficientsRK(name, order_diff, a_diff, b_diff, c_diff),
                    CoefficientsRK(name, order_diff2, a_diff2, b_diff2, c_diff2))
end

function TableauSERK(name::Symbol, order_drift::Int, a_drift::Matrix{T}, b_drift::Vector{T}, c_drift::Vector{T},
                                    order_diff::Int, a_diff::Matrix{T}, b_diff::Vector{T}, c_diff::Vector{T}) where {T}
    TableauSERK{T}(name, CoefficientsRK(name, order_drift, a_drift, b_drift, c_drift), CoefficientsRK(name, order_diff, a_diff, b_diff, c_diff),
                    CoefficientsRK(:NULL, 0, zero(a_drift), zero(b_drift), zero(c_drift)))
end



"Stochastic Explicit Runge-Kutta integrator."
struct IntegratorSERK{DT, TT, ET <: SDE{DT,TT}} <: StochasticIntegrator{DT,TT}
    equation::ET
    tableau::TableauSERK{TT}
    Δt::TT

    Δy::Vector{DT}

    Q::Vector{Vector{DT}}     # Q[j][k] - the k-th component of the j-th internal stage
    V::Vector{Vector{DT}}     # V[j][k] - the k-th component of v(Q[j])
    B::Vector{Matrix{DT}}     # B[j]    - the diffusion matrix B(Q[j])

    function IntegratorSERK{DT,TT}(equation::ET, tableau, Δt::TT) where {DT, TT, ET <: SDE{DT,TT}}
        D = equation.d
        M = equation.m
        NS= max(equation.ns,equation.ni)
        S = tableau.s

        # create internal stage vectors
        Q = create_internal_stage_vector(DT, D, S)
        V = create_internal_stage_vector(DT, D, S)
        B = create_internal_stage_vector(DT, D, M, S)

        new{DT,TT,ET}(equation, tableau, Δt, zeros(DT,M), Q, V, B)
    end
end

function IntegratorSERK(equation::SDE{DT,TT}, tableau::TableauSERK{TT}, Δt::TT) where {DT,TT}
    IntegratorSERK{DT,TT}(equation, tableau, Δt)
end

@inline equation(integrator::IntegratorSERK) = integrator.equation
@inline timestep(integrator::IntegratorSERK) = integrator.Δt
@inline tableau(integrator::IntegratorSERK)  = integrator.tableau
@inline nstages(integrator::IntegratorSERK)  = nstages(tableau(integrator))
@inline eachstage(integrator::IntegratorSERK) = 1:nstages(integrator)
@inline Base.eltype(integrator::IntegratorSERK{DT}) where {DT} = DT


"""
Integrate SDE with explicit Runge-Kutta integrator.
  Calculating the n-th time step of the explicit integrator for the sample path m
"""
function integrate_step!(int::IntegratorSERK{DT,TT}, sol::AtomicSolutionSDE{DT,TT}) where {DT,TT}
    local tᵢ::TT
    local ydrift::DT

    # reset cache
    reset!(sol, timestep(int))

    # calculates v(t,tQ) and assigns to the i-th column of V
    int.equation.v(sol.t̅, sol.q̅, int.V[1])
    # calculates B(t,tQ) and assigns to the matrix BQ[1][:,:]
    int.equation.B(sol.t̅, sol.q̅, int.B[1])


    @inbounds for i in 2:int.tableau.s
        for k in eachindex(int.Q[i])
            # contribution from the drift part
            ydrift = 0
            for j = 1:i-1
                ydrift += int.tableau.qdrift.a[i,j] * int.V[j][k]
            end

            # ΔW contribution from the diffusion part
            int.Δy .= 0
            for j = 1:i-1
                for l = 1:equation(int).m
                    int.Δy[l] += int.tableau.qdiff.a[i,j] * int.B[j][k,l]
                end
            end

            int.Q[i][k] = sol.q̅[k] + int.Δt * ydrift + dot(int.Δy,sol.ΔW)

            # ΔZ contribution from the diffusion part
            if int.tableau.qdiff2.name ≠ :NULL
                int.Δy .= 0
                for j = 1:i-1
                    for l = 1:equation(int).m
                        int.Δy[l] += int.tableau.qdiff2.a[i,j] * int.B[j][k,l]
                    end
                end

                int.Q[i][k] += dot(int.Δy,sol.ΔZ)/int.Δt
            end

        end
        tᵢ = sol.t̅ + int.Δt * int.tableau.qdrift.c[i]
        int.equation.v(tᵢ, int.Q[i], int.V[i])
        int.equation.B(tᵢ, int.Q[i], int.B[i])
    end

    # compute final update
    if int.tableau.qdiff2.name == :NULL
        update_solution!(sol.q, int.V, int.B, int.tableau.qdrift.b, int.tableau.qdiff.b, int.Δt, sol.ΔW, int.Δy)
    else
        update_solution!(sol.q, int.V, int.B, int.tableau.qdrift.b, int.tableau.qdiff.b, int.tableau.qdiff2.b, int.Δt, sol.ΔW, sol.ΔZ, int.Δy)
    end
end
