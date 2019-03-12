
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
    TableauSERK{T}(name, qdrift, qdiff, CoefficientsRK(T, :NULL, 0, zeros(qdrift.a), zeros(qdrift.b), zeros(qdrift.c)) )
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
                    CoefficientsRK(:NULL, 0, zeros(a_drift), zeros(b_drift), zeros(c_drift)))
end




"Stochastic Explicit Runge-Kutta integrator."
struct IntegratorSERK{DT,TT,FT} <: StochasticIntegrator{DT,TT}
    equation::SDE{DT,TT,FT}
    tableau::TableauSERK{TT}
    Δt::TT
    ΔW::Vector{DT}
    ΔZ::Vector{DT}

    q ::Matrix{Array{DT,1}}
    VQ::Array{DT,2}
    BQ::Array{DT,3}
    tQ::Array{DT,1}
    tV::Array{DT,1}
    tB::Array{DT,2}


    function IntegratorSERK{DT,TT,FT}(equation, tableau, Δt) where {DT,TT,FT}
        D = equation.d
        M = equation.m
        NS= equation.ns
        NI= equation.n
        S = tableau.s

        # create solution vectors
        q = create_solution_vector_double_double(DT, D, NS, NI)

        new(equation, tableau, Δt, zeros(DT,M), zeros(DT,M),
            q, zeros(DT,D,S), zeros(DT,D,M,S),
            zeros(DT,D), zeros(DT,D), zeros(DT,D,M))
    end
end

function IntegratorSERK(equation::SDE{DT,TT,FT}, tableau::TableauSERK{TT}, Δt::TT) where {DT,TT,FT}
    IntegratorSERK{DT,TT,FT}(equation, tableau, Δt)
end

function initialize!(int::IntegratorSERK, sol::SolutionSDE, k::Int, m::Int)
    @assert m ≥ 1
    @assert m ≤ sol.ni
    @assert k ≥ 1
    @assert k ≤ sol.ns


    # copy initial conditions from solution
    get_initial_conditions!(sol, int.q[k,m], k, m)
end

"Integrate SDE with explicit Runge-Kutta integrator."
# Calculating the n-th time step of the explicit integrator for the sample path r and the initial condition m
function integrate_step!(int::IntegratorSERK{DT,TT,FT}, sol::SolutionSDE{DT,TT,NQ,NW}, r::Int, m::Int, n::Int) where {DT,TT,FT,NQ,NW}
    local tᵢ::TT
    local ydrift::DT
    local ydiff = zeros(DT,sol.nm)

    # copy the increments of the Brownian Process
    if NW==1
        #1D Brownian motion, 1 sample path
        int.ΔW[1] = sol.W.ΔW[n-1]
        int.ΔZ[1] = sol.W.ΔZ[n-1]
    elseif NW==2
        #Multidimensional Brownian motion, 1 sample path
        int.ΔW .= sol.W.ΔW[n-1]
        int.ΔZ .= sol.W.ΔZ[n-1]
    elseif NW==3
        #1D or Multidimensional Brownian motion, r-th sample path
        int.ΔW .= sol.W.ΔW[n-1,r]
        int.ΔZ .= sol.W.ΔZ[n-1,r]
    end

    # calculates v(t,tQ) and assigns to the matrix tV
    int.equation.v(sol.t[0] + (n-1)*int.Δt, int.q[r,m], int.tV)
    # copies tV into the i-th column of VQ
    simd_copy_yx_first!(int.tV, int.VQ, 1)
    # calculates B(t,tQ) and assigns to the matrix tB
    int.equation.B(sol.t[0] + (n-1)*int.Δt, int.q[r,m], int.tB)
    # copies the matrix tB to BQ[:,:,i]
    simd_copy_yx_first!(int.tB, int.BQ, 1)


    for i in 2:int.tableau.s
        @inbounds for k in eachindex(int.tQ)
            # contribution from the drift part
            ydrift = 0.
            for j = 1:i-1
                ydrift += int.tableau.qdrift.a[i,j] * int.VQ[k,j]
            end

            # ΔW contribution from the diffusion part
            ydiff .= zeros(DT, sol.nm)
            for j = 1:i-1
                ydiff += int.tableau.qdiff.a[i,j] * int.BQ[k,:,j]
            end

            int.tQ[k] = int.q[r,m][k] + int.Δt * ydrift + dot(ydiff,int.ΔW)

            # ΔZ contribution from the diffusion part
            if int.tableau.qdiff2.name ≠ :NULL
                ydiff .= zeros(DT, sol.nm)
                for j = 1:i-1
                    ydiff += int.tableau.qdiff2.a[i,j] * int.BQ[k,:,j]
                end

                int.tQ[k] = int.tQ[k] + dot(ydiff,int.ΔZ)/int.Δt
            end

        end
        tᵢ = sol.t[0] + (n-1)*int.Δt + int.Δt * int.tableau.qdrift.c[i]
        int.equation.v(tᵢ, int.tQ, int.tV)
        simd_copy_yx_first!(int.tV, int.VQ, i)
        int.equation.B(tᵢ, int.tQ, int.tB)
        simd_copy_yx_first!(int.tB, int.BQ, i)
    end

    # compute final update
    if int.tableau.qdiff2.name == :NULL
        update_solution!(int.q[r,m], int.VQ, int.BQ, int.tableau.qdrift.b, int.tableau.qdiff.b, int.Δt, int.ΔW)
    else
        update_solution!(int.q[r,m], int.VQ, int.BQ, int.tableau.qdrift.b, int.tableau.qdiff.b, int.tableau.qdiff2.b, int.Δt, int.ΔW, int.ΔZ)
    end

    # take care of periodic solutions
    cut_periodic_solution!(int.q[r,m], int.equation.periodicity)

    # copy to solution
    copy_solution!(sol, int.q[r,m], n, r, m)
end
