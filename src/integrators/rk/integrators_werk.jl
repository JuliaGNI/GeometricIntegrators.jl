"""
 Holds the tableau of a weak explicit Runge-Kutta method.

   According to Andreas Rossler, "Second order Runge-Kutta methods for Stratonovich stochastic differential equations",
   BIT Numerical Mathematics (2007) 47, equation (5.1)
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


    # Order of the tableau is not included, because unlike in the deterministic
    # setting, it depends on the properties of the noise (e.g., the dimension of
    # the Wiener process and the commutativity properties of the diffusion matrix)
    #
    # Orders stored in qdrift, qdiff are irrelevant and set to 0
    #
    # qdrift0, qdrift1, qdrift2 correspond to A0, A1, A2 in the paper
    # qdiff0, qdiff1, qdiff2, qdiff3 correspond to B0, B1, B2, B3
    # qdrift0.b = alpha
    # qdiff0.b  = beta1
    # qdiff3.b  = beta2
    # qdrift0.c = c0
    # qdrift1.c = c1
    # qdrift2.c = c2

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
struct IntegratorWERK{DT,TT,FT} <: StochasticIntegrator{DT,TT}
    equation::SDE{DT,TT,FT}
    tableau::TableauWERK{TT}
    Δt::TT
    ΔW::Vector{DT}
    ΔZ::Vector{DT}

    q ::Matrix{Vector{DT}}
    Q0::Vector{DT}           # Q0[1..n]      - the internal stage H^(0)_i for a given i
    Q1::Vector{Vector{DT}}   # Q1[1..n,1..m] - the internal stages H^(l)_i for a given i
    Q2::Vector{Vector{DT}}   # Q2[1..n,1..m] - the internal stages \hat H^(l)_i for a given i
    V ::Vector{Vector{DT}}
    B1::Vector{Matrix{DT}}   # B1[i] holds the values of the diffusion matrix such that BQ1[i][:,l] is evaluated at H^(l)_i
    B2::Vector{Matrix{DT}}   # B2[i] holds the values of the diffusion matrix such that BQ2[i][:,l] is evaluated at \hat H^(l)_i
    # tV::Array{DT,1}          # the value of the drift term evaluated at the internal stage H^(0)_i or the l-th column of the diffusion term for H^(l)_i for a given i
    # tB::Array{DT,2}          # the value of the diffusion matrix; here used only at the first stage of the step
    tB::Vector{DT}


    function IntegratorWERK{DT,TT,FT}(equation, tableau, Δt) where {DT,TT,FT}
        D = equation.d
        M = equation.m
        NS= equation.ns
        NI= equation.n
        S = tableau.s

        # create solution vectors
        q = create_solution_vector(DT, D, NS, NI)

        # create internal stage vectors
        Q0 = zeros(DT, D)
        Q1 = create_internal_stage_vector(DT, D, M)
        Q2 = create_internal_stage_vector(DT, D, M)
        V  = create_internal_stage_vector(DT, D, S)
        B1 = create_internal_stage_vector(DT, D, M, S)
        B2 = create_internal_stage_vector(DT, D, M, S)

        new(equation, tableau, Δt, zeros(DT,M), zeros(DT,M), q, Q0, Q1, Q2, V, B1, B2, zeros(DT,D))
    end
end

function IntegratorWERK(equation::SDE{DT,TT,FT}, tableau::TableauWERK{TT}, Δt::TT) where {DT,TT,FT}
    IntegratorWERK{DT,TT,FT}(equation, tableau, Δt)
end

function initialize!(int::IntegratorWERK, sol::SolutionSDE, k::Int, m::Int)
    @assert m ≥ 1
    @assert m ≤ sol.ni
    @assert k ≥ 1
    @assert k ≤ sol.ns

    # copy initial conditions from solution
    get_initial_conditions!(sol, int.q[k,m], k, m)
end

"""
Integrate SDE with explicit Runge-Kutta integrator.
 Calculating the n-th time step of the explicit integrator for the sample path r and the initial condition m
"""
function integrate_step!(int::IntegratorWERK{DT,TT,FT}, sol::SolutionSDE{DT,TT,NQ,NW}, r::Int, m::Int, n::Int) where {DT,TT,FT,NQ,NW}
    local tᵢ::TT
    local ydrift::DT
    local ydiff  = zeros(DT,sol.nm)
    local ydiff2::DT

    # copy the random variables \hat I and \tilde I representing the Wiener process
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

    # THE FIRST INTERNAL STAGES ARE ALL EQUAL TO THE PREVIOUS STEP SOLUTION
    # calculates v(t,tQ0) and assigns to the 1st column of V
    int.equation.v(sol.t[0] + (n-1)*int.Δt, int.q[r,m], int.V[1])
    # calculates B(t,Q) and assigns to the matrix tB
    # no need to calculate the columns of B separately, because all internal stages
    # for i=1 are equal to int.q[r,m]
    int.equation.B(sol.t[0] + (n-1)*int.Δt, int.q[r,m], int.B1[1])
    # copy B1[i] to B2[i] (they're equal for i=1)
    int.B2[1] .= int.B1[1]


    for i in 2:int.tableau.s

        # Calculating the internal stage H^(0)_i
        @inbounds for k in eachindex(int.Q0)
            # contribution from the drift part
            ydrift = 0.
            for j = 1:i-1
                ydrift += int.tableau.qdrift0.a[i,j] * int.V[j][k]
            end

            # ΔW contribution from the diffusion part
            ydiff .= zeros(DT, sol.nm)
            for j = 1:i-1
                ydiff += int.tableau.qdiff0.a[i,j] * int.B1[j][k,:]
            end

            int.Q0[k] = int.q[r,m][k] + int.Δt * ydrift + dot(ydiff,int.ΔW)

        end


        # Calculating the internal stages H^(l)_i for l=1..sol.nm
        @inbounds for k in eachindex(int.Q1[1][:])
            # contribution from the drift part (same for all noises)
            ydrift = 0.
            for j = 1:i-1
                ydrift += int.tableau.qdrift1.a[i,j] * int.V[j][k]
            end

            # contribution from the diffusion part is different for each noise
            @inbounds for noise_idx in eachindex(int.Q1)

                # ΔW contribution from the diffusion part
                ydiff .= zeros(DT, sol.nm)
                for j = 1:i-1
                    for l = 1:sol.nm
                        if l==noise_idx
                            ydiff[l] += int.tableau.qdiff1.a[i,j] * int.B1[j][k,l]
                        else
                            ydiff[l] += int.tableau.qdiff3.a[i,j] * int.B1[j][k,l]
                        end
                    end
                end

                int.Q1[noise_idx][k] = int.q[r,m][k] + int.Δt * ydrift + dot(ydiff,int.ΔW)
            end
        end


        # Calculating the internal stages \hat H^(l)_i for l=1..sol.nm
        @inbounds for k in eachindex(int.Q2)
            # contribution from the drift part (same for all noises)
            ydrift = 0.
            for j = 1:i-1
                ydrift += int.tableau.qdrift2.a[i,j] * int.V[j][k]
            end

            # contribution from the diffusion part is different for each noise
            @inbounds for noise_idx in eachindex(int.Q2)

                # ΔW contribution from the diffusion part
                ydiff2 = 0.0

                for j = 1:i-1
                    for l = 1:sol.nm
                        # calculating the terms with the random variables \hat I_{noise_idx,l}
                        if l<noise_idx
                            ydiff2 += int.tableau.qdiff2.a[i,j] * int.B1[j][k,l] * int.ΔW[noise_idx] * int.ΔZ[l]
                        elseif l>noise_idx
                            ydiff2 += -int.tableau.qdiff2.a[i,j] * int.B1[j][k,l] * int.ΔW[l] * int.ΔZ[noise_idx]
                        end
                    end
                end

                int.Q2[noise_idx][k] = int.q[r,m][k] + int.Δt * ydrift + ydiff2/sqrt(int.Δt)
            end
        end

        # CALCULATING THE NEW VALUES OF V
        tᵢ = sol.t[0] + (n-1)*int.Δt + int.Δt * int.tableau.qdrift0.c[i]
        int.equation.v(tᵢ, int.Q0, int.V[i])

        # CALCULATING THE NEW VALUES OF BQ1
        # each column of B evaluated at a different internal stage
        tᵢ = sol.t[0] + (n-1)*int.Δt + int.Δt * int.tableau.qdrift1.c[i]

        for l = 1:sol.nm
            # here int.tV holds the l-th column of B
            # TODO check for correctness
            int.equation.B(tᵢ, int.Q1[l], int.tB, col=l)
            simd_copy_yx_first!(int.tB, int.B1[i], l)
        end


        #CALCULATING THE NEW VALUES OF BQ2
        # each column of B evaluated at a different internal stage
        tᵢ = sol.t[0] + (n-1)*int.Δt + int.Δt * int.tableau.qdrift2.c[i]

        for l = 1:sol.nm
            # here int.tV holds the l-th column of B
            # TODO check for correctness
            int.equation.B(tᵢ, int.Q2[l], int.tB, col=l)
            simd_copy_yx_first!(int.tB, int.B2[i], l)
        end

    end

    # compute final update
    update_solution!(int.q[r,m], int.V, int.B1, int.B2, int.tableau.qdrift0.b, int.tableau.qdiff0.b, int.tableau.qdiff3.b, int.Δt, int.ΔW)

    # take care of periodic solutions
    cut_periodic_solution!(int.q[r,m], int.equation.periodicity)

    # copy to solution
    copy_solution!(sol, int.q[r,m], n, r, m)
end
