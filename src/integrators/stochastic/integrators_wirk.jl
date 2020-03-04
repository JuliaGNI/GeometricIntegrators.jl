"""
 Holds the tableau of a weak implicit Runge-Kutta method.

   According to Wang, Hong, Xu, "Construction of Symplectic Runge-Kutta Methods for Stochastic Hamiltonian Systems",
   Commun. Comput. Phys. 21(1), 2017
"""
struct TableauWIRK{T} <: AbstractTableauIRK{T}
    name::Symbol
    s::Int
    qdrift0::CoefficientsRK{T}
    qdrift1::CoefficientsRK{T}
    qdiff0::CoefficientsRK{T}
    qdiff1::CoefficientsRK{T}
    qdiff3::CoefficientsRK{T}

    # Order of the tableau is not included, because unlike in the deterministic
    # setting, it depends on the properties of the noise (e.g., the dimension of
    # the Wiener process and the commutativity properties of the diffusion matrix)
    #
    # Orders stored in qdrift and qdiff are understood as the classical orders of these methods.
    #
    # qdrift0, qdrift1 correspond to A0, A1 in the paper
    # qdiff0, qdiff1, qdiff3 correspond to B0, B1, B3
    # qdrift0.b = alpha
    # qdiff0.b  = beta
    # qdrift0.c = c0
    # qdrift1.c = c1


    function TableauWIRK{T}(name, s, qdrift0, qdrift1, qdiff0, qdiff1, qdiff3) where {T}
        # THE COMMENTED OUT PART WAS FOR TableauFIRK. MAY IMPLEMENT SOMETHING
        # SIMILAR FOR TableauSFIRK LATER.

        # if (q.s > 1 && istrilstrict(q.a)) || (q.s==1 && q.a[1,1] == 0)
        #     warn("Initializing TableauFIRK with explicit tableau ", q.name, ".\n",
        #          "You might want to use TableauERK instead.")
        # elseif q.s > 1 && istril(q.a)
        #     warn("Initializing TableauFIRK with diagonally implicit tableau ", q.name, ".\n",
        #          "You might want to use TableauDIRK instead.")
        # end

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
    TableauWIRK(name, CoefficientsRK(name, 0, A0, α, c0), CoefficientsRK(name, 0, A1, α, c1),
                         CoefficientsRK(name, 0, B0, β1, c0), CoefficientsRK(name, 0, B1, β1, c1), CoefficientsRK(name, 0, B3, β1, c1))
end



"Parameters for right-hand side function of weak implicit Runge-Kutta methods."
mutable struct ParametersWIRK{DT, TT, ET <: SDE{DT,TT}, D, M, S} <: Parameters{DT,TT}
    equ::ET
    tab::TableauWIRK{TT}
    Δt::TT
    ΔW::Vector{DT}
    ΔZ::Vector{DT}

    t::TT
    q::Vector{DT}
end

function ParametersWIRK(equ::ET, tab::TableauWIRK{TT}, Δt::TT, ΔW::Vector{DT}, ΔZ::Vector{DT}) where {DT, TT, ET <: SDE{DT,TT}}
    @assert equ.m == length(ΔW) == length(ΔZ)
    ParametersWIRK{DT, TT, ET, equ.d, equ.m, tab.s}(equ, tab, Δt, ΔW, ΔZ, 0, zeros(DT, equ.d))
end


struct IntegratorCacheWIRK{DT}
    # Structure for holding the internal stages Q0, and Q1 the values of the drift vector
    # and the diffusion matrix evaluated at the internal stages V=v(Q0), B=B(Q1),
    # and the increments Y = Δt*a_drift*v(Q) + a_diff*B(Q)*ΔW
    Q0::Vector{Vector{DT}}     # Q0[i][k]   - the k-th component of the internal stage Q^(0)_i
    Q1::Vector{Matrix{DT}}     # Q1[i][k,l] - the k-th component of the internal stage Q^(l)_i
    V ::Vector{Vector{DT}}     # V [i][k]   - the k-th component of the drift vector v(Q0[i][:])
    B ::Vector{Matrix{DT}}     # B [i][:,:] - the diffusion matrix at the stage i such that the l-th column B[i][:,l] is evaluated at Q^(l)_i
    Y0::Vector{Vector{DT}}     # Y0[i][:]   - the increment of the internal stage Q0[i][:]
    Y1::Vector{Matrix{DT}}     # Y1[i][:,l] - the increment of the internal stage Q1[i][:,l]

    v::Vector{DT}
    b::Matrix{DT}
    y::Vector{DT}

    Δy::Vector{DT}

    function IntegratorCacheWIRK{DT}(D, M, S) where {DT}
        # create internal stage vectors
        Q0 = create_internal_stage_vector(DT, D, S)
        Q1 = create_internal_stage_vector(DT, D, M, S)
        V  = create_internal_stage_vector(DT, D, S)
        B  = create_internal_stage_vector(DT, D, M, S)
        Y0 = create_internal_stage_vector(DT, D, S)
        Y1 = create_internal_stage_vector(DT, D, M, S)

        # create velocity and update vector
        v = zeros(DT,D)
        b = zeros(DT,D,M)
        y = zeros(DT,D)

        # create temporary vectors
        Δy  = zeros(DT,M)

        new(Q0, Q1, V, B, Y0, Y1, v, b, y, Δy)
    end
end

"Stochastic implicit Runge-Kutta integrator."
struct IntegratorWIRK{DT, TT, PT <: ParametersWIRK{DT,TT},
                              ST <: NonlinearSolver{DT}} <: StochasticIntegrator{DT,TT}
    params::PT
    solver::ST
    cache::IntegratorCacheWIRK{DT}
end


function IntegratorWIRK(equation::SDE{DT,TT}, tableau::TableauWIRK{TT}, Δt::TT) where {DT,TT}
    D = equation.d
    M = equation.m
    NS= max(equation.ns,equation.ni)
    S = tableau.s

    # create params
    params = ParametersWIRK(equation, tableau, Δt, zeros(DT,M), zeros(DT,M))

    # create solver
    solver = create_nonlinear_solver(DT, D*S*(M+1), params)

    # create cache for internal stage vectors and update vectors
    cache = IntegratorCacheWIRK{DT}(D, M, S)

    # create integrator
    IntegratorWIRK{DT, TT, typeof(params), typeof(solver)}(params, solver, cache)
end

@inline equation(integrator::IntegratorWIRK) = integrator.params.equ
@inline timestep(integrator::IntegratorWIRK) = integrator.params.Δt
@inline tableau(integrator::IntegratorWIRK) = integrator.params.tab
@inline nstages(integrator::IntegratorWIRK)  = nstages(tableau(integrator))
@inline eachstage(integrator::IntegratorWIRK) = 1:nstages(integrator)
@inline noisedims(integrator::IntegratorWIRK) = integrator.params.equ.m
@inline Base.eltype(integrator::IntegratorWIRK{DT}) where {DT} = DT


function update_params!(int::IntegratorWIRK, sol::AtomicSolutionSDE)
    # set time for nonlinear solver and copy previous solution
    int.params.t  = sol.t
    int.params.q .= sol.q
    int.params.ΔW .= sol.ΔW
    int.params.ΔZ .= sol.ΔZ
end


"""
Unpacks the data stored in
x = (Y0[1][1], Y0[1][2]], ... Y0[1][D], ... Y0[S][D], Y1[1][1,1], Y1[1][2,1], ... Y1[1][D,1], Y1[1][1,2], Y1[1][2,2], ... Y1[1][D,2], ... Y1[S][D,M]  )
into Y0::Vector{Vector} and Y1::Vector{Matrix}, calculates the internal stages Q0 and Q1, the values of the RHS
of the SDE ( v(Q0) and B(Q1) ), and assigns them to V and B.
Unlike for FIRK, here Y = Δt a v(Q) + ̃a B(Q) ΔW
"""
@generated function compute_stages!(x::Vector{ST}, Q0::Vector{Vector{ST}}, Q1::Vector{Matrix{ST}}, V::Vector{Vector{ST}},
                                                   B::Vector{Matrix{ST}}, Y0::Vector{Vector{ST}}, Y1::Vector{Matrix{ST}},
                                                   params::ParametersWIRK{DT,TT,ET,D,M,S}) where {ST,DT,TT,ET,D,M,S}

    tQ::Vector{ST} = zeros(ST,D)
    tV::Vector{ST} = zeros(ST,D)
    tB::Vector{ST} = zeros(ST,D)    # This was tB::Matrix{ST} = zeros(ST,D,M) for SIRK, but here we calculate B column by column.
                                    # tB is actually superfluous, could just use tV; keeping tB for clarity though

    quote
        local tᵢ::TT

        @assert S == length(Q0) == length(Q1) == length(V) == length(B)

        # copy x to Y0 and Y1, and calculate Q0 and Q1
        for i in eachindex(Q0)

            @assert D == length(Q0[i]) == size(Q1[i],1) == length(V[i]) == size(B[i],1)
            @assert M == size(Q1[i],2) == size(B[i],2)

            for k in eachindex(Q0[i])
                Y0[i][k] = x[D*(i-1)+k]
                Q0[i][k] = params.q[k] + Y0[i][k]
            end
        end

        for i in eachindex(Q1)
            for l in 1:size(Q1[i],2)
                for k in 1:size(Q1[i],1)
                    Y1[i][k,l] = x[ D*S + D*M*(i-1) + D*(l-1) + k ]
                    Q1[i][k,l] = params.q[k] + Y1[i][k,l]
                end
            end
        end


        # compute V = v(Q0) and B=B(Q1)
        for i in 1:S
            # time point with the coefficient c0
            tᵢ = params.t + params.Δt * params.tab.qdrift0.c[i]
            # calculates v(t,tQ) and assigns to the i-th column of V
            params.equ.v(tᵢ, Q0[i], V[i])

            # time point with the coefficient c1
            tᵢ = params.t + params.Δt * params.tab.qdrift1.c[i]

            for l=1:M
                # copies Q1[i][:,l] to the vector tQ
                simd_copy_xy_first!($tQ, Q1[i], l)
                # calculates the l-th column of B(t,tQ) and assigns to the vector tB
                params.equ.B(tᵢ, $tQ, $tB)
                # copies the vector tB to B[i][:,l]
                simd_copy_yx_first!($tB, B[i], l)
            end
        end
    end
end

"Compute stages of weak implicit Runge-Kutta methods."
@generated function Integrators.function_stages!(x::Vector{ST}, b::Vector{ST}, params::ParametersWIRK{DT,TT,ET,D,M,S}) where {ST,DT,TT,ET,D,M,S}

    cache = IntegratorCacheWIRK{ST}(D, M, S)

    quote
        compute_stages!(x, $cache.Q0, $cache.Q1, $cache.V, $cache.B, $cache.Y0, $cache.Y1, params)

        local y1::ST
        local y2::ST

        # compute b = - (Y-AV)
        # the terms corresponding to Y0
        for i in 1:S
            for k in 1:D
                y1 = 0
                y2 = 0
                for j in 1:S
                    y1 += params.tab.qdrift0.a[i,j] * $cache.V[j][k] * params.Δt
                    y2 += params.tab.qdrift0.â[i,j] * $cache.V[j][k] * params.Δt
                    for l in 1:M
                        y1 += params.tab.qdiff0.a[i,j] * $cache.B[j][k,l] * params.ΔW[l]
                        y2 += params.tab.qdiff0.â[i,j] * $cache.B[j][k,l] * params.ΔW[l]
                    end
                end
                b[D*(i-1)+k] = - $cache.Y0[i][k] + (y1 + y2)
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
                        y1 += params.tab.qdrift1.a[i,j] * $cache.V[j][k] * params.Δt
                        y2 += params.tab.qdrift1.â[i,j] * $cache.V[j][k] * params.Δt

                        # The noise terms, calculated either with the B1 or B3 tableau
                        for noise_idx in 1:M
                            if noise_idx==l
                                y1 += params.tab.qdiff1.a[i,j] * $cache.B[j][k,noise_idx] * params.ΔW[noise_idx]
                                y2 += params.tab.qdiff1.â[i,j] * $cache.B[j][k,noise_idx] * params.ΔW[noise_idx]
                            else
                                y1 += params.tab.qdiff3.a[i,j] * $cache.BQ[j][k,noise_idx] * params.ΔW[noise_idx]
                                y2 += params.tab.qdiff3.â[i,j] * $cache.BQ[j][k,noise_idx] * params.ΔW[noise_idx]
                            end
                        end
                    end
                    b[D*S + D*M*(i-1) + D*(l-1) + k] = - $cache.Y1[i][k,l] + (y1 + y2)
                end
            end
        end
    end
end


"""
This function computes initial guesses for Y and assigns them to int.solver.x
For WIRK we are NOT IMPLEMENTING an InitialGuess.
"""
function initial_guess!(int::IntegratorWIRK{DT,TT}, sol::AtomicSolutionSDE{DT,TT}) where {DT,TT}
    # The simplest initial guess for Y is 0
    int.solver.x .= 0

    # Using an explicit integrator to predict the next step's value (like in SIRK)
    # does not seem to be a good idea here, because the integrators are convergent
    # in the weak sense only, and there is no guarantee that the explicit integrator
    # will produce anything close to the desired solution...
end


"""
Integrate SDE with a stochastic implicit Runge-Kutta integrator.
  Integrating the m-th sample path
"""
function Integrators.integrate_step!(int::IntegratorWIRK{DT,TT}, sol::AtomicSolutionSDE{DT,TT}) where {DT,TT}
    # update nonlinear solver parameters from cache
    update_params!(int, sol)

    # compute initial guess and assign to int.solver.x
    initial_guess!(int, sol)

    # reset cache
    reset!(sol, timestep(int))

    # call nonlinear solver
    solve!(int.solver)

    # print solver status
    print_solver_status(int.solver.status, int.solver.params)

    # check if solution contains NaNs or error bounds are violated
    check_solver_status(int.solver.status, int.solver.params)

    # compute the drift vector field and the diffusion matrix at internal stages
    compute_stages!(int.solver.x, int.cache.Q0, int.cache.Q1, int.cache.V, int.cache.B, int.cache.Y0, int.cache.Y1, int.params)

    # compute final update (same update function as for SIRK)
    update_solution!(sol.q, int.cache.V, int.cache.B, int.params.tab.qdrift0.b, int.params.tab.qdrift0.b̂, int.params.tab.qdiff0.b, int.params.tab.qdiff0.b̂, int.params.Δt, int.params.ΔW, int.cache.Δy)
end
