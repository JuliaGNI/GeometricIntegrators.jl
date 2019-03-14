"""
 Holds the tableau of a stochastic fully implicit split partitioned Runge-Kutta method.
 qdrift, pdrift1, pdrift2 hold the RK coefficients for the drift parts,
 and qdiff, pdiff1, pdiff2 hold the RK coefficients for the diffusion part of the SDE.
"""
struct TableauSFISPRK{T} <: AbstractTableauIRK{T}
    name::Symbol
    s::Int
    qdrift ::CoefficientsRK{T}
    qdiff  ::CoefficientsRK{T}
    pdrift1::CoefficientsRK{T}
    pdrift2::CoefficientsRK{T}
    pdiff1 ::CoefficientsRK{T}
    pdiff2 ::CoefficientsRK{T}

    # Order of the tableau is not included, because unlike in the deterministic
    # setting, it depends on the properties of the noise (e.g., the dimension of
    # the Wiener process and the commutativity properties of the diffusion matrix)
    #
    # Orders stored in qdrift and qdiff are understood as the classical orders of these methods.


    function TableauSFISPRK{T}(name, s, qdrift, qdiff, pdrift1, pdrift2, pdiff1, pdiff2) where {T}
        # THE COMMENTED OUT PART WAS FOR TableauFIRK. MAY IMPLEMENT SOMETHING
        # SIMILAR FOR TableauSFIRK LATER.

        # if (q.s > 1 && istrilstrict(q.a)) || (q.s==1 && q.a[1,1] == 0)
        #     warn("Initializing TableauFIRK with explicit tableau ", q.name, ".\n",
        #          "You might want to use TableauERK instead.")
        # elseif q.s > 1 && istril(q.a)
        #     warn("Initializing TableauFIRK with diagonally implicit tableau ", q.name, ".\n",
        #          "You might want to use TableauDIRK instead.")
        # end

        @assert s == qdrift.s == qdiff.s == pdrift1.s == pdrift2.s == pdiff1.s == pdiff2.s

        new(name, s, qdrift, qdiff, pdrift1, pdrift2, pdiff1, pdiff2)
    end
end

function TableauSFISPRK(name::Symbol, qdrift::CoefficientsRK{T},
                                      qdiff::CoefficientsRK{T},
                                      pdrift1::CoefficientsRK{T}, pdrift2::CoefficientsRK{T},
                                      pdiff1::CoefficientsRK{T}, pdiff2::CoefficientsRK{T}) where {T}
    TableauSFISPRK{T}(name, qdrift.s, qdrift, qdiff, pdrift1, pdrift2, pdiff1, pdiff2)
end

function TableauSFISPRK(name::Symbol, qorder_drift::Int, qa_drift::Matrix{T}, qb_drift::Vector{T}, qc_drift::Vector{T},
                                      qorder_diff::Int, qa_diff::Matrix{T}, qb_diff::Vector{T}, qc_diff::Vector{T},
                                      porder1_drift::Int, pa1_drift::Matrix{T}, pb1_drift::Vector{T}, pc1_drift::Vector{T},
                                      porder2_drift::Int, pa2_drift::Matrix{T}, pb2_drift::Vector{T}, pc2_drift::Vector{T},
                                      porder1_diff::Int, pa1_diff::Matrix{T}, pb1_diff::Vector{T}, pc1_diff::Vector{T},
                                      porder2_diff::Int, pa2_diff::Matrix{T}, pb2_diff::Vector{T}, pc2_diff::Vector{T}) where {T}
    TableauSFISPRK{T}(name, length(qc_drift), CoefficientsRK(name, qorder_drift, qa_drift, qb_drift, qc_drift),
                                              CoefficientsRK(name, qorder_diff, qa_diff, qb_diff, qc_diff),
                                              CoefficientsRK(name, porder1_drift, pa1_drift, pb1_drift, pc1_drift),
                                              CoefficientsRK(name, porder2_drift, pa2_drift, pb2_drift, pc2_drift),
                                              CoefficientsRK(name, porder1_diff, pa1_diff, pb1_diff, pc1_diff),
                                              CoefficientsRK(name, porder2_diff, pa2_diff, pb2_diff, pc2_diff))
end

# TODO function readTableauSFIRKFromFile(dir::AbstractString, name::AbstractString)


# "Parameters for right-hand side function of fully implicit Runge-Kutta methods."
#  A - if positive, the upper bound of the Wiener process increments; if A=0.0, no truncation
mutable struct ParametersSFISPRK{DT, TT, ET <: SPSDE{DT,TT}, D, M, S} <: Parameters{DT,TT}
    equ::ET
    tab::TableauSFISPRK{TT}
    Δt::TT
    ΔW::Vector{DT}
    ΔZ::Vector{DT}
    A::DT

    t::TT
    q::Vector{DT}
    p::Vector{DT}
end

function ParametersSFISPRK(equ::ET, tab::TableauSFISPRK{TT}, Δt::TT, ΔW::Vector{DT}, ΔZ::Vector{DT}, A::DT) where {DT, TT, ET <: SPSDE{DT,TT}}
    @assert equ.m == length(ΔW) == length(ΔZ)
    ParametersSFISPRK{DT, TT, ET, equ.d, equ.m, tab.s}(equ, tab, Δt, ΔW, ΔZ, A, 0, zeros(DT, equ.d), zeros(DT, equ.d))
end


struct NonlinearFunctionCacheSFISPRK{DT}
    # Structure for holding the internal stages Q, the values of the drift vector
    # and the diffusion matrix evaluated at the internal stages VQ=v(Q), BQ=B(Q),
    # and the increments Y = Δt*a_drift*v(Q) + a_diff*B(Q)*ΔW
    Q::Matrix{DT}
    P::Matrix{DT}
    VQP ::Matrix{DT}
    FQP1::Matrix{DT}
    FQP2::Matrix{DT}
    BQP ::Array{DT,3}
    GQP1::Array{DT,3}
    GQP2::Array{DT,3}
    Y::Matrix{DT}
    Z::Matrix{DT}

    vqp ::Vector{DT}
    fqp1::Vector{DT}
    fqp2::Vector{DT}
    bqp ::Matrix{DT}
    gqp1::Matrix{DT}
    gqp2::Matrix{DT}
    y::Vector{DT}
    z::Vector{DT}

    function NonlinearFunctionCacheSFISPRK{DT}(d, m, s) where {DT}

        # create internal stage vectors
        Q  = zeros(DT,d,s)
        P  = zeros(DT,d,s)
        VQP  = zeros(DT,d,s)
        FQP1 = zeros(DT,d,s)
        FQP2 = zeros(DT,d,s)
        BQP  = zeros(DT,d,m,s)
        GQP1 = zeros(DT,d,m,s)
        GQP2 = zeros(DT,d,m,s)
        Y  = zeros(DT,d,s)
        Z  = zeros(DT,d,s)

        # create velocity and update vector
        vqp  = zeros(DT,d)
        fqp1 = zeros(DT,d)
        fqp2 = zeros(DT,d)
        bqp  = zeros(DT,d,m)
        gqp1 = zeros(DT,d,m)
        gqp2 = zeros(DT,d,m)
        y  = zeros(DT,d)
        z  = zeros(DT,d)

        new(Q, P, VQP, FQP1, FQP2, BQP, GQP1, GQP2, Y, Z, vqp, fqp1, fqp2, bqp, gqp1, gqp2, y, z)
    end
end

# Unpacks the data stored in x = (Y[1,1], Y[2,1], ... Y[D,1], Y[1,2], ..., Z[1,1], Z[2,1], ... Z[D,1], Z[1,2], ...)
# into the matrix Y, Z, calculates the internal stages Q, P, the values of the RHS
# of the SDE ( vi(Q,P), fi(Q,P), Bi(Q,P) and Gi(Q,P) ), and assigns them to VQPi, FQPi, BQPi and GQPi.
# Unlike for FIRK, here
# Y = Δt a_drift v(Q,P) + a_diff B(Q,P) ΔW,
# Z = Δt ̃a1_drift f1(Q,P) + Δt ̃a2_drift f2(Q,P) + ̃a1_diff G1(Q,P) ΔW + ̃a2_diff G2(Q,P) ΔW.
@generated function compute_stages!(x::Vector{ST}, Q::Matrix{ST}, P::Matrix{ST},
                                                    VQP::Matrix{ST},
                                                    FQP1::Matrix{ST}, FQP2::Matrix{ST},
                                                    BQP::Array{ST,3},
                                                    GQP1::Array{ST,3}, GQP2::Array{ST,3},
                                                    Y::Matrix{ST}, Z::Matrix{ST},
                                                    params::ParametersSFISPRK{DT,TT,ET,D,M,S}) where {ST,DT,TT,ET,D,M,S}

    tQ ::Vector{ST} = zeros(ST,D)
    tP ::Vector{ST} = zeros(ST,D)
    tV ::Vector{ST} = zeros(ST,D)
    tF1::Vector{ST} = zeros(ST,D)
    tF2::Vector{ST} = zeros(ST,D)
    tB ::Matrix{ST} = zeros(ST,D,M)
    tG1::Matrix{ST} = zeros(ST,D,M)
    tG2::Matrix{ST} = zeros(ST,D,M)

    quote
        local tqᵢ::TT       #times for the q internal stages
        local tpᵢ::TT       #times for the p internal stages

        @assert D == size(Q,1) == size(VQP,1) == size(BQP,1) == size(P,1) == size(FQP1,1) == size(FQP2,1) == size(GQP1,1) == size(GQP2,1)
        @assert S == size(Q,2) == size(VQP,2) == size(BQP,3) == size(P,2) == size(FQP1,2) == size(FQP2,2) == size(GQP1,3) == size(GQP2,3)
        @assert M == size(BQP,2) == size(GQP1,2) == size(GQP2,2)

        # copy x to Y, Z and calculate Q, P
        for i in 1:S
            for k in 1:D
                Y[k,i] = x[D*(i-1)+k]
                Q[k,i] = params.q[k] + Y[k,i]
                Z[k,i] = x[D*(S+i-1)+k]
                P[k,i] = params.p[k] + Z[k,i]
            end
        end

        # compute VQi = vi(Q) and BQi=Bi(Q)
        for i in 1:S
            tqᵢ = params.t + params.Δt * params.tab.qdrift.c[i]
            # Not sure what about pdrift2.c[i] --- should there be another time point?
            tpᵢ = params.t + params.Δt * params.tab.pdrift1.c[i]
            # copies the i-th column of Q, P to the vectors tQ, tP
            simd_copy_xy_first!($tQ, Q, i)
            simd_copy_xy_first!($tP, P, i)
            # calculates v(t,tQ) and assigns to tV
            params.equ.v(tqᵢ, $tQ, $tP, $tV)
            params.equ.f1(tpᵢ, $tQ, $tP, $tF1)
            params.equ.f2(tpᵢ, $tQ, $tP, $tF2)
            # copies tV into the i-th column of VQ
            simd_copy_yx_first!($tV, VQP, i)
            simd_copy_yx_first!($tF1, FQP1, i)
            simd_copy_yx_first!($tF2, FQP2, i)
            # calculates B(t,tQ) and assigns to the matrix tB
            params.equ.B(tqᵢ, $tQ, $tP, $tB)
            params.equ.G1(tpᵢ, $tQ, $tP, $tG1)
            params.equ.G2(tpᵢ, $tQ, $tP, $tG2)
            # copies the matrix tB to BQ[:,:,i]
            simd_copy_yx_first!($tB, BQP, i)
            simd_copy_yx_first!($tG1, GQP1, i)
            simd_copy_yx_first!($tG2, GQP2, i)
        end
    end
end

# "Compute stages of fully implicit Runge-Kutta methods."
@generated function function_stages!(x::Vector{ST}, b::Vector{ST}, params::ParametersSFISPRK{DT,TT,ET,D,M,S}) where {ST,DT,TT,ET,D,M,S}

    cache = NonlinearFunctionCacheSFISPRK{ST}(D, M, S)

    quote
        compute_stages!(x, $cache.Q, $cache.P, $cache.VQP, $cache.FQP1, $cache.FQP2, $cache.BQP, $cache.GQP1, $cache.GQP2, $cache.Y, $cache.Z, params)

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
                    y1 += params.tab.qdrift.a[i,j] * $cache.VQP[k,j] * params.Δt + params.tab.qdiff.a[i,j] * dot($cache.BQP[k,:,j], params.ΔW)
                    y2 += params.tab.qdrift.â[i,j] * $cache.VQP[k,j] * params.Δt + params.tab.qdiff.â[i,j] * dot($cache.BQP[k,:,j], params.ΔW)
                    z1 += params.tab.pdrift1.a[i,j] * $cache.FQP1[k,j] * params.Δt + params.tab.pdrift2.a[i,j] * $cache.FQP2[k,j] * params.Δt
                          + params.tab.pdiff1.a[i,j] * dot($cache.GQP1[k,:,j], params.ΔW) + params.tab.pdiff2.a[i,j] * dot($cache.GQP2[k,:,j], params.ΔW)
                    z2 += params.tab.pdrift1.â[i,j] * $cache.FQP1[k,j] * params.Δt + params.tab.pdrift2.â[i,j] * $cache.FQP2[k,j] * params.Δt
                          + params.tab.pdiff1.â[i,j] * dot($cache.GQP1[k,:,j], params.ΔW) + params.tab.pdiff2.â[i,j] * dot($cache.GQP2[k,:,j], params.ΔW)
                end
                b[D*(i-1)+k] = - $cache.Y[k,i] + (y1 + y2)
                b[D*(S+i-1)+k] = - $cache.Z[k,i] + (z1 + z2)
            end
        end
    end
end


"Stochastic fully implicit partitioned Runge-Kutta integrator."
# InitialGuessPSDE not implemented for SFIPRK
struct IntegratorSFISPRK{DT, TT, PT <: ParametersSFISPRK{DT,TT},
                              ST <: NonlinearSolver{DT}, N} <: StochasticIntegrator{DT,TT}
    params::PT
    solver::ST
    #Not implementing InitialGuessSDE
    #iguess::IT
    fcache::NonlinearFunctionCacheSFISPRK{DT}

    q::Matrix{Vector{TwicePrecision{DT}}}
    p::Matrix{Vector{TwicePrecision{DT}}}
end

# K - the integer in the bound A = √(2 K Δt |log Δt|) due to Milstein & Tretyakov; K=0 no truncation
function IntegratorSFISPRK(equation::SPSDE{DT,TT,VT,F1T,F2T,BT,G1T,G2T,N}, tableau::TableauSFISPRK{TT}, Δt::TT; K::Int=0) where {DT,TT,VT,F1T,F2T,BT,G1T,G2T,N}
    D = equation.d
    M = equation.m
    NS= equation.ns
    NI= equation.n
    S = tableau.s

    # create params
    K==0 ? A = 0.0 : A = sqrt( 2*K*Δt*abs(log(Δt)) )
    params = ParametersSFISPRK(equation, tableau, Δt, zeros(DT,M), zeros(DT,M), A)

    # create solver
    solver = create_nonlinear_solver(DT, 2*D*S, params)

    # Not implementing InitialGuessSDE
    # create initial guess
    #iguess = InitialGuessODE(get_config(:ig_interpolation), equation, Δt)

    # create cache for internal stage vectors and update vectors
    fcache = NonlinearFunctionCacheSFISPRK{DT}(D, M, S)

    # create solution vectors
    q = create_solution_vector(DT, D, NS, NI)
    p = create_solution_vector(DT, D, NS, NI)

    # create integrator
    IntegratorSFISPRK{DT, TT, typeof(params), typeof(solver), N}(params, solver, fcache, q, p)
end

equation(integrator::IntegratorSFISPRK) = integrator.params.equ
timestep(integrator::IntegratorSFISPRK) = integrator.params.Δt
tableau(integrator::IntegratorSFISPRK) = integrator.params.tab
dims(integrator::IntegratorSFISPRK) = integrator.params.equ.d
Base.eltype(integrator::IntegratorSFISPRK{DT, TT, PT, ST, N}) where {DT, TT, PT, ST, N} = DT


function initialize!(int::IntegratorSFISPRK{DT,TT}, sol::SolutionPSDE, k::Int, m::Int) where {DT,TT}
    @assert m ≥ 1
    @assert m ≤ sol.ni
    @assert k ≥ 1
    @assert k ≤ sol.ns

    # copy the m-th initial condition for the k-th sample path
    get_initial_conditions!(sol, int.q[k,m], int.p[k,m], k, m)

    # Not implementing InitialGuessSDE
    # # initialise initial guess
    # initialize!(int.iguess, m, sol.t[0], int.q[m])
end

# NOT IMPLEMENTING InitialGuessSDE
# This function computes initial guesses for Y, Z and assigns them to int.solver.x
# The prediction is calculated using an explicit integrator.
function initial_guess!(int::IntegratorSFISPRK{DT,TT}) where {DT,TT}

    # SIMPLE SOLUTION
    # The simplest initial guess for Y, Z is 0
    int.solver.x .= zeros(eltype(int), 2*int.params.tab.s*dims(int))


    # FIX NEEDED !!!
    # APPARENTLY PREDICTION USING AN EXPLICIT METHOD DOESNT WORK WELL HERE
    # PROBABLY BECAUSE THE pdrift2.c[i] POINTS ARE NOT INCLUDED


    # # USING AN EXPLICIT INTEGRATOR TO COMPUTE AN INITIAL GUESS
    # # Below we use the R2 method of Burrage & Burrage to calculate
    # # the internal stages at the times c[1]...c[s].
    # # This approach seems to give very good approximations if the time step
    # # and magnitude of noise are not too large. If the noise intensity is too big,
    # # one may have to perform a few iterations of the explicit method with a smaller
    # # time step, use a higher-order explicit method (e.g. CL or G5), or use
    # # the simple solution above.
    #
    # local tV1  = zeros(DT,int.params.equ.d)
    # local tV2  = zeros(DT,int.params.equ.d)
    # local tF1  = zeros(DT,int.params.equ.d)
    # local tF2  = zeros(DT,int.params.equ.d)
    # local tmpF = zeros(DT,int.params.equ.d)
    # local tB1  = zeros(DT,int.params.equ.d, int.params.equ.m)
    # local tB2  = zeros(DT,int.params.equ.d, int.params.equ.m)
    # local tG1  = zeros(DT,int.params.equ.d, int.params.equ.m)
    # local tG2  = zeros(DT,int.params.equ.d, int.params.equ.m)
    # local tmpG = zeros(DT,int.params.equ.d, int.params.equ.m)
    # local Q    = zeros(DT,int.params.equ.d)
    # local P    = zeros(DT,int.params.equ.d)
    # local t2::TT
    # local Δt_local::TT
    # local ΔW_local = zeros(DT,int.params.equ.m)
    #
    # # When calling this function, int.params should contain the data:
    # # int.params.q - the q solution at the previous time step
    # # int.params.p - the p solution at the previous time step
    # # int.params.t - the time of the previous step
    # # int.params.ΔW- the increment of the Brownian motion for the current step
    #
    # #Evaluating the functions v and B at t,q - same for all stages
    # int.params.equ.v(int.params.t, int.params.q, int.params.p, tV1)
    # int.params.equ.B(int.params.t, int.params.q, int.params.p, tB1)
    # int.params.equ.f1(int.params.t, int.params.q, int.params.p, tF1)
    # int.params.equ.f2(int.params.t, int.params.q, int.params.p, tmpF)
    # #tF1 = tF1 + tmpF
    # int.params.equ.G1(int.params.t, int.params.q, int.params.p, tG1)
    # int.params.equ.G2(int.params.t, int.params.q, int.params.p, tmpG)
    # #tG1 = tG1 + tmpG
    #
    #
    # # Calculating the positions q at the points qdrift.c[i]
    # # if qdrift.c==pdrift.c, then also calculating the momenta
    # for i in 1:int.params.tab.s
    #
    #     # Taking the c[i] from the qdrift tableau.
    #     Δt_local  = int.params.tab.qdrift.c[i]*int.params.Δt
    #     ΔW_local .= int.params.tab.qdrift.c[i]*int.params.ΔW
    #
    #     Q = int.params.q + 2./3. * Δt_local * tV1 + 2./3. * tB1 * ΔW_local
    #     P = int.params.p + 2./3. * Δt_local * tF1 + 2./3. * tG1 * ΔW_local
    #
    #     t2 = int.params.t + 2./3.*Δt_local
    #
    #     int.params.equ.v(t2, Q, P, tV2)
    #     int.params.equ.B(t2, Q, P, tB2)
    #     int.params.equ.f1(t2, Q, P, tF2)
    #     int.params.equ.f2(t2, Q, P, tmpF)
    #     #tF2 = tF2 + tmpF
    #     int.params.equ.G1(t2, Q, P, tG2)
    #     int.params.equ.G2(t2, Q, P, tmpG)
    #     #tG2 = tG2 + tmpG
    #
    #     #Calculating the Y's and assigning them to the array int.solver.x as initial guesses
    #     for j in 1:int.params.equ.d
    #         int.solver.x[(i-1)*int.params.equ.d+j] =  Δt_local*(1./4.*tV1[j] + 3./4.*tV2[j]) + dot( (1./4.*tB1[j,:] + 3./4.*tB2[j,:]), ΔW_local )
    #     end
    #
    #     # if the collocation points are the same for both q and p parts
    #     # not sure what about pdrift2.c
    #     if int.params.tab.qdrift.c==int.params.tab.pdrift1.c
    #         for j in 1:int.params.equ.d
    #             int.solver.x[(int.params.tab.s+i-1)*int.params.equ.d+j] =  Δt_local*(1./4.*tF1[j] + 3./4.*tF2[j]) + dot( (1./4.*tG1[j,:] + 3./4.*tG2[j,:]), ΔW_local )
    #         end
    #     end
    # end
    #
    #
    # # If qdrift.c != pdrift1.c, then calculating the momenta p at the points pdrift1.c[i]
    # # not sure what about pdrift2.c
    # if int.params.tab.qdrift.c != int.params.tab.pdrift1.c
    #
    #     for i in 1:int.params.tab.s
    #
    #         # Taking the c[i] from the pdrift1 tableau.
    #         Δt_local  = int.params.tab.pdrift1.c[i]*int.params.Δt
    #         ΔW_local .= int.params.tab.pdrift1.c[i]*int.params.ΔW
    #
    #         Q = int.params.q + 2./3. * Δt_local * tV1 + 2./3. * tB1 * ΔW_local
    #         P = int.params.p + 2./3. * Δt_local * tF1 + 2./3. * tG1 * ΔW_local
    #
    #         t2 = int.params.t + 2./3.*Δt_local
    #
    #         int.params.equ.v(t2, Q, P, tV2)
    #         int.params.equ.B(t2, Q, P, tB2)
    #         int.params.equ.f1(t2, Q, P, tF2)
    #         int.params.equ.f2(t2, Q, P, tmpF)
    #         #tF2 = tF2 + tmpF
    #         int.params.equ.G1(t2, Q, P, tG2)
    #         int.params.equ.G2(t2, Q, P, tmpG)
    #         #tG2 = tG2 + tmpG
    #
    #         # Calculating the Z's and assigning them to the array int.solver.x as initial guesses
    #         # The guesses for the Y's have already been written to x above
    #         for j in 1:int.params.equ.d
    #             int.solver.x[(int.params.tab.s+i-1)*int.params.equ.d+j] =  Δt_local*(1./4.*tF1[j] + 3./4.*tF2[j]) + dot( (1./4.*tG1[j,:] + 3./4.*tG2[j,:]), ΔW_local )
    #         end
    #     end
    # end

end


"Integrate PSDE with a stochastic fully implicit partitioned Runge-Kutta integrator."
# Integrating the k-th sample path for the m-th initial condition
function integrate_step!(int::IntegratorSFISPRK{DT,TT}, sol::SolutionPSDE{DT,TT,NQ,NW}, k::Int, m::Int, n::Int) where {DT,TT,NQ,NW}

    @assert k ≥ 1
    @assert k ≤ sol.ns

    @assert m ≥ 1
    @assert m ≤ sol.ni

    @assert n ≥ 1
    @assert n ≤ sol.ntime


    # set time for nonlinear solver
    int.params.t  = sol.t[0] + (n-1)*int.params.Δt
    int.params.q .= int.q[k, m]
    int.params.p .= int.p[k, m]


    # copy the increments of the Brownian Process
    if NW==1
        #1D Brownian motion, 1 sample path
        int.params.ΔW[1] = sol.W.ΔW[n-1]
        int.params.ΔZ[1] = sol.W.ΔZ[n-1]
    elseif NW==2
        #Multidimensional Brownian motion, 1 sample path
        int.params.ΔW .= sol.W.ΔW[n-1]
        int.params.ΔZ .= sol.W.ΔZ[n-1]
    elseif NW==3
        #1D or Multidimensional Brownian motion, k-th sample path
        int.params.ΔW .= sol.W.ΔW[n-1,k]
        int.params.ΔZ .= sol.W.ΔZ[n-1,k]
    end


    # truncate the increments ΔW with A
    if int.params.A>0
        for i in 1:length(int.params.ΔW)
            if int.params.ΔW[i]<-int.params.A
                int.params.ΔW[i] = -int.params.A
            elseif int.params.ΔW[i]>int.params.A
                int.params.ΔW[i] = int.params.A
            end
        end
    end


    # compute initial guess and assign to int.solver.x
    initial_guess!(int)

    # call nonlinear solver
    solve!(int.solver)

    # print solver status
    print_solver_status(int.solver.status, int.solver.params, n)

    # check if solution contains NaNs or error bounds are violated
    check_solver_status(int.solver.status, int.solver.params, n)

    # compute the drift vector field and the diffusion matrix at internal stages
    compute_stages!(int.solver.x, int.fcache.Q, int.fcache.P, int.fcache.VQP, int.fcache.FQP1, int.fcache.FQP2, int.fcache.BQP, int.fcache.GQP1, int.fcache.GQP2, int.fcache.Y, int.fcache.Z, int.params)

    # compute final update
    update_solution!(int.q[k,m], int.p[k,m], int.fcache.VQP, int.fcache.FQP1, int.fcache.FQP2, int.fcache.BQP, int.fcache.GQP1, int.fcache.GQP2,
                    int.params.tab.qdrift.b, int.params.tab.qdrift.b̂, int.params.tab.qdiff.b, int.params.tab.qdiff.b̂,
                    int.params.tab.pdrift1.b, int.params.tab.pdrift1.b̂, int.params.tab.pdrift2.b, int.params.tab.pdrift2.b̂,
                    int.params.tab.pdiff1.b, int.params.tab.pdiff1.b̂, int.params.tab.pdiff2.b, int.params.tab.pdiff2.b̂,
                    int.params.Δt, int.params.ΔW)

    # # NOT IMPLEMENTING InitialGuessSDE
    # # # copy solution to initial guess
    # # update!(int.iguess, m, sol.t[0] + n*int.params.Δt, int.q[m])

    # take care of periodic solutions
    cut_periodic_solution!(int.q[k,m], int.params.equ.periodicity)

    # # copy to solution
    copy_solution!(sol, int.q[k,m], int.p[k,m], n, k, m)
end
