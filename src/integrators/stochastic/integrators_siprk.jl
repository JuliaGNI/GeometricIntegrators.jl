"""
 Holds the tableau of a stochastic implicit partitioned Runge-Kutta method.
 qdrift, pdrift hold the RK coefficients for the drift part,
 and qdiff, pdiff hold the RK coefficients for the diffusion part of the SDE.
"""
struct TableauSIPRK{T} <: AbstractTableauIRK{T}
    name::Symbol
    s::Int
    qdrift::CoefficientsRK{T}
    qdiff::CoefficientsRK{T}
    pdrift::CoefficientsRK{T}
    pdiff::CoefficientsRK{T}

    # Order of the tableau is not included, because unlike in the deterministic
    # setting, it depends on the properties of the noise (e.g., the dimension of
    # the Wiener process and the commutativity properties of the diffusion matrix)
    #
    # Orders stored in qdrift and qdiff are understood as the classical orders of these methods.


    function TableauSIPRK{T}(name, s, qdrift, qdiff, pdrift, pdiff) where {T}
        # THE COMMENTED OUT PART WAS FOR TableauFIRK. MAY IMPLEMENT SOMETHING
        # SIMILAR FOR TableauSFIRK LATER.

        # if (q.s > 1 && istrilstrict(q.a)) || (q.s==1 && q.a[1,1] == 0)
        #     warn("Initializing TableauFIRK with explicit tableau ", q.name, ".\n",
        #          "You might want to use TableauERK instead.")
        # elseif q.s > 1 && istril(q.a)
        #     warn("Initializing TableauFIRK with diagonally implicit tableau ", q.name, ".\n",
        #          "You might want to use TableauDIRK instead.")
        # end

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


"""
Parameters for right-hand side function of implicit Runge-Kutta methods.
  A - if positive, the upper bound of the Wiener process increments; if A=0.0, no truncation
"""
mutable struct ParametersSIPRK{DT, TT, ET <: PSDE{DT,TT}, D, M, S} <: Parameters{DT,TT}
    equ::ET
    tab::TableauSIPRK{TT}
    Δt::TT
    ΔW::Vector{DT}
    ΔZ::Vector{DT}
    A::DT

    t::TT
    q::Vector{DT}
    p::Vector{DT}
end

function ParametersSIPRK(equ::ET, tab::TableauSIPRK{TT}, Δt::TT, ΔW::Vector{DT}, ΔZ::Vector{DT}, A::DT) where {DT, TT, ET <: PSDE{DT,TT}}
    @assert equ.m == length(ΔW) == length(ΔZ)
    ParametersSIPRK{DT, TT, ET, equ.d, equ.m, tab.s}(equ, tab, Δt, ΔW, ΔZ, A, 0, zeros(DT, equ.d), zeros(DT, equ.d))
end


"""
Structure for holding the internal stages Q, the values of the drift vector
and the diffusion matrix evaluated at the internal stages VQ=v(Q), BQ=B(Q),
and the increments Y = Δt*a_drift*v(Q) + a_diff*B(Q)*ΔW
"""
struct NonlinearFunctionCacheSIPRK{DT}
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

    function NonlinearFunctionCacheSIPRK{DT}(d, m, s) where {DT}

        # create internal stage vectors
        Q = create_internal_stage_vector(DT, d, s)
        P = create_internal_stage_vector(DT, d, s)
        V = create_internal_stage_vector(DT, d, s)
        F = create_internal_stage_vector(DT, d, s)
        B = create_internal_stage_vector(DT, d, m, s)
        G = create_internal_stage_vector(DT, d, m, s)
        Y = create_internal_stage_vector(DT, d, s)
        Z = create_internal_stage_vector(DT, d, s)

        # create velocity and update vector
        v = zeros(DT,d)
        f = zeros(DT,d)
        b = zeros(DT,d,m)
        g = zeros(DT,d,m)
        y = zeros(DT,d)
        z = zeros(DT,d)

        new(Q, P, V, F, B, G, Y, Z, v, f, b, g, y, z)
    end
end

"""
Unpacks the data stored in x = (Y[1][1], Y[1][2], ... Y[1][D], Y[2][1], ..., Z[1][1], Z[1][2], ... Z[1][D], Z[2][1], ...)
into Y, Z::Vector{Vector}, calculates the internal stages Q, P, the values of the RHS
of the SDE ( v(Q,P), f(Q,P), B(Q,P) and G(Q,P) ), and assigns them to V, F, B and G.
Unlike for FIRK, here
Y = Δt a_drift v(Q,P) + a_diff B(Q,P) ΔW,
Z = Δt ̃a_drift v(Q,P) + ̃a_diff B(Q,P) ΔW.
"""
function compute_stages!(x::Vector{ST}, Q::Vector{Vector{ST}}, P::Vector{Vector{ST}},
                                        V::Vector{Vector{ST}}, F::Vector{Vector{ST}},
                                        B::Vector{Matrix{ST}}, G::Vector{Matrix{ST}},
                                        Y::Vector{Vector{ST}}, Z::Vector{Vector{ST}},
                                        params::ParametersSIPRK{DT,TT,ET,D,M,S}) where {ST,DT,TT,ET,D,M,S}

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
        params.equ.v(tqᵢ, Q[i], P[i], V[i])
        params.equ.f(tpᵢ, Q[i], P[i], F[i])
        # calculates B(t,Q,P) and assigns to the matrices B[i] and G[i]
        params.equ.B(tqᵢ, Q[i], P[i], B[i])
        params.equ.G(tpᵢ, Q[i], P[i], G[i])
    end
end

"Compute stages of implicit Runge-Kutta methods."
@generated function function_stages!(x::Vector{ST}, b::Vector{ST}, params::ParametersSIPRK{DT,TT,ET,D,M,S}) where {ST,DT,TT,ET,D,M,S}

    cache = NonlinearFunctionCacheSIPRK{ST}(D, M, S)

    quote
        compute_stages!(x, $cache.Q, $cache.P, $cache.V, $cache.F, $cache.B, $cache.G, $cache.Y, $cache.Z, params)

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
                    y1 += params.tab.qdrift.a[i,j] * $cache.V[j][k] * params.Δt
                    y2 += params.tab.qdrift.â[i,j] * $cache.V[j][k] * params.Δt
                    z1 += params.tab.pdrift.a[i,j] * $cache.F[j][k] * params.Δt
                    z2 += params.tab.pdrift.â[i,j] * $cache.F[j][k] * params.Δt
                    for l in 1:M
                        y1 += params.tab.qdiff.a[i,j] * $cache.B[j][k,l] * params.ΔW[l]
                        y2 += params.tab.qdiff.â[i,j] * $cache.B[j][k,l] * params.ΔW[l]
                        z1 += params.tab.pdiff.a[i,j] * $cache.G[j][k,l] * params.ΔW[l]
                        z2 += params.tab.pdiff.â[i,j] * $cache.G[j][k,l] * params.ΔW[l]
                    end
                end
                b[D*(  i-1)+k] = - $cache.Y[i][k] + (y1 + y2)
                b[D*(S+i-1)+k] = - $cache.Z[i][k] + (z1 + z2)
            end
        end
    end
end


"Stochastic implicit partitioned Runge-Kutta integrator."
struct IntegratorSIPRK{DT, TT, PT <: ParametersSIPRK{DT,TT},
                               ST <: NonlinearSolver{DT}, N} <: StochasticIntegrator{DT,TT}
    params::PT
    solver::ST
    # InitialGuessPSDE not implemented for SIPRK
    #iguess::IT
    fcache::NonlinearFunctionCacheSIPRK{DT}

    Q::Vector{DT}
    P::Vector{DT}
    tV1::Vector{DT}
    tV2::Vector{DT}
    tF1::Vector{DT}
    tF2::Vector{DT}
    tB1::Matrix{DT}
    tB2::Matrix{DT}
    tG1::Matrix{DT}
    tG2::Matrix{DT}
    tΔW::Vector{DT}
    ΔQ::Vector{DT}
    ΔQ1::Vector{DT}
    ΔQ2::Vector{DT}
    ΔP::Vector{DT}
    ΔP1::Vector{DT}
    ΔP2::Vector{DT}
    Δy::Vector{DT}
    Δz::Vector{DT}

    q::Matrix{Vector{TwicePrecision{DT}}}
    p::Matrix{Vector{TwicePrecision{DT}}}
end

# K - the integer in the bound A = √(2 K Δt |log Δt|) due to Milstein & Tretyakov; K=0 no truncation
function IntegratorSIPRK(equation::PSDE{DT,TT,VT,BT,N}, tableau::TableauSIPRK{TT}, Δt::TT; K::Int=0) where {DT,TT,VT,BT,N}
    D = equation.d
    M = equation.m
    NS= equation.ns
    NI= equation.ni
    S = tableau.s

    # create params
    K==0 ? A = 0.0 : A = sqrt( 2*K*Δt*abs(log(Δt)) )
    params = ParametersSIPRK(equation, tableau, Δt, zeros(DT,M), zeros(DT,M), A)

    # create solver
    solver = create_nonlinear_solver(DT, 2*D*S, params)

    # Not implementing InitialGuessSDE
    # create initial guess
    #iguess = InitialGuessODE(get_config(:ig_interpolation), equation, Δt)

    # create cache for internal stage vectors and update vectors
    fcache = NonlinearFunctionCacheSIPRK{DT}(D, M, S)

    # create temporary arrays
    Q   = zeros(DT,D)
    P   = zeros(DT,D)
    tV1 = zeros(DT,D)
    tV2 = zeros(DT,D)
    tF1 = zeros(DT,D)
    tF2 = zeros(DT,D)
    tB1 = zeros(DT,D,M)
    tB2 = zeros(DT,D,M)
    tG1 = zeros(DT,D,M)
    tG2 = zeros(DT,D,M)
    tΔW = zeros(DT,M)
    ΔQ  = zeros(DT,D)
    ΔQ1 = zeros(DT,D)
    ΔQ2 = zeros(DT,D)
    ΔP  = zeros(DT,D)
    ΔP1 = zeros(DT,D)
    ΔP2 = zeros(DT,D)
    Δy  = zeros(DT,M)
    Δz  = zeros(DT,M)

    # create solution vectors
    q = create_solution_vector(DT, D, NS, NI)
    p = create_solution_vector(DT, D, NS, NI)

    # create integrator
    IntegratorSIPRK{DT, TT, typeof(params), typeof(solver), N}(params, solver,
                fcache, Q, P, tV1, tV2, tF1, tF2, tB1, tB2, tG1, tG2, tΔW, ΔQ, ΔQ1, ΔQ2, ΔP, ΔP1, ΔP2, Δy, Δz, q, p)
end

equation(integrator::IntegratorSIPRK) = integrator.params.equ
timestep(integrator::IntegratorSIPRK) = integrator.params.Δt
tableau(integrator::IntegratorSIPRK) = integrator.params.tab
Base.eltype(integrator::IntegratorSIPRK{DT, TT, PT, ST, N}) where {DT, TT, PT, ST, N} = DT


function initialize!(int::IntegratorSIPRK{DT,TT}, sol::SolutionPSDE, k::Int, m::Int) where {DT,TT}
    check_solution_dimension_asserts(sol, k, m)

    # copy the m-th initial condition for the k-th sample path
    get_initial_conditions!(sol, int.q[k,m], int.p[k,m], k, m)

    # Not implementing InitialGuessSDE
    # # initialise initial guess
    # initialize!(int.iguess, m, sol.t[0], int.q[m])
end

"""
This function computes initial guesses for Y, Z and assigns them to int.solver.x
The prediction is calculated using an explicit integrator.
"""
function initial_guess!(int::IntegratorSIPRK{DT,TT}) where {DT,TT}

    # NOT IMPLEMENTING InitialGuessSDE

    # SIMPLE SOLUTION
    # The simplest initial guess for Y, Z is 0
    # int.solver.x .= zeros(eltype(int), 2*int.params.tab.s*ndims(int))

    # USING AN EXPLICIT INTEGRATOR TO COMPUTE AN INITIAL GUESS
    # Below we use the R2 method of Burrage & Burrage to calculate
    # the internal stages at the times c[1]...c[s].
    # This approach seems to give very good approximations if the time step
    # and magnitude of noise are not too large. If the noise intensity is too big,
    # one may have to perform a few iterations of the explicit method with a smaller
    # time step, use a higher-order explicit method (e.g. CL or G5), or use
    # the simple solution above.

    local t2::TT
    local Δt_local::TT

    # When calling this function, int.params should contain the data:
    # int.params.q - the q solution at the previous time step
    # int.params.p - the p solution at the previous time step
    # int.params.t - the time of the previous step
    # int.params.ΔW- the increment of the Brownian motion for the current step

    #Evaluating the functions v and B at t,q - same for all stages
    int.params.equ.v(int.params.t, int.params.q, int.params.p, int.tV1)
    int.params.equ.B(int.params.t, int.params.q, int.params.p, int.tB1)
    int.params.equ.f(int.params.t, int.params.q, int.params.p, int.tF1)
    int.params.equ.G(int.params.t, int.params.q, int.params.p, int.tG1)


    # Calculating the positions q at the points qdrift.c[i]
    # if qdrift.c==pdrift.c, then also calculating the momenta
    for i in 1:int.params.tab.s

        # Taking the c[i] from the qdrift tableau.
        Δt_local = int.params.tab.qdrift.c[i]  * int.params.Δt
        int.tΔW .= int.params.tab.qdrift.c[i] .* int.params.ΔW

        simd_mult!(int.ΔQ, int.tB1, int.tΔW)
        simd_mult!(int.ΔP, int.tG1, int.tΔW)
        @. int.Q = int.params.q + 2. / 3. * Δt_local * int.tV1 + 2. / 3. * int.ΔQ
        @. int.P = int.params.p + 2. / 3. * Δt_local * int.tF1 + 2. / 3. * int.ΔP

        t2 = int.params.t + 2. / 3. * Δt_local

        int.params.equ.v(t2, int.Q, int.P, int.tV2)
        int.params.equ.B(t2, int.Q, int.P, int.tB2)
        int.params.equ.f(t2, int.Q, int.P, int.tF2)
        int.params.equ.G(t2, int.Q, int.P, int.tG2)

        #Calculating the Y's and assigning them to the array int.solver.x as initial guesses
        simd_mult!(int.ΔQ1, int.tB1, int.tΔW)
        simd_mult!(int.ΔQ2, int.tB2, int.tΔW)
        for j in 1:int.params.equ.d
            int.solver.x[(i-1)*int.params.equ.d+j] =  Δt_local*(1. / 4. * int.tV1[j] + 3. / 4. * int.tV2[j]) + 1. / 4. * int.ΔQ1[j] + 3. / 4. * int.ΔQ2[j]
        end

        # if the collocation points are the same for both q and p parts
        simd_mult!(int.ΔP1, int.tG1, int.tΔW)
        simd_mult!(int.ΔP2, int.tG2, int.tΔW)
        if int.params.tab.qdrift.c==int.params.tab.pdrift.c
            for j in 1:int.params.equ.d
                int.solver.x[(int.params.tab.s+i-1)*int.params.equ.d+j] =  Δt_local*(1. / 4. * int.tF1[j] + 3. / 4. * int.tF2[j]) + 1. / 4. * int.ΔP1[j] + 3. / 4. * int.ΔP2[j]
            end
        end
    end


    # If qdrift.c != pdrift.c, then calculating the momenta p at the points pdrift.c[i]
    if int.params.tab.qdrift.c != int.params.tab.pdrift.c

        for i in 1:int.params.tab.s

            # Taking the c[i] from the pdrift tableau.
            Δt_local = int.params.tab.pdrift.c[i]  * int.params.Δt
            int.tΔW .= int.params.tab.pdrift.c[i] .* int.params.ΔW

            simd_mult!(int.ΔQ, int.tB1, int.tΔW)
            simd_mult!(int.ΔP, int.tG1, int.tΔW)
            @. int.Q = int.params.q + 2. / 3. * Δt_local * int.tV1 + 2. / 3. * int.ΔQ
            @. int.P = int.params.p + 2. / 3. * Δt_local * int.tF1 + 2. / 3. * int.ΔP

            t2 = int.params.t + 2. / 3. * Δt_local

            int.params.equ.v(t2, int.Q, int.P, int.tV2)
            int.params.equ.B(t2, int.Q, int.P, int.tB2)
            int.params.equ.f(t2, int.Q, int.P, int.tF2)
            int.params.equ.G(t2, int.Q, int.P, int.tG2)

            # Calculating the Z's and assigning them to the array int.solver.x as initial guesses
            # The guesses for the Y's have already been written to x above
            simd_mult!(int.ΔP1, int.tG1, int.tΔW)
            simd_mult!(int.ΔP2, int.tG2, int.tΔW)
            for j in 1:int.params.equ.d
                int.solver.x[(int.params.tab.s+i-1)*int.params.equ.d+j] =  Δt_local*(1. / 4. * int.tF1[j] + 3. / 4. * int.tF2[j]) + 1. / 4. * int.ΔP1[j] + 3. / 4. * int.ΔP2[j]
            end
        end
    end

end


"""
Integrate PSDE with a stochastic implicit partitioned Runge-Kutta integrator.
Integrating the k-th sample path for the m-th initial condition
"""
function integrate_step!(int::IntegratorSIPRK{DT,TT}, sol::SolutionPSDE{DT,TT,NQ,NW}, k::Int, m::Int, n::Int) where {DT,TT,NQ,NW}
    check_solution_dimension_asserts(sol, k, m, n)

    # set time for nonlinear solver
    int.params.t  = sol.t[0] + (n-1)*int.params.Δt
    int.params.q .= int.q[k, m]
    int.params.p .= int.p[k, m]

    # copy the increments of the Brownian Process
    if NW==1
        #1D Brownian motion, 1 sample path
        int.params.ΔW[1] = sol.W.ΔW[n]
        int.params.ΔZ[1] = sol.W.ΔZ[n]
    elseif NW==2
        #Multidimensional Brownian motion, 1 sample path
        for l = 1:sol.nm
            int.params.ΔW[l] = sol.W.ΔW[l,n]
            int.params.ΔZ[l] = sol.W.ΔZ[l,n]
        end
    elseif NW==3
        #1D or Multidimensional Brownian motion, k-th sample path
        for l = 1:sol.nm
            int.params.ΔW[l] = sol.W.ΔW[l,n,k]
            int.params.ΔZ[l] = sol.W.ΔZ[l,n,k]
        end
    end


    # truncate the increments ΔW with A
    if int.params.A>0
        for i in eachindex(int.params.ΔW)
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
    print_solver_status(int.solver.status, int.solver.params)

    # check if solution contains NaNs or error bounds are violated
    check_solver_status(int.solver.status, int.solver.params)

    # compute the drift vector field and the diffusion matrix at internal stages
    compute_stages!(int.solver.x, int.fcache.Q, int.fcache.P, int.fcache.V, int.fcache.F, int.fcache.B, int.fcache.G, int.fcache.Y, int.fcache.Z, int.params)

    # compute final update
    update_solution!(int.q[k,m], int.p[k,m], int.fcache.V, int.fcache.F, int.fcache.B, int.fcache.G,
                    int.params.tab.qdrift.b, int.params.tab.qdrift.b̂, int.params.tab.qdiff.b, int.params.tab.qdiff.b̂,
                    int.params.tab.pdrift.b, int.params.tab.pdrift.b̂, int.params.tab.pdiff.b, int.params.tab.pdiff.b̂,
                    int.params.Δt, int.params.ΔW, int.Δy, int.Δz)

    # # NOT IMPLEMENTING InitialGuessSDE
    # # # copy solution to initial guess
    # # update!(int.iguess, m, sol.t[0] + n*int.params.Δt, int.q[m])

    # take care of periodic solutions
    cut_periodic_solution!(int.q[k,m], int.params.equ.periodicity)

    # # copy to solution
    set_solution!(sol, int.q[k,m], int.p[k,m], n, k, m)
end
