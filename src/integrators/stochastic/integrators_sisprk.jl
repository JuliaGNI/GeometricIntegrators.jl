"""
Holds the tableau of a stochastic implicit split partitioned Runge-Kutta method.

qdrift, pdrift1, pdrift2 hold the RK coefficients for the drift parts,
and qdiff, pdiff1, pdiff2 hold the RK coefficients for the diffusion part of the SDE.

Order of the tableau is not included, because unlike in the deterministic
setting, it depends on the properties of the noise (e.g., the dimension of
the Wiener process and the commutativity properties of the diffusion matrix)

Orders stored in qdrift and qdiff are understood as the classical orders of these methods.

"""
struct TableauSISPRK{T} <: AbstractTableauIRK{T}
    name::Symbol
    s::Int
    qdrift ::CoefficientsRK{T}
    qdiff  ::CoefficientsRK{T}
    pdrift1::CoefficientsRK{T}
    pdrift2::CoefficientsRK{T}
    pdiff1 ::CoefficientsRK{T}
    pdiff2 ::CoefficientsRK{T}

    function TableauSISPRK{T}(name, s, qdrift, qdiff, pdrift1, pdrift2, pdiff1, pdiff2) where {T}
        @assert s == qdrift.s == qdiff.s == pdrift1.s == pdrift2.s == pdiff1.s == pdiff2.s
        new(name, s, qdrift, qdiff, pdrift1, pdrift2, pdiff1, pdiff2)
    end
end

function TableauSISPRK(name::Symbol, qdrift::CoefficientsRK{T},
                                     qdiff::CoefficientsRK{T},
                                     pdrift1::CoefficientsRK{T}, pdrift2::CoefficientsRK{T},
                                     pdiff1::CoefficientsRK{T}, pdiff2::CoefficientsRK{T}) where {T}
    TableauSISPRK{T}(name, qdrift.s, qdrift, qdiff, pdrift1, pdrift2, pdiff1, pdiff2)
end

function TableauSISPRK(name::Symbol, qorder_drift::Int, qa_drift::Matrix{T}, qb_drift::Vector{T}, qc_drift::Vector{T},
                                      qorder_diff::Int, qa_diff::Matrix{T}, qb_diff::Vector{T}, qc_diff::Vector{T},
                                      porder1_drift::Int, pa1_drift::Matrix{T}, pb1_drift::Vector{T}, pc1_drift::Vector{T},
                                      porder2_drift::Int, pa2_drift::Matrix{T}, pb2_drift::Vector{T}, pc2_drift::Vector{T},
                                      porder1_diff::Int, pa1_diff::Matrix{T}, pb1_diff::Vector{T}, pc1_diff::Vector{T},
                                      porder2_diff::Int, pa2_diff::Matrix{T}, pb2_diff::Vector{T}, pc2_diff::Vector{T}) where {T}
    TableauSISPRK{T}(name, length(qc_drift), CoefficientsRK(name, qorder_drift, qa_drift, qb_drift, qc_drift),
                                              CoefficientsRK(name, qorder_diff, qa_diff, qb_diff, qc_diff),
                                              CoefficientsRK(name, porder1_drift, pa1_drift, pb1_drift, pc1_drift),
                                              CoefficientsRK(name, porder2_drift, pa2_drift, pb2_drift, pc2_drift),
                                              CoefficientsRK(name, porder1_diff, pa1_diff, pb1_diff, pc1_diff),
                                              CoefficientsRK(name, porder2_diff, pa2_diff, pb2_diff, pc2_diff))
end

# TODO function readTableauSFIRKFromFile(dir::AbstractString, name::AbstractString)


"""
Parameters for right-hand side function of implicit Runge-Kutta methods.
  A - if positive, the upper bound of the Wiener process increments; if A=0.0, no truncation
"""
mutable struct ParametersSISPRK{DT, TT, ET <: SPSDE{DT,TT}, D, M, S} <: Parameters{DT,TT}
    equ::ET
    tab::TableauSISPRK{TT}
    Δt::TT
    ΔW::Vector{DT}
    ΔZ::Vector{DT}
    A::DT

    t::TT
    q::Vector{DT}
    p::Vector{DT}
end

function ParametersSISPRK(equ::ET, tab::TableauSISPRK{TT}, Δt::TT, ΔW::Vector{DT}, ΔZ::Vector{DT}, A::DT) where {DT, TT, ET <: SPSDE{DT,TT}}
    @assert equ.m == length(ΔW) == length(ΔZ)
    ParametersSISPRK{DT, TT, ET, equ.d, equ.m, tab.s}(equ, tab, Δt, ΔW, ΔZ, A, 0, zeros(DT, equ.d), zeros(DT, equ.d))
end


"""
Structure for holding the internal stages Q, the values of the drift vector
and the diffusion matrix evaluated at the internal stages V=v(Q), B=B(Q),
and the increments Y = Δt*a_drift*v(Q) + a_diff*B(Q)*ΔW
"""
struct IntegratorCacheSISPRK{DT}
    Q ::Vector{Vector{DT}}
    P ::Vector{Vector{DT}}
    V ::Vector{Vector{DT}}
    F1::Vector{Vector{DT}}
    F2::Vector{Vector{DT}}
    B ::Vector{Matrix{DT}}
    G1::Vector{Matrix{DT}}
    G2::Vector{Matrix{DT}}
    Y ::Vector{Vector{DT}}
    Z ::Vector{Vector{DT}}

    v ::Vector{DT}
    f1::Vector{DT}
    f2::Vector{DT}
    b ::Matrix{DT}
    g1::Matrix{DT}
    g2::Matrix{DT}
    y ::Vector{DT}
    z ::Vector{DT}

    Δy::Vector{DT}
    Δz::Vector{DT}

    function IntegratorCacheSISPRK{DT}(D, M, S) where {DT}

        # create internal stage vectors
        Q  = create_internal_stage_vector(DT, D, S)
        P  = create_internal_stage_vector(DT, D, S)
        V  = create_internal_stage_vector(DT, D, S)
        F1 = create_internal_stage_vector(DT, D, S)
        F2 = create_internal_stage_vector(DT, D, S)
        B  = create_internal_stage_vector(DT, D, M, S)
        G1 = create_internal_stage_vector(DT, D, M, S)
        G2 = create_internal_stage_vector(DT, D, M, S)
        Y  = create_internal_stage_vector(DT, D, S)
        Z  = create_internal_stage_vector(DT, D, S)

        # create velocity and update vector
        v  = zeros(DT,D)
        f1 = zeros(DT,D)
        f2 = zeros(DT,D)
        b  = zeros(DT,D,M)
        g1 = zeros(DT,D,M)
        g2 = zeros(DT,D,M)
        y  = zeros(DT,D)
        z  = zeros(DT,D)

        # create temporary arrays
        Δy  = zeros(DT,M)
        Δz  = zeros(DT,M)

        new(Q, P, V, F1, F2, B, G1, G2, Y, Z, v, f1, f2, b, g1, g2, y, z, Δy, Δz)
    end
end


"Stochastic implicit partitioned Runge-Kutta integrator."
struct IntegratorSISPRK{DT, TT, PT <: ParametersSISPRK{DT,TT},
                                ST <: NonlinearSolver{DT}} <: StochasticIntegrator{DT,TT}
    params::PT
    solver::ST
    cache::IntegratorCacheSISPRK{DT}
end

function IntegratorSISPRK(equation::SPSDE{DT,TT}, tableau::TableauSISPRK{TT}, Δt::TT; K::Int=0) where {DT,TT}
    # K - the integer in the bound A = √(2 K Δt |log Δt|) due to Milstein & Tretyakov; K=0 no truncation

    D = equation.d
    M = equation.m
    NS= max(equation.ns,equation.ni)
    S = tableau.s

    # create params
    K==0 ? A = 0.0 : A = sqrt( 2*K*Δt*abs(log(Δt)) )
    params = ParametersSISPRK(equation, tableau, Δt, zeros(DT,M), zeros(DT,M), A)

    # create solver
    solver = create_nonlinear_solver(DT, 2*D*S, params)

    # create cache for internal stage vectors and update vectors
    cache = IntegratorCacheSISPRK{DT}(D, M, S)

    # create integrator
    IntegratorSISPRK{DT, TT, typeof(params), typeof(solver)}(params, solver, cache)
end

@inline equation(integrator::IntegratorSISPRK) = integrator.params.equ
@inline timestep(integrator::IntegratorSISPRK) = integrator.params.Δt
@inline tableau(integrator::IntegratorSISPRK) = integrator.params.tab
@inline nstages(integrator::IntegratorSISPRK)  = nstages(tableau(integrator))
@inline eachstage(integrator::IntegratorSISPRK) = 1:nstages(integrator)
@inline Base.eltype(integrator::IntegratorSISPRK{DT}) where {DT} = DT


function update_params!(int::IntegratorSISPRK, sol::AtomicSolutionPSDE)
    # set time for nonlinear solver and copy previous solution
    int.params.t  = sol.t
    int.params.q .= sol.q
    int.params.p .= sol.p
    int.params.ΔW .= sol.ΔW
    int.params.ΔZ .= sol.ΔZ

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
end


"""
Unpacks the data stored in x = (Y[1][1], Y[1][2], ... Y[1][D], Y[2][1], ..., Z[1][1], Z[1][2], ... Z[1][D], Z[2][1], ...)
into Y, Z, calculates the internal stages Q, P, the values of the RHS
of the SDE ( vi(Q,P), fi(Q,P), Bi(Q,P) and Gi(Q,P) ), and assigns them to V[i], F[i], B[i] and G[i].
Unlike for FIRK, here
Y = Δt a_drift v(Q,P) + a_diff B(Q,P) ΔW,
Z = Δt â1_drift f1(Q,P) + Δt â2_drift f2(Q,P) + â1_diff G1(Q,P) ΔW + â2_diff G2(Q,P) ΔW.
"""
function compute_stages!(x::Vector{ST}, cache::IntegratorCacheSISPRK{ST},
                         params::ParametersSISPRK{DT,TT,ET,D,M,S}) where {ST,DT,TT,ET,D,M,S}

    local tqᵢ::TT       #times for the q internal stages
    local tpᵢ₁::TT      #times for the p internal stages
    local tpᵢ₂::TT      #times for the p internal stages

    Q  = cache.Q
    P  = cache.P
    V  = cache.V
    F1 = cache.F1
    F2 = cache.F2
    B  = cache.B
    G1 = cache.G1
    G2 = cache.G2
    Y  = cache.Y
    Z  = cache.Z

    # copy x to Y, Z and calculate Q, P
    for i in 1:S
        for k in 1:D
            Y[i][k] = x[D*(  i-1)+k]
            Z[i][k] = x[D*(S+i-1)+k]
        end
        Q[i] .= params.q .+ Y[i]
        P[i] .= params.p .+ Z[i]
    end

    # compute Vi = v(Qi,Pi), Fi = f(Qi,Pi), Bi=B(Qi,Pi) and Gi=G(Qi,Pi)
    for i in 1:S
        tqᵢ  = params.t + params.Δt * params.tab.qdrift.c[i]
        tpᵢ₁ = params.t + params.Δt * params.tab.pdrift1.c[i]
        tpᵢ₂ = params.t + params.Δt * params.tab.pdrift2.c[i]

        # calculates v(t,Q,P), f1(t,Q,P), f2(t,Q,P) and assigns to V, F1, F2
        params.equ.v(tqᵢ, Q[i], P[i], V[i])
        params.equ.f1(tpᵢ₁, Q[i], P[i], F1[i])
        params.equ.f2(tpᵢ₂, Q[i], P[i], F2[i])

        # calculates B(t,Q,P), G1(t,Q,P), G2(t,Q,P) and assigns to the matrices B[i], G1[i], G2[i]
        params.equ.B(tqᵢ, Q[i], P[i], B[i])
        params.equ.G1(tpᵢ₁, Q[i], P[i], G1[i])
        params.equ.G2(tpᵢ₂, Q[i], P[i], G2[i])
    end
end

"Compute stages of stochastic implicit split partitioned Runge-Kutta methods."
@generated function Integrators.function_stages!(x::Vector{ST}, b::Vector{ST}, params::ParametersSISPRK{DT,TT,ET,D,M,S}) where {ST,DT,TT,ET,D,M,S}

    cache = IntegratorCacheSISPRK{ST}(D, M, S)

    quote
        compute_stages!(x, $cache, params)

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
                    y1 += params.tab.qdrift.a[i,j]  * $cache.V[j][k]  * params.Δt
                    y2 += params.tab.qdrift.â[i,j]  * $cache.V[j][k]  * params.Δt
                    z1 += params.tab.pdrift1.a[i,j] * $cache.F1[j][k] * params.Δt
                       +  params.tab.pdrift2.a[i,j] * $cache.F2[j][k] * params.Δt
                    z2 += params.tab.pdrift1.â[i,j] * $cache.F1[j][k] * params.Δt
                       +  params.tab.pdrift2.â[i,j] * $cache.F2[j][k] * params.Δt
                    for l in 1:M
                        y1 += params.tab.qdiff.a[i,j]  * $cache.B[j][k,l]  * params.ΔW[l]
                        y2 += params.tab.qdiff.â[i,j]  * $cache.B[j][k,l]  * params.ΔW[l]
                        z1 += params.tab.pdiff1.a[i,j] * $cache.G1[j][k,l] * params.ΔW[l]
                           +  params.tab.pdiff2.a[i,j] * $cache.G2[j][k,l] * params.ΔW[l]
                        z2 += params.tab.pdiff1.â[i,j] * $cache.G1[j][k,l] * params.ΔW[l]
                           +  params.tab.pdiff2.â[i,j] * $cache.G2[j][k,l] * params.ΔW[l]
                    end
                end
                b[D*(  i-1)+k] = - $cache.Y[i][k] + (y1 + y2)
                b[D*(S+i-1)+k] = - $cache.Z[i][k] + (z1 + z2)
            end
        end
    end
end


"""
This function computes initial guesses for Y, Z and assigns them to int.solver.x

For SISPRK we are NOT IMPLEMENTING an InitialGuess.

SIMPLE SOLUTION
The simplest initial guess for Y, Z is 0
"""
function initial_guess!(int::IntegratorSISPRK{DT,TT}, sol::AtomicSolutionPSDE{DT,TT}) where {DT,TT}
    int.solver.x .= 0
end


function update_solution!(sol::AtomicSolutionPSDE{T},
                          cache::IntegratorCacheSISPRK{T}, tab::TableauSISPRK{T},
                          Δt::T, ΔW::Vector{T}) where {T}

    update_solution!(sol, cache.V, cache.F1, cache.F2, cache.B, cache.G1, cache.G2,
                     tab.qdrift.b, tab.qdiff.b, tab.pdrift1.b, tab.pdrift2.b, tab.pdiff1.b, tab.pdiff2.b,
                     Δt, ΔW, cache.Δy, cache.Δz)
    update_solution!(sol, cache.V, cache.F1, cache.F2, cache.B, cache.G1, cache.G2,
                     tab.qdrift.b̂, tab.qdiff.b̂, tab.pdrift1.b̂, tab.pdrift2.b̂, tab.pdiff1.b̂, tab.pdiff2.b̂,
                     Δt, ΔW, cache.Δy, cache.Δz)
end


"""
Integrate PSDE with a stochastic implicit partitioned Runge-Kutta integrator.
"""
function Integrators.integrate_step!(int::IntegratorSISPRK{DT,TT}, sol::AtomicSolutionPSDE{DT,TT}) where {DT,TT}
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
    compute_stages!(int.solver.x, int.cache, int.params)

    # compute final update
    update_solution!(sol, int.cache, int.params.tab, int.params.Δt, int.params.ΔW)
end
