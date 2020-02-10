
# TODO Add VPRKpSecondary tableau.

"Parameters for right-hand side function of variational partitioned Runge-Kutta methods."
mutable struct ParametersVPRKpSecondary{DT, TT, ET <: VODE{DT,TT}, D, S} <: AbstractParametersVPRK{DT,TT,ET,D,S}
    equ::ET
    tab::TableauVPRK{TT}
    Δt::TT

    ω::Matrix{TT}
    R::Vector{TT}

    t̅::TT
    q̅::Vector{DT}
    p̅::Vector{DT}
end

function ParametersVPRKpSecondary(equ::ET, tab::TableauVPRK{TT}, Δt::TT) where {DT, TT, ET <: VODE{DT,TT}}
    D = equ.d
    S = tab.s

    # compute reduction matrix
    ω = zeros(TT, S-1, S)

    for i in 1:(S-1)
        for j in 1:S
            ω[i,j] = tab.q.b[j] * tab.q.c[j]^(i-1)
        end
    end

    R = convert(Vector{TT}, [1, tab.R∞])
    q = zeros(DT,D)
    p = zeros(DT,D)

    ParametersVPRKpSecondary{DT,TT,ET,D,S}(equ, tab, Δt, ω, R, 0, q, p)
end


@doc raw"""
Variational partitioned Runge-Kutta integrator with projection on secondary constraint.

The VPRK integrator solves the following system of equations for the internal stages,
```math
\begin{align*}
Q_{n,i} &= q_{n} + h \sum \limits_{j=1}^{s} a_{ij} \, \big( V_{n,j} + \Lambda_{n,j} \big) , \\
P_{n,i} &= p_{n} + h \sum \limits_{i=1}^{s} \bar{a}_{ij} \, \big( F_{n,j} + \nabla \vartheta (Q_{n,j})^{T} \cdot \Lambda_{n,j} \big) - d_i \lambda , \\
0 &= \sum \limits_{i=1}^{s} d_i V_i , \\
0 &= \sum \limits_{j=1}^{s} \omega_{ij} \Psi_{n,j} ,
\end{align*}
```
with definitions
```math
\begin{align*}
P_{n,i} &= \dfrac{\partial L}{\partial v} (Q_{n,i}, V_{n,i}) , \\
F_{n,i} &= \dfrac{\partial L}{\partial q} (Q_{n,i}, V_{n,i}) , \\
\Psi_{n,i} &= \psi(Q_{n,i}, V_{n,i}, P_{n,i}, F_{n,i}) ,
\end{align*}
```
and update rule
```math
\begin{align*}
q_{n+1} &= q_{n} + h \sum \limits_{i=1}^{s} b_{i} \, \big( V_{n,i} + \Lambda_{n,i} \big) , \\
p_{n+1} &= p_{n} + h \sum \limits_{i=1}^{s} \bar{b}_{i} \, \big( F_{n,i} + \nabla \vartheta (Q_{n,j})^{T} \cdot \Lambda_{n,j} \big) , \\
0 &= \phi (q_{n+1}, p_{n+1}) ,
\end{align*}
```
satisfying the symplecticity conditions
```math
\begin{align*}
b_{i} \bar{a}_{ij} + \bar{b}_{j} a_{ji} &= b_{i} \bar{b}_{j} , &
\bar{b}_i &= b_i ,
\end{align*}
```
the primary constraint,
```math
\begin{align*}
\phi(q,p) = p - \vartheta (q) = 0 ,
\end{align*}
```
at the final solution ``(q_{n+1}, p_{n+1})``,
and super positions of the secondary constraints,
```math
\begin{align*}
\psi(q,\dot{q},p,\dot{p})
= \dot{p} - \dot{q} \cdot \nabla \vartheta (q)
= \big( \nabla \vartheta (q) - \nabla \vartheta^{T} (q) \big) \cdot \dot{q} - \nabla H (q)
= 0,
\end{align*}
```
which, evaluated at the internal stages, read
```math
\begin{align*}
\Psi_{n,j} = \big( \nabla \vartheta (Q_{n,j}) - \nabla \vartheta^{T} (Q_{n,j}) \big) \cdot V_{n,j} - \nabla H (Q_{n,j}) .
\end{align*}
```
Here, ``\omega`` is a ``(s-1) \times s`` matrix, chosen such that the resulting
method has optimal order.
The vector ``d`` is zero for Gauss-Legendre methods and needs to be chosen
appropriately for Gauss-Lobatto methods (for details see documentation of
VPRK methods).
"""
struct IntegratorVPRKpSecondary{DT, TT, PT <: ParametersVPRKpSecondary{DT,TT},
                                        ST <: NonlinearSolver{DT},
                                        IT <: InitialGuessPODE{DT,TT}, D, S} <: AbstractIntegratorVPRK{DT,TT}
    params::PT
    solver::ST
    iguess::IT
    cache::IntegratorCacheVPRK{DT,D,S}
end

function IntegratorVPRKpSecondary(equation::ET, tableau::TableauVPRK{TT}, Δt::TT) where {DT, TT, ET <: VODE{DT,TT}}
    D = equation.d
    M = equation.n
    S = tableau.s

    # create params
    params = ParametersVPRKpSecondary(equation, tableau, Δt)

    # create solver
    solver = create_nonlinear_solver(DT, 2*D*S, params)

    # create initial guess
    iguess = InitialGuessPODE(get_config(:ig_interpolation), equation, Δt)

    # create cache for internal stage vectors and update vectors
    cache = IntegratorCacheVPRK{DT,D,S}(true)

    # create integrator
    IntegratorVPRKpSecondary{DT, TT, typeof(params), typeof(solver), typeof(iguess), D, S}(
                params, solver, iguess, cache)
end


function compute_stages_vprk!(x, q, p, Q, V, Λ, P, F, R, Φ, params)
    # copy x to V
    compute_stages_v_vprk!(x, V, params)

    # copy x to Λ
    compute_stages_λ_vprk!(x, Λ, params)

    # compute Q
    compute_stages_q_vprk!(q, Q, V, Λ, params)

    # compute p̅, R and Ψ
    compute_projection_vprk!(q, p, Q, V, Λ, R, Φ, params)

    # compute P and F
    compute_stages_p_vprk!(Q, V, P, F, params)
end


function compute_stages_q_vprk!(q::Vector{ST}, Q::Vector{Vector{ST}},
                V::Vector{Vector{ST}}, Λ::Vector{Vector{ST}},
                params::ParametersVPRKpSecondary{DT,TT,ET,D,S}) where {ST,DT,TT,ET,D,S}

    # @assert D == size(Q,1) == size(V,1) == size(Λ,1) == length(q̅)
    @assert S == length(Q) == length(V) == length(Λ)

    local y1::ST
    local y2::ST
    local y3::ST
    local y4::ST

    # compute Q
    for i in 1:S
        for k in 1:D
            y1 = y2 = y3 = y4 = 0
            for j in 1:S
                y1 += params.tab.q.a[i,j] * V[j][k]
                y2 += params.tab.q.â[i,j] * V[j][k]
                y3 += params.tab.q.a[i,j] * Λ[j][k]
                y4 += params.tab.q.â[i,j] * Λ[j][k]
            end
            Q[i][k] = params.q̅[k] + params.Δt * (y1 + y2 + y3 + y4)
        end
    end

    # compute q
    for k in 1:D
        y1 = y2 = y3 = y4 = 0
        for j in 1:S
            y1 += params.tab.q.b[j] * V[j][k]
            y2 += params.tab.q.b̂[j] * V[j][k]
            y3 += params.tab.q.b[j] * Λ[j][k]
            y4 += params.tab.q.b̂[j] * Λ[j][k]
        end
        q[k] = params.q̅[k] + params.Δt * (y1 + y2 + y3 + y4)
    end
end


@generated function compute_projection_vprk!(q::Vector{ST}, p::Vector{ST},
                Q::Vector{Vector{ST}}, V::Vector{Vector{ST}}, Λ::Vector{Vector{ST}},
                R::Vector{Vector{ST}}, Ψ::Vector{Vector{ST}},
                params::ParametersVPRKpSecondary{DT,TT,ET,D,S}
            ) where {ST,DT,TT,ET,D,S}

    # create temporary vectors
    dH = zeros(ST,D)
    Ω = zeros(ST,D,D)

    quote
        local t₀::TT = params.t̅
        local t₁::TT = params.t̅ + params.Δt
        local tᵢ::TT

        # compute p=ϑ(q)
        params.equ.ϑ(t₁, q, p)

        for i in 1:S
            tᵢ = t₀ + params.Δt * params.tab.p.c[i]

            params.equ.g(tᵢ, Q[i], Λ[i], R[i])
            params.equ.Ω(tᵢ, Q[i], $Ω)
            params.equ.∇H(tᵢ, Q[i], $dH)

            # TODO Check if ω() returns Ω or Ω^T -> this decides the sign on dH below
            #      Seems to be Ω^T -> then dH should be subtracted

            Ψ[i] .= $Ω * V[i] .- $dH
        end
    end
end


function compute_rhs_vprk!(b::Vector{ST}, P::Vector{Vector{ST}}, F::Vector{Vector{ST}}, R::Vector{Vector{ST}},
                params::ParametersVPRKpSecondary{DT,TT,ET,D,S}) where {ST,DT,TT,ET,D,S}

    local z1::ST
    local z2::ST
    local z3::ST
    local z4::ST

    # compute b = - [(P-AF-AR)]
    for i in 1:S
        for k in 1:D
            z1 = z2 = z3 = z4 = 0
            for j in 1:S
                z1 += params.tab.p.a[i,j] * F[j][k]
                z2 += params.tab.p.â[i,j] * F[j][k]
                z3 += params.tab.p.a[i,j] * R[j][k]
                z4 += params.tab.p.â[i,j] * R[j][k]
            end
            b[D*(i-1)+k] = (P[i][k] - params.p̅[k]) - params.Δt * (z1 + z2 + z3 + z4)
        end
    end
end

function compute_rhs_vprk_projection!(b::Vector{ST}, p::Vector{ST},
                F::Vector{Vector{ST}}, R::Vector{Vector{ST}}, Ψ::Vector{Vector{ST}}, offset::Int,
                params::ParametersVPRKpSecondary{DT,TT,ET,D,S}
            ) where {ST,DT,TT,ET,D,S}

    local z1::ST
    local z2::ST
    local z3::ST
    local z4::ST

    for i in 1:S
        for k in 1:D
            b[offset+D*(i-1)+k] = Ψ[i][k]
        end
    end

    # TODO # reactivate averaging of secondary constraints and primary constraint of solution

    # for i in 1:(S-1)
    #     for k in 1:D
    #         b[offset+D*(i-1)+k] = 0
    #         for j in 1:S
    #             b[offset+D*(i-1)+k] += params.ω[i,j] * Ψ[j][k]
    #         end
    #     end
    # end
    #
    # for k in 1:D
    #     z1 = z2 = z3 = z4 = 0
    #     for j in 1:S
    #         z1 += params.tab.p.b[j] * F[j][k]
    #         z2 += params.tab.p.b̂[j] * F[j][k]
    #         z3 += params.tab.p.b[j] * R[j][k]
    #         z4 += params.tab.p.b̂[j] * R[j][k]
    #     end
    #     b[offset+D*(S-1)+k] = (p[k] - params.p̅[k]) - params.Δt * (z1 + z2 + z3 + z4)
    # end
end


"Compute stages of variational partitioned Runge-Kutta methods."
@generated function function_stages!(x::Vector{ST}, b::Vector{ST},
                params::ParametersVPRKpSecondary{DT,TT,ET,D,S}
            ) where {ST,DT,TT,ET,D,S}

    cache = IntegratorCacheVPRK{ST, D, S}(true)

    function_stages = quote
        compute_stages_vprk!(x, $cache.q̃, $cache.p̃,
                                $cache.Q, $cache.V, $cache.Λ,
                                $cache.P, $cache.F, $cache.R,
                                $cache.Φ, params)

        # compute b = [P-AF-AR]
        compute_rhs_vprk!(b, $cache.P, $cache.F, $cache.R, params)

        # compute b = Φ
        compute_rhs_vprk_projection!(b, $cache.p̃, $cache.F, $cache.R, $cache.Φ, D*S, params)

        compute_rhs_vprk_correction!(b, $cache.V, params)
    end

    return function_stages
end


function initial_guess!(int::IntegratorVPRKpSecondary, sol::AtomicSolutionPODE)
    for i in eachstage(int)
        evaluate!(int.iguess, sol.q, sol.p, sol.v, sol.f,
                              sol.q̅, sol.p̅, sol.v̅, sol.f̅,
                              int.cache.q̃, int.cache.ṽ,
                              tableau(int).q.c[i])

        for k in eachdim(int)
            int.solver.x[ndims(int)*(0*nstages(int)+i-1)+k] = int.cache.ṽ[k]
            int.solver.x[ndims(int)*(1*nstages(int)+i-1)+k] = 0
        end
    end
end


"Integrate ODE with variational partitioned Runge-Kutta integrator."
function integrate_step!(int::IntegratorVPRKpSecondary{DT,TT}, sol::AtomicSolutionPODE{DT,TT}) where {DT,TT}
    # update nonlinear solver parameters from cache
    update_params!(int.params, sol)

    # compute initial guess
    initial_guess!(int, sol)

    # reset solution
    reset!(sol, timestep(int))

    # call nonlinear solver
    solve!(int.solver)

    # print solver status
    print_solver_status(int.solver.status, int.solver.params)

    # check if solution contains NaNs or error bounds are violated
    check_solver_status(int.solver.status, int.solver.params)

    # compute final update
    compute_stages_vprk!(int.solver.x, int.cache.q̃, int.cache.p̃,
                          int.cache.Q, int.cache.V, int.cache.Λ,
                          int.cache.P, int.cache.F, int.cache.R,
                          int.cache.Φ, int.params)

    # compute unprojected solution
    update_solution!(sol.q, sol.q̃, int.cache.V, tableau(int).q.b, tableau(int).q.b̂, timestep(int))
    update_solution!(sol.p, sol.p̃, int.cache.F, tableau(int).p.b, tableau(int).p.b̂, timestep(int))

    # add projection to solution
    update_solution!(sol.q, sol.q̃, int.cache.Λ, tableau(int).q.b, tableau(int).q.b̂, timestep(int))
    update_solution!(sol.p, sol.p̃, int.cache.R, tableau(int).p.b, tableau(int).p.b̂, timestep(int))

    # copy solution to initial guess
    update!(int.iguess, sol.t, sol.q, sol.p, sol.v, sol.f)
end
