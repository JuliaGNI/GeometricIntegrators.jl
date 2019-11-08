
# TODO Add VPRKpSecondary tableau.

"Parameters for right-hand side function of variational partitioned Runge-Kutta methods."
mutable struct ParametersVPRKpSecondary{DT, TT, ET <: VODE{DT,TT}, D, S} <: AbstractParametersVPRK{DT,TT,ET,D,S}
    equ::ET
    tab::TableauVPRK{TT}
    Δt::TT

    ω::Matrix{TT}
    R::Vector{TT}

    t::TT
    q::Vector{DT}
    p::Vector{DT}
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


function compute_stages_vprk!(x, q̅, p̅, Q, V, Λ, P, F, R, Φ, params)
    # copy x to V
    compute_stages_v_vprk!(x, V, params)

    # copy x to Λ
    compute_stages_λ_vprk!(x, Λ, params)

    # compute Q
    compute_stages_q_vprk!(q̅, Q, V, Λ, params)

    # compute p̅, R and Ψ
    compute_projection_vprk!(q̅, p̅, Q, V, Λ, R, Φ, params)

    # compute P and F
    compute_stages_p_vprk!(Q, V, P, F, params)
end


function compute_stages_q_vprk!(q̅::Vector{ST}, Q::Matrix{ST}, V::Matrix{ST}, Λ::Matrix{ST},
                params::ParametersVPRKpSecondary{DT,TT,ET,D,S}) where {ST,DT,TT,ET,D,S}

    @assert D == size(Q,1) == size(V,1) == size(Λ,1) == length(q̅)
    @assert S == size(Q,2) == size(V,2) == size(Λ,2)

    local y1::ST
    local y2::ST
    local y3::ST
    local y4::ST

    # compute Q
    for i in 1:S
        for k in 1:D
            y1 = y2 = y3 = y4 = 0
            for j in 1:S
                y1 += params.tab.q.a[i,j] * V[k,j]
                y2 += params.tab.q.â[i,j] * V[k,j]
                y3 += params.tab.q.a[i,j] * Λ[k,j]
                y4 += params.tab.q.â[i,j] * Λ[k,j]
            end
            Q[k,i] = params.q[k] + params.Δt * (y1 + y2 + y3 + y4)
        end
    end

    # compute q̅
    for k in 1:D
        y1 = y2 = y3 = y4 = 0
        for j in 1:S
            y1 += params.tab.q.b[j] * V[k,j]
            y2 += params.tab.q.b̂[j] * V[k,j]
            y3 += params.tab.q.b[j] * Λ[k,j]
            y4 += params.tab.q.b̂[j] * Λ[k,j]
        end
        q̅[k] = params.q[k] + params.Δt * (y1 + y2 + y3 + y4)
    end
end


@generated function compute_projection_vprk!(q̅::Vector{ST}, p̅::Vector{ST},
                Q::Matrix{ST}, V::Matrix{ST}, Λ::Matrix{ST}, R::Matrix{ST}, Ψ::Matrix{ST},
                params::ParametersVPRKpSecondary{DT,TT,ET,D,S}
            ) where {ST,DT,TT,ET,D,S}

    # create temporary vectors
    tQ = zeros(ST,D)
    tV = zeros(ST,D)
    tΛ = zeros(ST,D)
    tR = zeros(ST,D)
    tΨ = zeros(ST,D)
    tH = zeros(ST,D)
    tΩ = zeros(ST,D,D)

    quote
        local t₀::TT = params.t
        local t₁::TT = params.t + params.Δt
        local tᵢ::TT
        local v = zeros(q̅)

        # compute p̅=ϑ(q̅)
        params.equ.ϑ(t₁, q̅, v, p̅)

        for i in 1:S
            simd_copy_xy_first!($tQ, Q, i)
            simd_copy_xy_first!($tV, V, i)
            simd_copy_xy_first!($tΛ, Λ, i)

            tᵢ = t₀ + params.Δt * params.tab.p.c[i]

            params.equ.g(tᵢ, $tQ, $tΛ, $tR)
            params.equ.ω(tᵢ, $tQ, $tΩ)
            params.equ.dH(tᵢ, $tQ, $tH)

            # TODO Check if ω() returns Ω or Ω^T -> this decides the sign on dH below
            #      Seems to be Ω^T -> then dH should be subtracted

            simd_mult!($tΨ, $tΩ, $tV)

            $tΨ .-= $tH

            simd_copy_yx_first!($tR, R, i)
            simd_copy_yx_first!($tΨ, Ψ, i)
        end
    end
end


function compute_rhs_vprk!(b::Vector{ST}, P::Matrix{ST}, F::Matrix{ST}, R::Matrix{ST},
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
                z1 += params.tab.p.a[i,j] * F[k,j]
                z2 += params.tab.p.â[i,j] * F[k,j]
                z3 += params.tab.p.a[i,j] * R[k,j]
                z4 += params.tab.p.â[i,j] * R[k,j]
            end
            b[D*(i-1)+k] = (P[k,i] - params.p[k]) - params.Δt * (z1 + z2 + z3 + z4)
        end
    end
end

function compute_rhs_vprk_projection!(b::Vector{ST}, p̅::Vector{ST},
                F::Matrix{ST}, R::Matrix{ST}, Ψ::Matrix{ST}, offset::Int,
                params::ParametersVPRKpSecondary{DT,TT,ET,D,S}
            ) where {ST,DT,TT,ET,D,S}

    local ψ::ST
    local z1::ST
    local z2::ST
    local z3::ST
    local z4::ST

    # for i in 1:S
    #     for k in 1:D
    #         b[offset+D*(i-1)+k] = Ψ[k,i]
    #     end
    # end

    for i in 1:(S-1)
        for k in 1:D
            ψ = 0
            for j in 1:S
                ψ += params.ω[i,j] * Ψ[k,j]
            end
            b[offset+D*(i-1)+k] = ψ
        end
    end

    for k in 1:D
        z1 = z2 = z3 = z4 = 0
        for j in 1:S
            z1 += params.tab.p.b[j] * F[k,j]
            z2 += params.tab.p.b̂[j] * F[k,j]
            z3 += params.tab.p.b[j] * R[k,j]
            z4 += params.tab.p.b̂[j] * R[k,j]
        end
        b[offset+D*(S-1)+k] = (p̅[k] - params.p[k]) - params.Δt * (z1 + z2 + z3 + z4)
    end
end


"Compute stages of variational partitioned Runge-Kutta methods."
@generated function function_stages!(x::Vector{ST}, b::Vector{ST},
                params::ParametersVPRKpSecondary{DT,TT,ET,D,S}
            ) where {ST,DT,TT,ET,D,S}

    scache = NonlinearFunctionCacheVPRK{ST}(D,S)
    pcache = NonlinearFunctionCacheVPRKprojection{ST}(D,S)

    function_stages = quote
        compute_stages_vprk!(x, $pcache.q̅, $pcache.p̅,
                                $scache.Q, $scache.V, $pcache.Λ,
                                $scache.P, $scache.F, $pcache.R,
                                $pcache.Φ, params)

        # compute b = [P-AF-AR]
        compute_rhs_vprk!(b, $scache.P, $scache.F, $pcache.R, params)

        # compute b = Φ
        compute_rhs_vprk_projection!(b, $pcache.p̅, $scache.F, $pcache.R, $pcache.Φ, D*S, params)

        compute_rhs_vprk_correction!(b, $scache.V, params)
    end

    return function_stages
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
                                        IT <: InitialGuessPODE{DT,TT}} <: AbstractIntegratorVPRK{DT,TT}
    params::PT
    solver::ST
    iguess::IT

    scache::NonlinearFunctionCacheVPRK{DT}
    pcache::NonlinearFunctionCacheVPRKprojection{DT}

    q::Vector{Vector{TwicePrecision{DT}}}
    p::Vector{Vector{TwicePrecision{DT}}}
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
    scache = NonlinearFunctionCacheVPRK{DT}(D,S)
    pcache = NonlinearFunctionCacheVPRKprojection{DT}(D,S)

    # create solution vectors
    q = create_solution_vector(DT, D, M)
    p = create_solution_vector(DT, D, M)

    # create integrator
    IntegratorVPRKpSecondary{DT, TT, typeof(params), typeof(solver), typeof(iguess)}(
                params, solver, iguess, scache, pcache, q, p)
end


function initial_guess!(int::IntegratorVPRKpSecondary{DT,TT}, m::Int) where {DT,TT}
    for i in 1:int.params.tab.s
        evaluate!(int.iguess, m, int.scache.y, int.scache.z, int.scache.v, int.params.tab.q.c[i], int.params.tab.p.c[i])
        for k in 1:int.params.equ.d
            int.solver.x[int.params.equ.d*(0*int.params.tab.s+i-1)+k] = int.scache.v[k]
            int.solver.x[int.params.equ.d*(1*int.params.tab.s+i-1)+k] = 0
        end
    end
end


"Integrate ODE with variational partitioned Runge-Kutta integrator."
function integrate_step!(int::IntegratorVPRKpSecondary{DT,TT}, sol::SolutionPDAE{DT,TT}, m::Int, n::Int) where {DT,TT}
    # check if m and n are compatible with solution dimensions
    check_solution_dimension_asserts(sol, m, n)

    # set time and solution for nonlinear solver
    int.params.t = sol.t[0] + (n-1)*int.params.Δt
    int.params.q .= int.q[m]
    int.params.p .= int.p[m]

    # compute initial guess
    initial_guess!(int, m)

    # call nonlinear solver
    solve!(int.solver)

    # print solver status
    print_solver_status(int.solver.status, int.solver.params, n)

    # check if solution contains NaNs or error bounds are violated
    check_solver_status(int.solver.status, int.solver.params, n)

    # compute final update
    compute_stages_vprk!(int.solver.x, int.pcache.q̅, int.pcache.p̅,
                         int.scache.Q, int.scache.V, int.pcache.Λ,
                         int.scache.P, int.scache.F, int.pcache.R,
                         int.pcache.Φ, int.params)

    # compute unprojected solution
    update_solution!(int.q[m], int.scache.V, int.params.tab.q.b, int.params.tab.q.b̂, int.params.Δt)
    update_solution!(int.p[m], int.scache.F, int.params.tab.p.b, int.params.tab.p.b̂, int.params.Δt)

    # add projection to solution
    update_solution!(int.q[m], int.pcache.Λ, int.params.tab.q.b, int.params.tab.q.b̂, int.params.Δt)
    update_solution!(int.p[m], int.pcache.R, int.params.tab.p.b, int.params.tab.p.b̂, int.params.Δt)

    # copy solution to initial guess for next time step
    update!(int.iguess, m, sol.t[0] + n*int.params.Δt, int.q[m], int.p[m])

    # take care of periodic solutions
    cut_periodic_solution!(int.q[m], int.params.equ.periodicity)

    # copy to solution
    copy_solution!(sol, int.q[m], int.p[m], int.pcache.λ, n, m)
end
