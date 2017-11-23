
"Parameters for right-hand side function of variational partitioned Runge-Kutta methods."
mutable struct NonlinearFunctionParametersVPRKpSecondary{DT, TT, ET <: VODE{DT,TT}, D, S} <: AbstractNonlinearFunctionParametersVPRK{DT,TT,ET,D,S}
    equ::ET
    tab::TableauVPRK{TT}
    Δt::TT

    ω::Matrix{TT}
    R::Vector{TT}

    t::TT
    q::Vector{DT}
    p::Vector{DT}
end

function NonlinearFunctionParametersVPRKpSecondary(equ::ET, tab::TableauVPRK{TT}, Δt::TT) where {DT, TT, ET <: VODE{DT,TT}}
    # compute reduction matrix
    ω = zeros(TT, S-1, S)

    for i in 1:(S-1)
        for j in 1:S
            ω[i,j] = tab.b[j] * tab.c[j]^(i-1)
        end
    end

    R = convert(Vector{TT}, [1, tab.R∞])
    q = zeros(DT, equ.d)
    p = zeros(DT, equ.d)

    NonlinearFunctionParametersVPRKpMidpoint{DT, TT, ET, equ.d, tab.s}(equ, tab, Δt, ω, R, 0, q, p)
end


function compute_stages_vprk!(x, q̅, p̅, λ, Q, V, Λ, U, P, F, R, G, Φ, params)
    # copy x to V
    compute_stages_v_vprk!(x, V, params)

    # copy x to Λ
    compute_stages_λ_vprk!(x, Λ, params)

    # compute U, G and p̅
    compute_projection_vprk!(x, q̅, p̅, λ, V, U, G, params)

    # compute Q
    compute_stages_q_vprk!(Q, V, U, params)

    # compute R and Φ
    compute_internal_projection_vprk!(x, Q, V, Λ, R, Φ, params)

    # compute P and F
    compute_stages_p_vprk!(Q, V, P, F, params)
end


@generated function compute_projection_vprk!(x::Vector{ST},
                q̅::Vector{ST}, p̅::Vector{ST}, λ::Vector{ST},
                V::Matrix{ST}, U::Matrix{ST}, G::Matrix{ST},
                params::NonlinearFunctionParametersVPRKpSecondary{DT,TT,ET,D,S}
            ) where {ST,DT,TT,ET,D,S}

    # create temporary vectors
    tG = zeros(ST,D)

    compute_projection_vprk = quote
        local t₀::TT = params.t
        local t₁::TT = params.t + params.Δt
        local y::ST

        # compute q̅
        for k in 1:D
            y = params.R[1] * U[k,1] + params.R[2] * U[k,2]
            for j in 1:S
                y += params.tab.q.b[j] * V[k,j]
            end
            q̅[k] = params.q[k] + params.Δt * y
        end

        # compute U=λ at tₙ
        params.equ.v(t₀, params.q, params.p, λ)
        simd_copy_yx_first!(λ, U, 1)

        # compute G=g(q,λ) at tₙ
        params.equ.g(t₀, params.q, λ, $tG)
        simd_copy_yx_first!($tG, G, 1)

        # compute U=λ at tₙ+₁
        params.equ.v(t₁, q̅, params.p, λ)
        simd_copy_yx_first!(λ, U, 2)

        # compute G=g(q,λ) at tₙ+₁
        params.equ.g(t₁, q̅, λ, $tG)
        simd_copy_yx_first!($tG, G, 2)

        # compute p̅=α(q̅)
        params.equ.α(t₁, q̅, λ, p̅)
    end

    return compute_projection_vprk
end


@generated function compute_internal_projection_vprk!(
                x::Vector{ST}, Q::Matrix{ST}, V::Matrix{ST},
                Λ::Matrix{ST}, R::Matrix{ST}, Φ::Matrix{ST},
                params::NonlinearFunctionParametersVPRKpSecondary{DT,TT,ET,D,S}
            ) where {ST,DT,TT,ET,D,S}

    # create temporary vectors
    tQ = zeros(ST,D)
    tΛ = zeros(ST,D)
    tV = zeros(ST,D)
    tR = zeros(ST,D)
    tΦ = zeros(ST,D)
    tH = zeros(ST,D)
    tΩ = zeros(ST,D,D)

    compute_internal_projection_vprk = quote
        for i in 1:S
            simd_copy_xy_first!($tQ, Q, i)
            simd_copy_xy_first!($tV, V, i)
            simd_copy_xy_first!($tΛ, Λ, i)

            params.equ.ω(params.t, $tQ, $tΩ)
            params.equ.dH(params.t, $tQ, $tH)

            simd_mult!($tR, $tΩ, $tΛ)
            simd_mult!($tΦ, $tΩ, $tV)

            $tΦ .-= $tH
            # simd_axpy!(-one(ST), $tH, $tΦ)

            simd_copy_yx_first!($tR, R, i)
            simd_copy_yx_first!($tΦ, Φ, i)
        end
    end

    return compute_internal_projection_vprk
end


function compute_rhs_vprk_symplectic_projection!(b::Vector{ST}, p̅::Vector{ST},
                F::Matrix{ST}, R::Matrix{ST}, G::Matrix{ST}, Φ::Matrix{ST},
                offset::Int, params::AbstractNonlinearFunctionParametersVPRK{DT,TT,ET,D,S}
            ) where {ST,DT,TT,ET,D,S}

    local z::ST

    # for i in 1:S
    #     for k in 1:D
    #         b[offset+D*(i-1)+k] = 0
    #         for j in 1:S
    #             b[offset+D*(i-1)+k] -= Φ[k,j]
    #         end
    #     end
    # end
    #
    # for i in 1:(S-1)
    #     for k in 1:D
    #         b[offset+D*(i-1)+k] = 0
    #         for j in 1:S
    #             b[offset+D*(i-1)+k] -= params.ω_λ[i,j] * Φ[k,j]
    #         end
    #     end
    # end

    for k in 1:D
        z = params.R[1] * G[k,1] + params.R[2] * G[k,2]
        for j in 1:S
            z += params.tab.p.b[j] * (F[k,j] + R[k,j])
        end
        b[offset+D*(S-1)+k] = - (p̅[k] - params.p[k]) + params.Δt * z
    end
end


"Compute stages of variational partitioned Runge-Kutta methods."
@generated function function_stages!(x::Vector{ST}, b::Vector{ST},
                params::NonlinearFunctionParametersVPRKpSecondary{DT,TT,ET,D,S}
            ) where {ST,DT,TT,ET,D,S}

    scache = NonlinearFunctionCacheVPRK{ST}(D,S)
    pcache = NonlinearFunctionCacheVPRKprojection{ST}(D,S)

    function_stages = quote
        compute_stages_vprk!(x, $pcache.q̅, $pcache.p̅, $pcache.λ,
                                $scache.Q, $scache.V, $pcache.Λ,
                                $pcache.U, $scache.P, $scache.F,
                                $pcache.R, $pcache.G, $pcache.Φ, params)

        # compute b = - [P-AF-U]
        compute_rhs_vprk!(b, $scache.P, $scache.F, $pcache.R, $pcache.G, params)

        # compute b = - [p-bF-G]
        compute_rhs_vprk_symplectic_projection!(b, $pcache.p̅, $scache.F, $pcache.R, $pcache.G, $pcache.Φ, D*S, params)

        compute_rhs_vprk_correction!(b, $scache.V, params)
    end

    return function_stages
end


"""
Variational partitioned Runge-Kutta integrator with projection on secondary constraint.

```math
\\begin{align*}
P_{n,i} &= \\dfrac{\\partial L}{\\partial v} (Q_{n,i}, V_{n,i}) , &
Q_{n,i} &= q_{n} + h \\sum \\limits_{j=1}^{s} a_{ij} \\, \\big( V_{n,j} + \\Lambda_{n,j} \\big) , &
q_{n+1} &= q_{n} + h \\sum \\limits_{i=1}^{s} b_{i} \\, \\big( V_{n,i} + \\Lambda_{n,i} \\big) , \\\\
F_{n,i} &= \\dfrac{\\partial L}{\\partial q} (Q_{n,i}, V_{n,i}) , &
P_{n,i} &= p_{n} + h \\sum \\limits_{i=1}^{s} \\bar{a}_{ij} \\, \\big( F_{n,j} + \\nabla \\vartheta (Q_{n,j}) \\cdot \\Lambda_{n,j} \\big) - d_i \\lambda , &
p_{n+1} &= p_{n} + h \\sum \\limits_{i=1}^{s} \\bar{b}_{i} \\, \\big( F_{n,i} + \\nabla \\vartheta (Q_{n,j}) \\cdot \\Lambda_{n,j} \\big) , \\\\
0 &= \\sum \\limits_{i=1}^{s} d_i V_i , &
0 &= \\sum \\limits_{j=1}^{s} \\omega_{ij} \\Psi_{n,j} , &
0 &= \\phi (q_{n+1}, p_{n+1}) ,
\\end{align*}
```
satisfying the symplecticity conditions
```math
\\begin{align*}
b_{i} \\bar{a}_{ij} + \\bar{b}_{j} a_{ji} &= b_{i} \\bar{b}_{j} , &
\\bar{b}_i &= b_i ,
\\end{align*}
```
the primary constraint,
```math
\\begin{align*}
\\phi(q,p) = p - \\vartheta (q) = 0 ,
\\end{align*}
```
at the final solution ``(q_{n+1}, p_{n+1})``,
and super positions of the secondary constraints,
```math
\\begin{align*}
\\psi(q,\\dot{q},p,\\dot{p})
= \\dot{p} - \\dot{q} \\cdot \\nabla \\vartheta (q)
= \\big( \\nabla \\vartheta (q) - \\nabla \\vartheta^{T} (q) \\big) \\cdot \\dot{q} - \\nabla H (q)
= 0,
\\end{align*}
```
at the internal stages,
```math
\\begin{align*}
\\Psi_{n,j} = \\big( \\nabla \\vartheta (Q_{n,j}) - \\nabla \\vartheta^{T} (Q_{n,j}) \\big) \\cdot V_{n,j} - \\nabla H (Q_{n,j}) .
\\end{align*}
```
Here, ``\\omega`` is a ``(s-1) \\times s`` matrix, chosen such that the resulting
method has optimal order.
The vector ``d`` is zero for Gauss-Legendre methods and needs to be chosen
appropriately for Gauss-Lobatto methods (for details see documentation of
VPRK methods).
"""
#
# TODO Fix this and correctly implement what is described above.
#
struct IntegratorVPRKpSecondary{DT, TT, PT <: NonlinearFunctionParametersVPRKpSecondary{DT,TT},
                                        ST <: NonlinearSolver{DT},
                                        IT <: InitialGuessPODE{DT,TT}} <: AbstractIntegratorVPRK{DT,TT}
    params::PT
    solver::ST
    iguess::IT

    scache::NonlinearFunctionCacheVPRK{DT}
    pcache::NonlinearFunctionCacheVPRKprojection{DT}

    q::Vector{Vector{Double{DT}}}
    p::Vector{Vector{Double{DT}}}
end

function IntegratorVPRKpSecondary(equation::ET, tableau::TableauVPRK{TT}, Δt::TT) where {DT, TT, ET <: VODE{DT,TT}}
    D = equation.d
    M = equation.n
    S = tableau.s

    # create params
    params = NonlinearFunctionParametersVPRKpSecondary(equation, tableau, Δt)

    # create solver
    solver = create_nonlinear_solver(DT, 2*D*S, sparams)

    # create initial guess
    iguess = InitialGuessPODE(get_config(:ig_interpolation), equation, Δt)

    # create cache for internal stage vectors and update vectors
    scache = NonlinearFunctionCacheVPRK{DT}(D,S)
    pcache = NonlinearFunctionCacheVPRKprojection{DT}(D,S)

    # create solution vectors
    q = create_solution_vector_double_double(DT, D, M)
    p = create_solution_vector_double_double(DT, D, M)

    # create integrator
    IntegratorVPRKpSecondary{DT, TT, typeof(params), typeof(solver), typeof(iguess)}(
                params, solver, iguess, scache, pcache, q, p)
end


function initial_guess!(int::IntegratorVPRKpSecondary{DT,TT}, m::Int) where {DT,TT}
    for i in 1:params.tab.s
        evaluate!(int.iguess, int.scache.y, int.scache.z, int.scache.v, int.params.tab.q.c[i], int.params.tab.p.c[i])
        for k in 1:int.params.equ.d
            int.solver.x[int.params.equ.d*(0*int.params.tab.s+i-1)+k] = int.scache.v[k]
            int.solver.x[int.params.equ.d*(1*int.params.tab.s+i-1)+k] = 0
        end
    end
end


"Integrate ODE with variational partitioned Runge-Kutta integrator."
function integrate!(int::IntegratorVPRKpSecondary{DT,TT}, sol::SolutionPDAE{DT,TT}, m1::Int, m2::Int) where {DT,TT}
    # check if m and n are compatible with solution dimensions
    check_solution_dimension_asserts(sol, m, n)

    # set time and solution for nonlinear solver
    int.sparams.t = sol.t[0] + (n-1)*int.Δt
    int.sparams.q .= int.q[m]
    int.sparams.p .= int.p[m]

    # compute initial guess
    initial_guess!(int, m)

    # call nonlinear solver
    solve!(int.solver)

    # print solver status
    print_solver_status(int.solver.status, int.solver.params, n)

    # check if solution contains NaNs or error bounds are violated
    check_solver_status(int.solver.status, int.solver.params, n)

    # compute final update
    compute_stages_vprk!(int.solver.x,
                         int.pcache.q̅, int.pcache.p̅, int.pcache.λ,
                         int.scache.Q, int.scache.V, int.pcache.Λ, int.pcache.U,
                         int.scache.P, int.scache.F, int.pcache.R, int.pcache.G,
                         int.pcache.Φ, int.params)

    # compute unprojected solution
    update_solution!(int.q[m], int.scache.V, int.params.tab.q.b, int.params.tab.q.b̂, int.Δt)
    update_solution!(int.p[m], int.scache.F, int.params.tab.p.b, int.params.tab.p.b̂, int.Δt)

    # add projection to solution
    update_solution!(int.q[m], int.pcache.U, int.params.R, int.params.Δt)
    update_solution!(int.p[m], int.pcache.G, int.params.R, int.params.Δt)

    # copy solution to initial guess for next time step
    update!(int.iguess, m, sol.t[0] + n*int.params.Δt, int.q[m], int.p[m])

    # take care of periodic solutions
    cut_periodic_solution!(int.q[m], int.equation.periodicity)

    # copy to solution
    copy_solution!(sol, int.q[m], int.p[m], int.pcache.λ, n, m)
end
