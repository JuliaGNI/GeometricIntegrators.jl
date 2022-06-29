
# TODO Add VPRKpSecondary tableau.

"Parameters for right-hand side function of Variational Partitioned Runge-Kutta methods."
const ParametersVPRKpSecondary = AbstractParametersVPRK{:vprk_psecondary}


@doc raw"""
Variational partitioned Runge-Kutta integrator with projection on secondary constraint.

The VPRK integrator solves the following system of equations for the internal stages,
```math
\begin{aligned}
Q_{n,i} &= q_{n} + h \sum \limits_{j=1}^{s} a_{ij} \, \big( V_{n,j} + \Lambda_{n,j} \big) , \\
P_{n,i} &= p_{n} + h \sum \limits_{i=1}^{s} \bar{a}_{ij} \, \big( F_{n,j} + \nabla \vartheta (Q_{n,j})^{T} \cdot \Lambda_{n,j} \big) - d_i \lambda , \\
0 &= \sum \limits_{i=1}^{s} d_i V_i , \\
0 &= \sum \limits_{j=1}^{s} \omega_{ij} \Psi_{n,j} ,
\end{aligned}
```
with definitions
```math
\begin{aligned}
P_{n,i} &= \dfrac{\partial L}{\partial v} (Q_{n,i}, V_{n,i}) , \\
F_{n,i} &= \dfrac{\partial L}{\partial q} (Q_{n,i}, V_{n,i}) , \\
\Psi_{n,i} &= \psi(Q_{n,i}, V_{n,i}, P_{n,i}, F_{n,i}) ,
\end{aligned}
```
and update rule
```math
\begin{aligned}
q_{n+1} &= q_{n} + h \sum \limits_{i=1}^{s} b_{i} \, \big( V_{n,i} + \Lambda_{n,i} \big) , \\
p_{n+1} &= p_{n} + h \sum \limits_{i=1}^{s} \bar{b}_{i} \, \big( F_{n,i} + \nabla \vartheta (Q_{n,j})^{T} \cdot \Lambda_{n,j} \big) , \\
0 &= \phi (q_{n+1}, p_{n+1}) ,
\end{aligned}
```
satisfying the symplecticity conditions
```math
\begin{aligned}
b_{i} \bar{a}_{ij} + \bar{b}_{j} a_{ji} &= b_{i} \bar{b}_{j} , &
\bar{b}_i &= b_i ,
\end{aligned}
```
the primary constraint,
```math
\begin{aligned}
\phi(q,p) = p - \vartheta (q) = 0 ,
\end{aligned}
```
at the final solution ``(q_{n+1}, p_{n+1})``,
and super positions of the secondary constraints,
```math
\begin{aligned}
\psi(q,\dot{q},p,\dot{p})
= \dot{p} - \dot{q} \cdot \nabla \vartheta (q)
= \big( \nabla \vartheta (q) - \nabla \vartheta^{T} (q) \big) \cdot \dot{q} - \nabla H (q)
= 0,
\end{aligned}
```
which, evaluated at the internal stages, read
```math
\begin{aligned}
\Psi_{n,j} = \big( \nabla \vartheta (Q_{n,j}) - \nabla \vartheta^{T} (Q_{n,j}) \big) \cdot V_{n,j} - \nabla H (Q_{n,j}) .
\end{aligned}
```
Here, ``\omega`` is a ``(s-1) \times s`` matrix, chosen such that the resulting
method has optimal order.
The vector ``d`` is zero for Gauss-Legendre methods and needs to be chosen
appropriately for Gauss-Lobatto methods (for details see documentation of
VPRK methods).
"""
struct IntegratorVPRKpSecondary{DT, TT, D, S,
                PT <: ParametersVPRKpSecondary{DT,TT},
                ST <: NonlinearSolver{DT},
                IT <: InitialGuessIODE{TT}} <: AbstractIntegratorVPRK{DT,TT,D,S}
    params::PT
    solver::ST
    iguess::IT
    caches::CacheDict{PT}

    function IntegratorVPRKpSecondary(params::ParametersVPRKpSecondary{DT,TT,D,S}, solver::ST, iguess::IT, caches) where {DT,TT,D,S,ST,IT}
        new{DT, TT, D, S, typeof(params), ST, IT}(params, solver, iguess, caches)
    end

    function IntegratorVPRKpSecondary{DT,D}(equations::NamedTuple, tableau::TableauVPRK{TT}, Δt::TT) where {DT,TT,D}
        # get number of stages
        S = tableau.s

        # compute reduction matrix
        ω = zeros(TT, S-1, S)
        for i in 1:(S-1)
            for j in 1:S
                ω[i,j] = tableau.q.b[j] * tableau.q.c[j]^(i-1)
            end
        end

        # create params
        R = convert(Vector{TT}, [1, tableau.R∞])
        params = ParametersVPRKpSecondary{DT,D}(equations, tableau, Δt, NamedTuple{(:R,:ω)}((R,ω)))

        # create cache dict
        caches = CacheDict(params)

        # create solver
        solver = create_nonlinear_solver(DT, 2*D*S, params, caches)

        # create initial guess
        iguess = InitialGuessIODE(get_config(:ig_extrapolation), equations[:v̄], equations[:f̄], Δt)

        # create integrator
        IntegratorVPRKpSecondary(params, solver, iguess, caches)
    end

    function IntegratorVPRKpSecondary(equation::LDAEProblem{DT}, tableau, Δt=tstep(equation); kwargs...) where {DT}
        @assert hassecondary(equation)
        IntegratorVPRKpSecondary{DT, ndims(equation)}(functions(equation), tableau, Δt; kwargs...)
    end
end


function initial_guess!(int::IntegratorVPRKpSecondary{DT}, sol::AtomicSolutionPDAE{DT},
                        cache::IntegratorCacheVPRK{DT}=int.caches[DT]) where {DT}
    for i in eachstage(int)
        evaluate!(int.iguess, sol.q̄, sol.p̄, sol.v̄, sol.f̄,
                              sol.q, sol.p, sol.v, sol.f,
                              cache.q̃, cache.ṽ,
                              tableau(int).q.c[i])

        for k in eachdim(int)
            int.solver.x[ndims(int)*(0*nstages(int)+i-1)+k] = cache.ṽ[k]
            int.solver.x[ndims(int)*(1*nstages(int)+i-1)+k] = 0
        end
    end
end


function compute_stages_vprk!(x, q, v, p, Q, V, Λ, P, F, R, Φ, params)
    # copy x to V
    compute_stages_v_vprk!(x, V, params)

    # copy x to Λ
    compute_stages_λ_vprk!(x, Λ, params)

    # compute Q
    compute_stages_q_vprk!(q, Q, V, Λ, params)

    # compute P and F
    compute_stages_p_vprk!(Q, V, P, F, params)

    # compute p̄, R and Ψ
    compute_projection_vprk!(q, v, p, Q, P, V, F, Λ, R, Φ, params)
end


function compute_stages_q_vprk!(q::Vector{ST}, Q::Vector{Vector{ST}},
                V::Vector{Vector{ST}}, Λ::Vector{Vector{ST}},
                params::ParametersVPRKpSecondary{DT,TT,D,S}) where {ST,DT,TT,D,S}

    for (Qᵢ,Vᵢ,Λᵢ) in zip(Q,V,Λ)
        @assert D == length(Qᵢ) == length(Vᵢ) == length(Λᵢ) == length(q)
    end

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
            Q[i][k] = params.q̄[k] + params.Δt * (y1 + y2) + params.Δt * (y3 + y4)
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
        q[k] = params.q̄[k] + params.Δt * (y1 + y2) + params.Δt * (y3 + y4)
    end
end


function compute_projection_vprk!(q::Vector{ST}, v::Vector{ST}, p::Vector{ST},
                Q::Vector{Vector{ST}}, P::Vector{Vector{ST}},
                V::Vector{Vector{ST}}, F::Vector{Vector{ST}},
                Λ::Vector{Vector{ST}}, R::Vector{Vector{ST}}, Ψ::Vector{Vector{ST}},
                params::ParametersVPRKpSecondary{DT,TT,D,S}) where {ST,DT,TT,D,S}

    # create temporary variables
    local t₀::TT = params.t̄
    local t₁::TT = params.t̄ + params.Δt
    local tᵢ::TT
    # local dH = zeros(ST,D)
    # local Ω  = zeros(ST,D,D)

    # compute p=ϑ(q)
    params.equ[:ϑ](t₁, q, v, p)

    for i in 1:S
        tᵢ = t₀ + params.Δt * params.tab.p.c[i]
        
        # params.equ[:ω](tᵢ, Q[i], V[i], Ω)
        # params.equ[:∇H](tᵢ, Q[i], dH)

        # Ψ[i] .= Ω * V[i] .- dH

        # params.equ[:u](tᵢ, Q[i], Λ[i], U[i])
        params.equ[:g](tᵢ, Q[i], Λ[i], R[i])
        params.equ[:ψ](tᵢ, Q[i], P[i], V[i], F[i], Ψ[i])
    end
end


function compute_rhs_vprk!(b::Vector{ST}, P::Vector{Vector{ST}}, F::Vector{Vector{ST}}, R::Vector{Vector{ST}},
                params::ParametersVPRKpSecondary{DT,TT,D,S}) where {ST,DT,TT,D,S}

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
            b[D*(i-1)+k] = (P[i][k] - params.p̄[k]) - params.Δt * (z1 + z2) - params.Δt * (z3 + z4)
        end
    end
end

function compute_rhs_vprk_projection!(b::Vector{ST}, p::Vector{ST},
                F::Vector{Vector{ST}}, R::Vector{Vector{ST}}, Ψ::Vector{Vector{ST}}, offset::Int,
                params::ParametersVPRKpSecondary{DT,TT,D,S}) where {ST,DT,TT,D,S}

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
    #             b[offset+D*(i-1)+k] += params.pparams[:ω][i,j] * Ψ[j][k]
    #         end
    #     end
    # end
    
    # for k in 1:D
    #     z1 = z2 = z3 = z4 = 0
    #     for j in 1:S
    #         z1 += params.tab.p.b[j] * F[j][k]
    #         z2 += params.tab.p.b̂[j] * F[j][k]
    #         z3 += params.tab.p.b[j] * R[j][k]
    #         z4 += params.tab.p.b̂[j] * R[j][k]
    #     end
    #     b[offset+D*(S-1)+k] = (p[k] - params.p̄[k]) - params.Δt * (z1 + z2) + params.Δt * (z3 + z4)
    # end
end


"Compute stages of variational partitioned Runge-Kutta methods."
function Integrators.function_stages!(x::Vector{ST}, b::Vector{ST},
                params::ParametersVPRKpSecondary{DT,TT,D,S},
                caches::CacheDict) where {ST,DT,TT,D,S}

    # get cache for internal stages
    cache = caches[ST]

    compute_stages_vprk!(x, cache.q̃, cache.ṽ, cache.p̃,
                            cache.Q, cache.V, cache.Λ,
                            cache.P, cache.F, cache.R,
                            cache.Φ, params)

    # compute b = [P-AF-AR]
    compute_rhs_vprk!(b, cache.P, cache.F, cache.R, params)

    # compute b = Φ
    compute_rhs_vprk_projection!(b, cache.p̃, cache.F, cache.R, cache.Φ, D*S, params)

    compute_rhs_vprk_correction!(b, cache.V, params)
end


function Integrators.integrate_step!(int::IntegratorVPRKpSecondary{DT,TT}, sol::AtomicSolutionPDAE{DT,TT},
                                     cache::IntegratorCacheVPRK{DT}=int.caches[DT]) where {DT,TT}
    # update nonlinear solver parameters from cache
    update_params!(int.params, sol)

    # compute initial guess
    initial_guess!(int, sol, cache)

    # reset solution
    reset!(sol, timestep(int))

    # call nonlinear solver
    solve!(int.solver)

    # print solver status
    print_solver_status(int.solver.status, int.solver.params)

    # check if solution contains NaNs or error bounds are violated
    check_solver_status(int.solver.status, int.solver.params)

    # compute final update
    compute_stages_vprk!(int.solver.x,
                         cache.q̃, cache.ṽ, cache.p̃,
                         cache.Q, cache.V, cache.Λ,
                         cache.P, cache.F, cache.R,
                         cache.Φ, int.params)

    # compute unprojected solution
    update_solution!(sol.q, sol.q̃, cache.V, tableau(int).q.b, tableau(int).q.b̂, timestep(int))
    update_solution!(sol.p, sol.p̃, cache.F, tableau(int).p.b, tableau(int).p.b̂, timestep(int))

    # add projection to solution
    update_solution!(sol.q, sol.q̃, cache.Λ, tableau(int).q.b, tableau(int).q.b̂, timestep(int))
    update_solution!(sol.p, sol.p̃, cache.R, tableau(int).p.b, tableau(int).p.b̂, timestep(int))

    # copy solution to initial guess
    update_vector_fields!(int.iguess, sol.t, sol.q, sol.p, sol.v, sol.f)
end
