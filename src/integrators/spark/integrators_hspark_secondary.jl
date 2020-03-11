
const TableauHSPARKsecondary{DT} = AbstractTableauSPARK{:hspark_secondary,DT}
const ParametersHSPARKsecondary = AbstractParametersHSPARK{:hspark_secondary}


@doc raw"""
Specialised Partitioned Additive Runge-Kutta integrator for Hamiltonian systems
subject to Dirac constraints with projection on secondary constraint.

This integrator solves the following system of equations for the internal stages,
```math
\begin{align}
Q_{n,i} &= q_{n} + h \sum \limits_{j=1}^{s} a_{ij} V_{n,j} + h \sum \limits_{j=1}^{r} \alpha_{ij} U_{n,j} , & i &= 1, ..., s , \\
P_{n,i} &= p_{n} + h \sum \limits_{j=1}^{s} a_{ij} F_{n,j} + h \sum \limits_{j=1}^{r} \alpha_{ij} G_{n,j} , & i &= 1, ..., s , \\
\tilde{Q}_{n,i} &= q_{n} + h \sum \limits_{j=1}^{s} \tilde{a}_{ij} V_{n,j} + h \sum \limits_{j=1}^{r} \tilde{\alpha}_{ij} U_{n,j} , & i &= 1, ..., r , \\
\tilde{P}_{n,i} &= p_{n} + h \sum \limits_{j=1}^{s} \tilde{a}_{ij} F_{n,j} + h \sum \limits_{j=1}^{r} \tilde{\alpha}_{ij} G_{n,j} , & i &= 1, ..., r , \\
0 &= \sum \limits_{j=1}^{r} \omega_{ij} \tilde{\Phi}_{n,j} , & i &= 1, ..., r-1 ,
\end{align}
```
with definitions
```math
\begin{align}
V_{n,i} &= \hphantom{-} \frac{\partial H}{\partial p} (Q_{n,i}, P_{n,i}) , & i &= 1, ..., s , \\
F_{n,i} &=           -  \frac{\partial H}{\partial q} (Q_{n,i}, P_{n,i}) , & i &= 1, ..., s , \\
U_{n,i} &= \hphantom{-} \frac{\partial \phi}{\partial p} (\tilde{Q}_{n,i}, \tilde{P}_{n,i})^{T} \Lambda_{n,i} , & i &= 1, ..., r , \\
G_{n,i} &=           -  \frac{\partial \phi}{\partial q} (\tilde{Q}_{n,i}, \tilde{P}_{n,i})^{T} \Lambda_{n,i} , & i &= 1, ..., r , \\
\tilde{\Phi}_{n,i} &= \phi(\tilde{Q}_{n,i}, \tilde{P}_{n,i}) , & i &= 1, ..., r ,
\end{align}
```
and update rule
```math
\begin{align}
q_{n+1} &= q_{n} + h \sum \limits_{i=1}^{s} b_{i} V_{n,i} + h \sum \limits_{i=1}^{r} \beta_{i} U_{n,i} , \\
p_{n+1} &= p_{n} + h \sum \limits_{i=1}^{s} b_{i} F_{n,i} + h \sum \limits_{i=1}^{r} \beta_{i} G_{n,i} , \\
0 &= \phi (q_{n+1}, p_{n+1}) .
\end{align}
```
"""
struct IntegratorHSPARKsecondary{DT, TT, ET <: PDAE{DT,TT},
                                         PT <: ParametersHSPARKsecondary{DT,TT},
                                         ST <: NonlinearSolver{DT},
                                         IT <: InitialGuessPODE{DT,TT}, D, S, Σ} <: AbstractIntegratorHSPARK{DT, TT}
    equation::ET
    tableau::TableauHSPARKsecondary{TT}

    params::PT
    solver::ST
    iguess::IT
end


function IntegratorHSPARKsecondary(equation::VDAE{DT,TT},
                                   tableau::TableauVSPARKsecondary{TT}, Δt::TT) where {DT,TT}
    D = equation.d
    S = tableau.s
    Σ = tableau.σ

    # println(tableau)

    N = 2*D*S + 4*D*Σ

    if isdefined(tableau, :d) && length(tableau.d) > 0
        N += D
    end

    # create params
    params = ParametersHSPARKsecondary{DT,D,S,Σ}(equation.ϑ, equation.f, equation.g, equation.g̅, equation.ϕ, equation.ψ, Δt, tableau)

    # create solver
    solver = create_nonlinear_solver(DT, N, params)

    # create initial guess
    iguess = InitialGuessPODE(get_config(:ig_interpolation), equation, Δt)

    # create integrator
    IntegratorHSPARKsecondary{DT, TT, typeof(equation), typeof(params), typeof(solver), typeof(iguess), D, S, Σ}(
                              equation, tableau, params, solver, iguess)
end

@inline pstages(int::IntegratorHSPARKsecondary) = int.tableau.σ
@inline Base.ndims(int::IntegratorHSPARKsecondary{DT,TT,ET,PT,ST,IT,D,S,Σ}) where {DT,TT,ET,PT,ST,IT,D,S,Σ} = D


function compute_stages!(x::Vector{ST}, cache::IntegratorCacheSPARK{ST,TT,D,S,Σ},
                                        params::ParametersHSPARKsecondary{DT,TT,D,S,Σ}) where {ST,DT,TT,D,S,Σ}
    local t::TT

    for i in 1:S
        for k in 1:D
            # copy x to Y, Z
            cache.Yi[i][k] = x[2*(D*(i-1)+k-1)+1]
            cache.Zi[i][k] = x[2*(D*(i-1)+k-1)+2]

            # compute Q and P
            cache.Qi[i][k] = params.q[k] + params.Δt * cache.Yi[i][k]
            cache.Pi[i][k] = params.p[k] + params.Δt * cache.Zi[i][k]
        end

        # compute f(X)
        t = params.t + params.Δt * params.tab.p.c[i]
        params.f_v(t, cache.Qi[i], cache.Pi[i], cache.Vi[i])
        params.f_f(t, cache.Qi[i], cache.Pi[i], cache.Fi[i])
    end

    for i in 1:Σ
        for k in 1:D
            # copy y to Y, Z and Λ
            cache.Yp[i][k] = x[2*D*S+4*(D*(i-1)+k-1)+1]
            cache.Zp[i][k] = x[2*D*S+4*(D*(i-1)+k-1)+2]
            cache.Up[i][k] = x[2*D*S+4*(D*(i-1)+k-1)+3]
            cache.Λp[i][k] = x[2*D*S+4*(D*(i-1)+k-1)+4]

            # compute Q and V
            cache.Qp[i][k] = params.q[k] + params.Δt * cache.Yp[i][k]
            cache.Pp[i][k] = params.p[k] + params.Δt * cache.Zp[i][k]
        end

        # compute f(X)
        t = params.t + params.Δt * params.tab.p̃.c[i]
        params.f_v(t, cache.Qp[i], cache.Pp[i], cache.Vp[i])
        params.f_f(t, cache.Qp[i], cache.Pp[i], cache.Fp[i])

        params.f_g(t, cache.Qp[i], cache.Pp[i], cache.Up[i], cache.Gp[i])
        params.f_g(t, cache.Qp[i], cache.Pp[i], cache.Λp[i], cache.G̅p[i])

        params.f_ϕ(t, cache.Qp[i], cache.Pp[i], cache.Φp[i])
        params.f_ψ(t, cache.Qp[i], cache.Pp[i], cache.Vp[i], cache.Fp[i], cache.Ψp[i])
    end

    # if length(params.tab.d) > 0
    #     for k in 1:D
    #         cache.μ[k] = x[2*D*S+4*D*Σ+k]
    #     end
    # end

    # compute q and p
    cache.q̃ .= params.q
    cache.p̃ .= params.p
    for i in 1:S
        cache.q̃ .+= params.Δt .* params.tab.q.b[1][i] .* cache.Vi[i]
        cache.p̃ .+= params.Δt .* params.tab.p.b[1][i] .* cache.Fi[i]
    end
    for i in 1:Σ
        cache.q̃ .+= params.Δt .* params.tab.q.b[2][i] .* cache.Up[i]
        cache.q̃ .+= params.Δt .* params.tab.q.b[3][i] .* cache.Λp[i]
        cache.p̃ .+= params.Δt .* params.tab.p.b[2][i] .* cache.Gp[i]
        cache.p̃ .+= params.Δt .* params.tab.p.b[3][i] .* cache.G̅p[i]
    end

    # compute ϕ(q,p)
    t = params.t + params.Δt
    params.f_ϕ(t, cache.q̃, cache.p̃, cache.ϕ̃)
end


"Compute stages of specialised partitioned additive Runge-Kutta methods for variational systems."
@generated function Integrators.function_stages!(y::Vector{ST}, b::Vector{ST}, params::ParametersHSPARKsecondary{DT,TT,D,S,Σ}) where {ST,DT,TT,D,S,Σ}
    cache = IntegratorCacheSPARK{ST,TT,D,S,Σ}()

    quote
        compute_stages!(y, $cache, params)

        # compute b = - [(Y-AV-AU), (Z-AF-AG), Φ]
        for i in 1:S
            for k in 1:D
                b[2*(D*(i-1)+k-1)+1] = - $cache.Yi[i][k]
                b[2*(D*(i-1)+k-1)+2] = - $cache.Zi[i][k]
                for j in 1:S
                    b[2*(D*(i-1)+k-1)+1] += params.tab.q.a[1][i,j] * $cache.Vi[j][k]
                    b[2*(D*(i-1)+k-1)+2] += params.tab.p.a[1][i,j] * $cache.Fi[j][k]
                end
                for j in 1:Σ
                    b[2*(D*(i-1)+k-1)+1] += params.tab.q.a[2][i,j] * $cache.Up[j][k]
                    b[2*(D*(i-1)+k-1)+1] += params.tab.q.a[3][i,j] * $cache.Λp[j][k]
                    b[2*(D*(i-1)+k-1)+2] += params.tab.p.a[2][i,j] * $cache.Gp[j][k]
                    b[2*(D*(i-1)+k-1)+2] += params.tab.p.a[3][i,j] * $cache.G̅p[j][k]
                end
            end
        end

        # compute b = - [(Y-AV-AU), (Z-AF-AG), Φ, ωΨ]
        for i in 1:Σ
            for k in 1:D
                b[2*D*S+4*(D*(i-1)+k-1)+1] = - $cache.Yp[i][k]
                b[2*D*S+4*(D*(i-1)+k-1)+2] = - $cache.Zp[i][k]
                b[2*D*S+4*(D*(i-1)+k-1)+3] = - $cache.Φp[i][k]
                # b[2*D*S+4*(D*(i-1)+k-1)+4] = - $cache.Ψp[i][k]
                b[2*D*S+4*(D*(i-1)+k-1)+4] = 0
                for j in 1:S
                    b[2*D*S+4*(D*(i-1)+k-1)+1] += params.tab.q̃.a[1][i,j] * $cache.Vi[j][k]
                    b[2*D*S+4*(D*(i-1)+k-1)+2] += params.tab.p̃.a[1][i,j] * $cache.Fi[j][k]
                end
                for j in 1:Σ
                    b[2*D*S+4*(D*(i-1)+k-1)+1] += params.tab.q̃.a[2][i,j] * $cache.Up[j][k]
                    b[2*D*S+4*(D*(i-1)+k-1)+1] += params.tab.q̃.a[3][i,j] * $cache.Λp[j][k]
                    b[2*D*S+4*(D*(i-1)+k-1)+2] += params.tab.p̃.a[2][i,j] * $cache.Gp[j][k]
                    b[2*D*S+4*(D*(i-1)+k-1)+2] += params.tab.p̃.a[3][i,j] * $cache.G̅p[j][k]
                end
                for j in 1:Σ
                    b[2*D*S+4*(D*(i-1)+k-1)+4] -= params.tab.ω[i,j] * $cache.Ψp[j][k]
                end
                b[2*D*S+4*(D*(i-1)+k-1)+4] -= params.tab.ω[i,Σ+1] * $cache.ϕ̃[k]
            end
        end

        # if length(params.tab.d) > 0
        #     for i in 1:Σ
        #         for k in 1:D
        #             b[2*D*S+4*(D*(i-1)+k-1)+2] -= $cache.μ[k] * params.tab.d[i] / params.tab.p.b[2][i]
        #         end
        #     end
        #
        #     for k in 1:D
        #         b[2*D*S+4*D*Σ+k] = 0
        #         for i in 1:Σ
        #             b[2*D*S+4*D*Σ+k] -= $cache.Vp[i][k] * params.tab.d[i]
        #         end
        #     end
        # end
    end
end


function initial_guess!(int::IntegratorHSPARKsecondary, cache::IntegratorCacheSPARK)
    for i in 1:nstages(int)
        evaluate!(int.iguess, cache.q, cache.p, cache.v, cache.f,
                              cache.q̅, cache.p̅, cache.v̅, cache.f̅,
                              cache.q̃, cache.p̃, cache.ṽ, cache.f̃,
                              tableau(int).q.c[i], tableau(int).p.c[i])

        for k in 1:ndims(int)
            int.solver.x[2*(ndims(int)*(i-1)+k-1)+1] = (cache.q̃[k] - cache.q[k])/timestep(int)
            int.solver.x[2*(ndims(int)*(i-1)+k-1)+2] = (cache.p̃[k] - cache.p[k])/timestep(int)
        end
    end

    for i in 1:pstages(int)
        evaluate!(int.iguess, cache.q, cache.p, cache.v, cache.f,
                              cache.q̅, cache.p̅, cache.v̅, cache.f̅,
                              cache.q̃, cache.p̃, cache.ṽ, cache.f̃,
                              tableau(int).q̃.c[i], tableau(int).p̃.c[i])

        for k in 1:ndims(int)
            int.solver.x[2*ndims(int)*nstages(int)+4*(ndims(int)*(i-1)+k-1)+1] = (cache.q̃[k] - cache.q[k])/timestep(int)
            int.solver.x[2*ndims(int)*nstages(int)+4*(ndims(int)*(i-1)+k-1)+2] = (cache.p̃[k] - cache.p[k])/timestep(int)
            int.solver.x[2*ndims(int)*nstages(int)+4*(ndims(int)*(i-1)+k-1)+3] = 0
            int.solver.x[2*ndims(int)*nstages(int)+4*(ndims(int)*(i-1)+k-1)+4] = 0
        end
    end

    if length(tableau(int).d) > 0
        for k in 1:ndims(int)
            int.solver.x[2*ndims(int)*nstages(int)+4*ndims(int)*pstages(int)+k] = 0
        end
    end
end


function update_solution!(int::IntegratorHSPARKsecondary{DT,TT}, sol::AtomicSolutionPDAE{DT,TT}) where {DT,TT}
    # compute final update
    update_solution!(sol.q, sol.q̃, int.cache.Vi, int.params.tab.q.b[1], timestep(int))
    update_solution!(sol.p, sol.p̃, int.cache.Fi, int.params.tab.p.b[1], timestep(int))

    # compute projection
    update_solution!(sol.q, sol.q̃, int.cache.Up, int.params.tab.q.b[2], timestep(int))
    update_solution!(sol.q, sol.q̃, int.cache.Λp, int.params.tab.q.b[3], timestep(int))
    update_solution!(sol.p, sol.p̃, int.cache.Gp, int.params.tab.p.b[2], timestep(int))
    update_solution!(sol.p, sol.p̃, int.cache.G̅p, int.params.tab.p.b[3], timestep(int))
end
