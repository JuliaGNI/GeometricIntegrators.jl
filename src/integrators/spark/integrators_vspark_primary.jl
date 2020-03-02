
const TableauVSPARKprimary = AbstractTableauSPARK{:vspark_primary}
const ParametersVSPARKprimary = AbstractParametersVSPARK{:vspark_primary}


@doc raw"""
Specialised Partitioned Additive Runge-Kutta integrator for Variational systems.

This integrator solves the following system of equations for the internal stages,
```math
\begin{align}
Q_{n,i} &= q_{n} + h \sum \limits_{j=1}^{s} a_{ij} V_{n,j} + h \sum \limits_{j=1}^{r} \alpha_{ij} U_{n,j} , & i &= 1, ..., s , \\
P_{n,i} &= p_{n} + h \sum \limits_{j=1}^{s} a_{ij} F_{n,j} + h \sum \limits_{j=1}^{r} \alpha_{ij} G_{n,j} , & i &= 1, ..., s , \\
\tilde{Q}_{n,i} &= q_{n} + h \sum \limits_{j=1}^{s} \tilde{a}_{ij} V_{n,j} + h \sum \limits_{j=1}^{r} \tilde{\alpha}_{ij} U_{n,j} , & i &= 1, ..., r , \\
\tilde{P}_{n,i} &= p_{n} + h \sum \limits_{j=1}^{s} \tilde{a}_{ij} F_{n,j} + h \sum \limits_{j=1}^{r} \tilde{\alpha}_{ij} G_{n,j} , & i &= 1, ..., r , \\
0 &= \sum \limits_{j=1}^{r} \omega_{ij} \tilde{\Phi}_{n,j} , & i &= 1, ..., r-1 , \\
0 &= \sum \limits_{i=1}^{r} \tilde{d}_i \, \Lambda_{n,i} ,
\end{align}
```
with definitions
```math
\begin{align}
P_{n,i} &= \frac{\partial L}{\partial v} (Q_{n,i}, V_{n,i}) , & i &= 1, ..., s , \\
F_{n,i} &= \frac{\partial L}{\partial q} (Q_{n,i}, V_{n,i}) , & i &= 1, ..., s , \\
U_{n,i} &= \hphantom{-} \frac{\partial \phi}{\partial p} (\tilde{Q}_{n,i}, \tilde{P}_{n,i})^{T} \Lambda_{n,i} , & i &= 1, ..., r , \\
G_{n,i} &=           -  \frac{\partial \phi}{\partial q} (\tilde{Q}_{n,i}, \tilde{P}_{n,i})^{T} \Lambda_{n,i} , & i &= 1, ..., r , \\
\tilde{\Phi}_{n,i} &= \phi(\tilde{Q}_{n,i}, \tilde{P}_{n,i}) , & i &= 1, ..., r ,
\end{align}
```
and update rule
```math
\begin{align}
q_{n+1} &= q_{n} + h \sum \limits_{i=1}^{s} b_{i} V_{n,i} + h \sum \limits_{i=1}^{r} \beta_{i} U_{n,i} , \\
p_{n+1} &= p_{n} + h \sum \limits_{i=1}^{s} b_{i} F_{n,i} + h \sum \limits_{i=1}^{r} \beta_{i} G_{n,i} .
\end{align}
```
"""
struct IntegratorVSPARKprimary{DT, TT, ET <: IDAE{DT,TT},
                                       PT <: ParametersVSPARKprimary{DT,TT},
                                       ST <: NonlinearSolver{DT},
                                       IT <: InitialGuessPODE{DT,TT}, D, S, R} <: AbstractIntegratorVSPARK{DT, TT}
    equation::ET
    tableau::TableauVSPARKprimary{TT}

    params::PT
    solver::ST
    iguess::IT
    cache::IntegratorCacheSPARK{DT,TT,D,S,R}
end

function IntegratorVSPARKprimary(equation::IDAE{DT,TT},
                                 tableau::TableauVSPARKprimary{TT}, Δt::TT) where {DT,TT}
    D = equation.d
    S = tableau.s
    R = tableau.r
    P = tableau.ρ

    # @assert tableau.ρ == tableau.r-1

    N = 2*D*S + 2*D*R

    if isdefined(tableau, :d) && length(tableau.d) > 0
        N += D
    end

    # create params
    params = ParametersVSPARKprimary{DT,D,S,R,P}(equation.f, equation.ϑ, equation.u, equation.g, equation.ϕ, Δt, tableau)

    # create solver
    solver = create_nonlinear_solver(DT, N, params)

    # create initial guess
    iguess = InitialGuessPODE(get_config(:ig_interpolation), equation, Δt)

    # create cache
    cache = IntegratorCacheSPARK{DT,TT,D,S,R}()

    # create integrator
    IntegratorVSPARKprimary{DT, TT, typeof(equation), typeof(params), typeof(solver), typeof(iguess), D, S, R}(
                                        equation, tableau, params, solver, iguess, cache)
end


@inline Base.ndims(int::IntegratorVSPARKprimary{DT,TT,ET,PT,ST,IT,D,S,R}) where {DT,TT,ET,PT,ST,IT,D,S,R} = D


function compute_stages!(x::Vector{ST}, cache::IntegratorCacheSPARK{ST,TT,D,S,R},
                                        params::ParametersVSPARKprimary{DT,TT,D,S,R,P}) where {ST,DT,TT,D,S,R,P}
    local tpᵢ::TT
    local tλᵢ::TT

    # copy x to Vi and Zi
    for i in 1:S
        for k in 1:D
            cache.Vi[i][k] = x[2*(D*(i-1)+k-1)+1]
            cache.Zi[i][k] = x[2*(D*(i-1)+k-1)+2]
        end
    end

    # copy x to Λp and Zp
    for i in 1:R
        for k in 1:D
            cache.Λp[i][k] = x[2*D*S+2*(D*(i-1)+k-1)+1]
            cache.Zp[i][k] = x[2*D*S+2*(D*(i-1)+k-1)+2]
        end
        # tλᵢ = params.t + params.Δt * params.tab.λ.c[i]
        # params.f_u(tλᵢ, cache.Qp[i], cache.Pp[i], cache.Λp[i], cache.Up[i])
        cache.Up[i] .= cache.Λp[i]
    end

    for i in 1:S
        # compute Y
        cache.Yi[i] .= 0
        for j in 1:S
            cache.Yi[i] .+= params.tab.q.a[i,j] .* cache.Vi[j]
        end
        for j in 1:R
            cache.Yi[i] .+= params.tab.q.α[i,j] .* cache.Up[j]
        end

        # compute Q and P
        cache.Qi[i] .= params.q .+ params.Δt .* cache.Yi[i]
        cache.Pi[i] .= params.p .+ params.Δt .* cache.Zi[i]

        # compute f(X)
        tpᵢ = params.t + params.Δt * params.tab.p.c[i]
        params.f_f(tpᵢ, cache.Qi[i], cache.Vi[i], cache.Fi[i])
        params.f_p(tpᵢ, cache.Qi[i], cache.Vi[i], cache.Φi[i])

        cache.Φi[i] .-= cache.Pi[i]
    end

    for i in 1:R
        # compute Y
        cache.Yp[i] .= 0
        for j in 1:S
            cache.Yp[i] .+= params.tab.q̃.a[i,j] .* cache.Vi[j]
        end
        for j in 1:R
            cache.Yp[i] .+= params.tab.q̃.α[i,j] .* cache.Up[j]
        end

        # compute Q and P
        cache.Qp[i] .= params.q .+ params.Δt .* cache.Yp[i]
        cache.Pp[i] .= params.p .+ params.Δt .* cache.Zp[i]

        # compute f(X)
        tλᵢ = params.t + params.Δt * params.tab.λ.c[i]
        params.f_g(tλᵢ, cache.Qp[i], cache.Pp[i], cache.Λp[i], cache.Gp[i])
        params.f_ϕ(tλᵢ, cache.Qp[i], cache.Pp[i], cache.Φp[i])
    end

    if length(params.tab.d) > 0
        for k in 1:D
            cache.μ[k] = x[2*D*S+2*D*R+k]
        end
    end

    # compute q and p
    cache.q̃ .= params.q
    cache.p̃ .= params.p
    for i in 1:S
        cache.q̃ .+= params.Δt .* params.tab.q.b[i] .* cache.Vi[i]
        cache.p̃ .+= params.Δt .* params.tab.p.b[i] .* cache.Fi[i]
    end
    for i in 1:R
        cache.q̃ .+= params.Δt .* params.tab.q.β[i] .* cache.Up[i]
        cache.p̃ .+= params.Δt .* params.tab.p.β[i] .* cache.Gp[i]
    end

    # compute ϕ(q,p)
    tλᵢ = params.t + params.Δt
    params.f_ϕ(tλᵢ, cache.q̃, cache.p̃, cache.ϕ̃)
end


"Compute stages of specialised partitioned additive Runge-Kutta methods for variational systems."
@generated function function_stages!(y::Vector{ST}, b::Vector{ST}, params::ParametersVSPARKprimary{DT,TT,D,S,R,P}) where {ST,DT,TT,D,S,R,P}
    cache = IntegratorCacheSPARK{ST,TT,D,S,R}()

    quote
        compute_stages!(y, $cache, params)

        # compute b = - [Φ, (Z-AF-AG)]
        for i in 1:S
            for k in 1:D
                b[2*(D*(i-1)+k-1)+1] = - $cache.Φi[i][k]
                b[2*(D*(i-1)+k-1)+2] = - $cache.Zi[i][k]
                for j in 1:S
                    b[2*(D*(i-1)+k-1)+2] += params.tab.p.a[i,j] * $cache.Fi[j][k]
                end
                for j in 1:R
                    b[2*(D*(i-1)+k-1)+2] += params.tab.p.α[i,j] * $cache.Gp[j][k]
                end
            end
        end

        # compute b = - [ωΦ, (Z-AF-AG)]
        for i in 1:R
            for k in 1:D
                b[2*D*S+2*(D*(i-1)+k-1)+2] = - $cache.Zp[i][k]
                for j in 1:S
                    b[2*D*S+2*(D*(i-1)+k-1)+2] += params.tab.p̃.a[i,j] * $cache.Fi[j][k]
                end
                for j in 1:R
                    b[2*D*S+2*(D*(i-1)+k-1)+2] += params.tab.p̃.α[i,j] * $cache.Gp[j][k]
                end
            end
        end
        for i in 1:R-P
            for k in 1:D
                b[2*D*S+2*(D*(i-1)+k-1)+1] = 0
                for j in 1:R
                    b[2*D*S+2*(D*(i-1)+k-1)+1] -= params.tab.ω[i,j] * $cache.Φp[j][k]
                end
                b[2*D*S+2*(D*(i-1)+k-1)+1] -= params.tab.ω[i,R+1] * $cache.ϕ̃[k]
            end
        end

        # compute b = d_λ ⋅ Λ
        for i in R-P+1:R
            for k in 1:D
                b[2*D*S+2*(D*(R-1)+k-1)+1] = 0
                for j in 1:R
                    b[2*D*S+2*(D*(i-1)+k-1)+1] -= params.tab.δ[j] * $cache.Λp[j][k]
                end
            end
        end

        if length(params.tab.d) > 0
            for i in 1:S
                for k in 1:D
                    b[2*(D*(i-1)+k-1)+2] -= $cache.μ[k] * params.tab.d[i]
                end
            end

            for k in 1:D
                b[2*D*S+2*D*R+k] = 0
                for i in 1:S
                    b[2*D*S+2*D*R+k] -= $cache.Vi[i][k] * params.tab.d[i]
                end
            end
        end
    end
end


function initial_guess!(int::IntegratorVSPARKprimary, sol::AtomicSolutionPDAE)
    for i in eachstage(int)
        evaluate!(int.iguess, sol.q, sol.p, sol.v, sol.f,
                              sol.q̅, sol.p̅, sol.v̅, sol.f̅,
                              int.cache.q̃, int.cache.p̃, int.cache.ṽ, int.cache.f̃,
                              tableau(int).q.c[i], tableau(int).p.c[i])

        for k in eachdim(int)
            int.solver.x[2*(ndims(int)*(i-1)+k-1)+1] = int.cache.ṽ[k]
            int.solver.x[2*(ndims(int)*(i-1)+k-1)+2] = (int.cache.p̃[k] - sol.p[k])/timestep(int)
        end
    end

    for i in 1:pstages(int)
        evaluate!(int.iguess, sol.q, sol.p, sol.v, sol.f,
                              sol.q̅, sol.p̅, sol.v̅, sol.f̅,
                              int.cache.q̃, int.cache.p̃, int.cache.ṽ, int.cache.f̃,
                              tableau(int).q̃.c[i], tableau(int).p̃.c[i])

        for k in eachdim(int)
            int.solver.x[2*ndims(int)*nstages(int)+2*(ndims(int)*(i-1)+k-1)+1] = 0
            int.solver.x[2*ndims(int)*nstages(int)+2*(ndims(int)*(i-1)+k-1)+2] = (int.cache.p̃[k] - sol.p[k])/timestep(int)
        end
    end

    if isdefined(tableau(int), :λ) && tableau(int).λ.c[1] == 0
        for k in eachdim(int)
            int.solver.x[2*ndims(int)*nstages(int)+2*(k-1)+1] = int.cache.λ[k]
        end
    end

    if isdefined(tableau(int), :d) && length(tableau(int).d) > 0
        for k in eachdim(int)
            int.solver.x[2*ndims(int)*nstages(int)+2*ndims(int)*pstages(int)+k] = 0
        end
    end
end
