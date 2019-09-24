
const TableauVSPARK{DT} = AbstractTableauVSPARK{:vspark,DT}
const ParametersVSPARK{DT,TT,D,S,R,FT,PT,UT,GT,ϕT} = AbstractParametersVSPARK{:vspark,DT,TT,D,S,R,FT,PT,UT,GT,ϕT}


@doc raw"""
Specialised Partitioned Additive Runge-Kutta integrator for Variational systems.

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
p_{n+1} &= p_{n} + h \sum \limits_{i=1}^{s} b_{i} F_{n,i} + h \sum \limits_{i=1}^{r} \beta_{i} G_{n,i} , \\
0 &= \phi (q_{n+1}, p_{n+1}) .
\end{align}
```
"""
struct IntegratorVSPARK{DT, TT, ET <: IDAE{DT,TT},
                                PT <: ParametersVSPARK{DT,TT},
                                ST <: NonlinearSolver{DT},
                                IT <: InitialGuessPODE{DT,TT}} <: AbstractIntegratorVSPARK{DT, TT}
    equation::ET
    tableau::TableauVSPARK{TT}

    params::PT
    solver::ST
    iguess::IT
end

function IntegratorVSPARK(equation::IDAE{DT,TT,FT,PT,UT,GT,ϕT,VT},
                         tableau::TableauVSPARK{TT}, Δt::TT) where {DT,TT,FT,PT,UT,GT,ϕT,VT}
    D = equation.d
    S = tableau.s
    R = tableau.r

    @assert tableau.ρ == tableau.r-1

    N = 3*D*S + 3*D*R

    if isdefined(tableau, :d)
        N += D
        d_v = tableau.d
    else
        d_v = DT[]
    end

    # create params
    params = ParametersVSPARK{DT,TT,D,S,R,FT,PT,UT,GT,ϕT}(
                                equation.f, equation.p, equation.u, equation.g, equation.ϕ, Δt,
                                tableau.q, tableau.p, tableau.q̃, tableau.p̃, tableau.λ, tableau.ωt_δ, d_v)

    # create solver
    solver = create_nonlinear_solver(DT, N, params)

    # create initial guess
    iguess = InitialGuessPODE(get_config(:ig_interpolation), equation, Δt)

    # create integrator
    IntegratorVSPARK{DT, TT, typeof(equation), typeof(params), typeof(solver), typeof(iguess)}(
                                        equation, tableau, params, solver, iguess)
end


function compute_stages!(x::Vector{ST}, cache::IntegratorCacheVSPARK{ST,TT,D,S,R},
                                        params::ParametersVSPARK{DT,TT,D,S,R}) where {ST,DT,TT,D,S,R}
    local tpᵢ::TT
    local tλᵢ::TT

    for i in 1:S
        for k in 1:D
            # copy x to Y, Z
            cache.Yi[i][k] = x[3*(D*(i-1)+k-1)+1]
            cache.Zi[i][k] = x[3*(D*(i-1)+k-1)+2]
            cache.Vi[i][k] = x[3*(D*(i-1)+k-1)+3]

            # compute Q and P
            cache.Qi[i][k] = params.q[k] + params.Δt * cache.Yi[i][k]
            cache.Pi[i][k] = params.p[k] + params.Δt * cache.Zi[i][k]
        end

        # compute f(X)
        tpᵢ = params.t + params.Δt * params.t_p.c[i]
        params.f_f(tpᵢ, cache.Qi[i], cache.Vi[i], cache.Fi[i])
        params.f_p(tpᵢ, cache.Qi[i], cache.Vi[i], cache.Φi[i])
        cache.Φi[i] .-= cache.Pi[i]
    end

    for i in 1:R
        for k in 1:D
            # copy y to Y, Z and Λ
            cache.Yp[i][k] = x[3*D*S+3*(D*(i-1)+k-1)+1]
            cache.Zp[i][k] = x[3*D*S+3*(D*(i-1)+k-1)+2]
            cache.Λp[i][k] = x[3*D*S+3*(D*(i-1)+k-1)+3]

            # compute Q and V
            cache.Qp[i][k] = params.q[k] + params.Δt * cache.Yp[i][k]
            cache.Pp[i][k] = params.p[k] + params.Δt * cache.Zp[i][k]
        end

        # compute f(X)
        tλᵢ = params.t + params.Δt * params.t_λ.c[i]
        params.f_u(tλᵢ, cache.Qp[i], cache.Pp[i], cache.Λp[i], cache.Up[i])
        params.f_g(tλᵢ, cache.Qp[i], cache.Pp[i], cache.Λp[i], cache.Gp[i])
        params.f_ϕ(tλᵢ, cache.Qp[i], cache.Pp[i], cache.Φp[i])
    end

    if length(params.d_v) > 0
        for k in 1:D
            cache.μ[k] = x[3*D*S+3*D*R+k]
        end
    end

    # compute q and p
    cache.q̃ .= params.q
    cache.p̃ .= params.p
    for i in 1:S
        cache.q̃ .+= params.Δt .* params.t_q.b[i] .* cache.Vi[i]
        cache.p̃ .+= params.Δt .* params.t_p.b[i] .* cache.Fi[i]
    end
    for i in 1:R
        cache.q̃ .+= params.Δt .* params.t_q.β[i] .* cache.Up[i]
        cache.p̃ .+= params.Δt .* params.t_p.β[i] .* cache.Gp[i]
    end

    # compute ϕ(q,p)
    tλᵢ = params.t + params.Δt
    params.f_ϕ(tλᵢ, cache.q̃, cache.p̃, cache.ϕ̃)
end


"Compute stages of specialised partitioned additive Runge-Kutta methods for variational systems."
@generated function function_stages!(y::Vector{ST}, b::Vector{ST}, params::ParametersVSPARK{DT,TT,D,S,R}) where {ST,DT,TT,D,S,R}
    cache = IntegratorCacheVSPARK{ST,TT,D,S,R}()

    quote
        compute_stages!(y, $cache, params)

        # compute b = - [(Y-AV-AU), (Z-AF-AG), Φ]
        for i in 1:S
            for k in 1:D
                b[3*(D*(i-1)+k-1)+1] = - $cache.Yi[i][k]
                b[3*(D*(i-1)+k-1)+2] = - $cache.Zi[i][k]
                b[3*(D*(i-1)+k-1)+3] = - $cache.Φi[i][k]
                for j in 1:S
                    b[3*(D*(i-1)+k-1)+1] += params.t_q.a[i,j] * $cache.Vi[j][k]
                    b[3*(D*(i-1)+k-1)+2] += params.t_p.a[i,j] * $cache.Fi[j][k]
                end
                for j in 1:R
                    b[3*(D*(i-1)+k-1)+1] += params.t_q.α[i,j] * $cache.Up[j][k]
                    b[3*(D*(i-1)+k-1)+2] += params.t_p.α[i,j] * $cache.Gp[j][k]
                end
            end
        end

        # compute b = - [(Y-AV-AU), (Z-AF-AG), ωΦ]
        for i in 1:R
            for k in 1:D
                b[3*D*S+3*(D*(i-1)+k-1)+1] = - $cache.Yp[i][k]
                b[3*D*S+3*(D*(i-1)+k-1)+2] = - $cache.Zp[i][k]
                for j in 1:S
                    b[3*D*S+3*(D*(i-1)+k-1)+1] += params.t_q̃.a[i,j] * $cache.Vi[j][k]
                    b[3*D*S+3*(D*(i-1)+k-1)+2] += params.t_p̃.a[i,j] * $cache.Fi[j][k]
                end
                for j in 1:R
                    b[3*D*S+3*(D*(i-1)+k-1)+1] += params.t_q̃.α[i,j] * $cache.Up[j][k]
                    b[3*D*S+3*(D*(i-1)+k-1)+2] += params.t_p̃.α[i,j] * $cache.Gp[j][k]
                end
            end
        end
        for i in 1:R-1
            for k in 1:D
                b[3*D*S+3*(D*(i-1)+k-1)+3] = 0
                for j in 1:R
                    b[3*D*S+3*(D*(i-1)+k-1)+3] -= params.t_ω[i,j] * $cache.Φp[j][k]
                end
            end
        end

        # compute b = -ϕ
        for k in 1:D
            b[3*D*S+3*(D*(R-1)+k-1)+3] = - $cache.ϕ̃[k]
        end

        if length(params.d_v) > 0
            for i in 1:S
                for k in 1:D
                    b[3*(D*(i-1)+k-1)+3] -= $cache.μ[k] * params.d_v[i]
                end
            end

            for k in 1:D
                b[3*D*S+3*D*R+k] = 0
                for i in 1:S
                    b[3*D*S+3*D*R+k] -= $cache.Vi[i][k] * params.d_v[i]
                end
            end
        end
    end
end
