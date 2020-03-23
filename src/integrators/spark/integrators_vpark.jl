
"Holds the tableau of an Variational Partitioned Additive Runge-Kutta method."
const TableauVPARK = AbstractTableauSPARK{:vpark}

"Parameters for right-hand side function of Variational Partitioned Additive Runge-Kutta methods."
const ParametersVPARK = AbstractParametersSPARK{:vpark}



@doc raw"""
Variational partitioned additive Runge-Kutta integrator.

This integrator solves the following system of equations for the internal stages,
```math
\begin{aligned}
Q_{n,i} &= q_{n} + h \sum \limits_{j=1}^{s} a_{ij} V_{n,j} + h \sum \limits_{j=1}^{r} \alpha_{ij} U_{n,j} , & i &= 1, ..., s , \\
P_{n,i} &= p_{n} + h \sum \limits_{j=1}^{s} a_{ij} F_{n,j} + h \sum \limits_{j=1}^{r} \alpha_{ij} G_{n,j} , & i &= 1, ..., s , \\
\tilde{Q}_{n,i} &= q_{n} + h \sum \limits_{j=1}^{s} \tilde{a}_{ij} V_{n,j} + h \sum \limits_{j=1}^{r} \tilde{\alpha}_{ij} U_{n,j} , & i &= 1, ..., r , \\
\tilde{P}_{n,i} &= p_{n} + h \sum \limits_{j=1}^{s} \tilde{a}_{ij} F_{n,j} + h \sum \limits_{j=1}^{r} \tilde{\alpha}_{ij} G_{n,j} , & i &= 1, ..., r , \\
\tilde{\Phi}_{n,i} &= 0 , & i &= 1, ..., r ,
\end{aligned}
```
with definitions
```math
\begin{aligned}
P_{n,i} &= \frac{\partial L}{\partial v} (Q_{n,i}, V_{n,i}) , & i &= 1, ..., s , \\
F_{n,i} &= \frac{\partial L}{\partial q} (Q_{n,i}, V_{n,i}) , & i &= 1, ..., s , \\
U_{n,i} &= \hphantom{-} \frac{\partial \phi}{\partial p} (\tilde{Q}_{n,i}, \tilde{P}_{n,i})^{T} \Lambda_{n,i} , & i &= 1, ..., r , \\
G_{n,i} &=           -  \frac{\partial \phi}{\partial q} (\tilde{Q}_{n,i}, \tilde{P}_{n,i})^{T} \Lambda_{n,i} , & i &= 1, ..., r , \\
\tilde{\Phi}_{n,i} &= \phi(\tilde{Q}_{n,i}, \tilde{P}_{n,i}) , & i &= 1, ..., r ,
\end{aligned}
```
and update rule
```math
\begin{aligned}
q_{n+1} &= q_{n} + h \sum \limits_{i=1}^{s} b_{i} V_{n,i} + h \sum \limits_{i=1}^{r} \beta_{i} U_{n,i} , \\
p_{n+1} &= p_{n} + h \sum \limits_{i=1}^{s} b_{i} F_{n,i} + h \sum \limits_{i=1}^{r} \beta_{i} G_{n,i} .
\end{aligned}
```
"""
struct IntegratorVPARK{DT, TT, D, S, R, PT <: ParametersVPARK{DT,TT,D,S,R},
                                        ST <: NonlinearSolver{DT},
                                        IT <: InitialGuessPODE{DT,TT}} <: AbstractIntegratorVSPARK{DT,TT,D,S,R}
    params::PT
    solver::ST
    iguess::IT
    caches::CacheDict{PT}

    function IntegratorVPARK(params::ParametersVPARK{DT,TT,D,S,R}, solver::ST, iguess::IT, caches) where {DT,TT,D,S,R,ST,IT}
        new{DT, TT, D, S, R, typeof(params), ST, IT}(params, solver, iguess, caches)
    end

    function IntegratorVPARK{DT,D}(equations::NamedTuple, tableau::TableauVPARK{TT}, Δt::TT) where {DT,TT,D}
        # get number of stages
        S = tableau.s
        R = tableau.r

        N = 3*D*S + 3*D*R

        if isdefined(tableau, :d) && length(tableau.d) > 0
            N += D
        end

        # create params
        params = ParametersVPARK{DT,D}(equations, tableau, Δt)

        # create cache dict
        caches = CacheDict(params)

        # create solver
        solver = create_nonlinear_solver(DT, N, params, caches)

        # create initial guess
        iguess = InitialGuessPODE{DT,D}(get_config(:ig_interpolation), equations[:v], equations[:f], Δt)

        # create integrator
        IntegratorVPARK(params, solver, iguess, caches)
    end

    function IntegratorVPARK(equation::IDAE{DT,TT}, tableau::TableauVPARK{TT}, Δt::TT; kwargs...) where {DT,TT}
        IntegratorVPARK{DT, equation.d}(get_function_tuple(equation), tableau, Δt; kwargs...)
    end
end


function compute_stages!(x::Vector{ST}, cache::IntegratorCacheSPARK{ST,D,S,R},
                                        params::ParametersVPARK{DT,TT,D,S,R}) where {ST,DT,TT,D,S,R}
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
        tpᵢ = params.t + params.Δt * params.tab.p.c[i]
        params.equs[:ϑ](tpᵢ, cache.Qi[i], cache.Vi[i], cache.Φi[i])
        params.equs[:f](tpᵢ, cache.Qi[i], cache.Vi[i], cache.Fi[i])

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
        tλᵢ = params.t + params.Δt * params.tab.λ.c[i]
        params.equs[:u](tλᵢ, cache.Qp[i], cache.Pp[i], cache.Λp[i], cache.Up[i])
        params.equs[:g](tλᵢ, cache.Qp[i], cache.Pp[i], cache.Λp[i], cache.Gp[i])
        params.equs[:ϕ](tλᵢ, cache.Qp[i], cache.Pp[i], cache.Φp[i])
    end

    if isdefined(params.tab, :d) && length(params.tab.d) > 0
        for k in 1:D
            cache.μ[k] = x[3*D*S+3*D*R+k]
        end
    end
end


"Compute stages of variational partitioned additive Runge-Kutta methods."
function Integrators.function_stages!(y::Vector{ST}, b::Vector{ST}, params::ParametersVPARK{DT,TT,D,S,R},
                                      caches::CacheDict) where {ST,DT,TT,D,S,R}

    # get cache for internal stages
    cache = caches[ST]

    compute_stages!(y, cache, params)

    # compute b = - [(Y-AV-AU), (Z-AF-AG), Φ]
    for i in 1:S
        for k in 1:D
            b[3*(D*(i-1)+k-1)+1] = - cache.Yi[i][k]
            b[3*(D*(i-1)+k-1)+2] = - cache.Zi[i][k]
            b[3*(D*(i-1)+k-1)+3] = - cache.Φi[i][k]
            for j in 1:S
                b[3*(D*(i-1)+k-1)+1] += params.tab.q.a[i,j] * cache.Vi[j][k]
                b[3*(D*(i-1)+k-1)+2] += params.tab.p.a[i,j] * cache.Fi[j][k]
            end
            for j in 1:R
                b[3*(D*(i-1)+k-1)+1] += params.tab.q.α[i,j] * cache.Up[j][k]
                b[3*(D*(i-1)+k-1)+2] += params.tab.p.α[i,j] * cache.Gp[j][k]
            end
        end
    end

    # compute b = - [(Y-AV-AU), (Z-AF-AG), Φ]
    for i in 1:R
        for k in 1:D
            b[3*D*S+3*(D*(i-1)+k-1)+1] = - cache.Yp[i][k]
            b[3*D*S+3*(D*(i-1)+k-1)+2] = - cache.Zp[i][k]
            b[3*D*S+3*(D*(i-1)+k-1)+3] = - cache.Φp[i][k]
            for j in 1:S
                b[3*D*S+3*(D*(i-1)+k-1)+1] += params.tab.q̃.a[i,j] * cache.Vi[j][k]
                b[3*D*S+3*(D*(i-1)+k-1)+2] += params.tab.p̃.a[i,j] * cache.Fi[j][k]
            end
            for j in 1:R
                b[3*D*S+3*(D*(i-1)+k-1)+1] += params.tab.q̃.α[i,j] * cache.Up[j][k]
                b[3*D*S+3*(D*(i-1)+k-1)+2] += params.tab.p̃.α[i,j] * cache.Gp[j][k]
            end
        end
    end

    # compute b = - [Λ₁-λ]
    if params.tab.λ.c[1] == 0
        for k in 1:D
            b[3*D*S+3*(k-1)+3] = - cache.Λp[1][k] + params.λ[k]
        end
    end

    if isdefined(params.tab, :d) && length(params.tab.d) > 0
        for i in 1:S
            for k in 1:D
                b[3*(D*(i-1)+k-1)+3] -= cache.μ[k] * params.tab.d[i]
            end
        end

        for k in 1:D
            b[3*D*S+3*D*R+k] = 0
            for i in 1:S
                b[3*D*S+3*D*R+k] -= cache.Vi[i][k] * params.tab.d[i]
            end
        end
    end
end


function update_solution!(int::IntegratorVPARK{DT,TT}, sol::AtomicSolutionPDAE{DT,TT},
                          cache::IntegratorCacheSPARK{DT}=int.caches[DT]) where {DT,TT}
    # compute final update
    update_solution!(sol.q, sol.q̃, cache.Vi, int.params.tab.q.b, timestep(int))
    update_solution!(sol.p, sol.p̃, cache.Fi, int.params.tab.p.b, timestep(int))

    # compute projection
    update_solution!(sol.q, sol.q̃, cache.Up, int.params.tab.q.β, timestep(int))
    update_solution!(sol.p, sol.p̃, cache.Gp, int.params.tab.p.β, timestep(int))
    update_multiplier!(sol.λ, cache.Λp, int.params.tab.λ.b)
end
