
"Holds the tableau of an Partitioned Additive Runge-Kutta method for Hamiltonian systems."
TableauHPARK = TableauVPARK


"Parameters for right-hand side function of Hamiltonian partitioned additive Runge-Kutta methods."
mutable struct ParametersHPARK{DT, TT, D, S, R, ET <: NamedTuple} <: AbstractParametersSPARK{DT,TT}
    equs::ET
    tab::TableauHPARK{TT}
    Δt::TT

    @ParametersSPARK

    function ParametersHPARK{DT,D}(equs::ET, tab::TableauVPARK{TT}, Δt::TT) where {DT,TT,D,S,R,ET <: NamedTuple}
        # create solution vectors
        q = zeros(DT,D)
        p = zeros(DT,D)
        λ = zeros(DT,D)

        new{DT,TT,D,tab.s,tab.r,ET}(equs, tab, Δt, zero(TT), q, p, λ)
    end
end


@doc raw"""
Partitioned Additive Runge-Kutta integrator for Hamiltonian systems subject
to Dirac constraints.

This integrator solves the following system of equations for the internal stages,
```math
\begin{align}
Q_{n,i} &= q_{n} + h \sum \limits_{j=1}^{s} a_{ij} V_{n,j} + h \sum \limits_{j=1}^{r} \alpha_{ij} U_{n,j} , & i &= 1, ..., s , \\
P_{n,i} &= p_{n} + h \sum \limits_{j=1}^{s} a_{ij} F_{n,j} + h \sum \limits_{j=1}^{r} \alpha_{ij} G_{n,j} , & i &= 1, ..., s , \\
\tilde{Q}_{n,i} &= q_{n} + h \sum \limits_{j=1}^{s} \tilde{a}_{ij} V_{n,j} + h \sum \limits_{j=1}^{r} \tilde{\alpha}_{ij} U_{n,j} , & i &= 1, ..., r , \\
\tilde{P}_{n,i} &= p_{n} + h \sum \limits_{j=1}^{s} \tilde{a}_{ij} F_{n,j} + h \sum \limits_{j=1}^{r} \tilde{\alpha}_{ij} G_{n,j} , & i &= 1, ..., r , \\
\tilde{\Phi}_{n,i} &= 0 , & i &= 1, ..., r ,
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
p_{n+1} &= p_{n} + h \sum \limits_{i=1}^{s} b_{i} F_{n,i} + h \sum \limits_{i=1}^{r} \beta_{i} G_{n,i} .
\end{align}
```
"""
struct IntegratorHPARK{DT, TT, D, S, R, PT <: ParametersHPARK{DT,TT},
                                        ST <: NonlinearSolver{DT},
                                        IT <: InitialGuessPODE{DT,TT}} <: AbstractIntegratorHSPARK{DT,TT}
    params::PT
    solver::ST
    iguess::IT
    cache::IntegratorCacheSPARK{DT, TT, D, S, R}

    function IntegratorHPARK(params::ParametersHPARK{DT,TT,D,S,R}, solver::ST, iguess::IT, cache) where {DT,TT,D,S,R,ST,IT}
        new{DT, TT, D, S, R, typeof(params), ST, IT}(params, solver, iguess, cache)
    end

    function IntegratorHPARK{DT,D}(equations::NamedTuple, tableau::TableauHPARK{TT}, Δt::TT) where {DT,TT,D}
        # get number of stages
        S = tableau.s
        R = tableau.r

        N = 2*D*S + 3*D*R

        # create params
        params = ParametersHPARK{DT,D}(equations, tableau, Δt)

        # create solver
        solver = create_nonlinear_solver(DT, N, params)

        # create initial guess
        iguess = InitialGuessPODE{DT,D}(get_config(:ig_interpolation), equations[:v], equations[:f], Δt)

        # create cache
        cache = IntegratorCacheSPARK{DT, TT, D, S, R}()

        # create integrator
        IntegratorHPARK(params, solver, iguess, cache)
    end

    function IntegratorHPARK(equation::PDAE{DT,TT}, tableau::TableauHPARK{TT}, Δt::TT; kwargs...) where {DT,TT}
        IntegratorHPARK{DT, equation.d}(get_function_tuple(equation), tableau, Δt; kwargs...)
    end
end


@inline equation(int::IntegratorHPARK, i::Symbol) = int.params.equs[i]
@inline equations(int::IntegratorHPARK) = int.params.equs
@inline tableau(int::IntegratorHPARK) = int.params.tab
@inline nstages(int::IntegratorHPARK{DT,TT,D,S,R}) where {DT,TT,D,S,R} = S
@inline pstages(int::IntegratorHPARK{DT,TT,D,S,R}) where {DT,TT,D,S,R} = R
@inline Base.ndims(int::IntegratorHPARK{DT,TT,D,S,R}) where {DT,TT,D,S,R} = D


function Integrators.initialize!(int::IntegratorHPARK, sol::AtomicSolutionPDAE)
    sol.t̅ = sol.t - timestep(int)

    equation(int, :v)(sol.t, sol.q, sol.p, sol.v)
    equation(int, :f)(sol.t, sol.q, sol.p, sol.f)

    initialize!(int.iguess, sol.t, sol.q, sol.p, sol.v, sol.f,
                            sol.t̅, sol.q̅, sol.p̅, sol.v̅, sol.f̅)
end


function compute_stages!(x::Vector{ST}, cache::IntegratorCacheSPARK{ST,TT,D,S,R},
                                        params::ParametersHPARK{DT,TT,D,S,R}) where {ST,DT,TT,D,S,R}
    local tqᵢ::TT
    local tpᵢ::TT
    local tλᵢ::TT

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
        tqᵢ = params.t + params.Δt * params.tab.q.c[i]
        tpᵢ = params.t + params.Δt * params.tab.p.c[i]
        params.equs[:v](tqᵢ, cache.Qi[i], cache.Pi[i], cache.Vi[i])
        params.equs[:f](tpᵢ, cache.Qi[i], cache.Pi[i], cache.Fi[i])
    end

    for i in 1:R
        for k in 1:D
            # copy y to Y, Z and Λ
            cache.Yp[i][k] = x[2*D*S+3*(D*(i-1)+k-1)+1]
            cache.Zp[i][k] = x[2*D*S+3*(D*(i-1)+k-1)+2]
            cache.Λp[i][k] = x[2*D*S+3*(D*(i-1)+k-1)+3]

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
end


"Compute stages of variational partitioned additive Runge-Kutta methods."
@generated function Integrators.function_stages!(y::Vector{ST}, b::Vector{ST}, params::ParametersHPARK{DT,TT,D,S,R}) where {ST,DT,TT,D,S,R}
    cache = IntegratorCacheSPARK{ST,TT,D,S,R}()

    quote
        compute_stages!(y, $cache, params)

        # compute b = - [(Y-AV-AU), (Z-AF-AG)]
        for i in 1:S
            for k in 1:D
                b[2*(D*(i-1)+k-1)+1] = - $cache.Yi[i][k]
                b[2*(D*(i-1)+k-1)+2] = - $cache.Zi[i][k]
                for j in 1:S
                    b[2*(D*(i-1)+k-1)+1] += params.tab.q.a[i,j] * $cache.Vi[j][k]
                    b[2*(D*(i-1)+k-1)+2] += params.tab.p.a[i,j] * $cache.Fi[j][k]
                end
                for j in 1:R
                    b[2*(D*(i-1)+k-1)+1] += params.tab.q.α[i,j] * $cache.Up[j][k]
                    b[2*(D*(i-1)+k-1)+2] += params.tab.p.α[i,j] * $cache.Gp[j][k]
                end
            end
        end

        # compute b = - [(Y-AV-AU), (Z-AF-AG), Φ]
        for i in 1:R
            for k in 1:D
                b[2*D*S+3*(D*(i-1)+k-1)+1] = - $cache.Yp[i][k]
                b[2*D*S+3*(D*(i-1)+k-1)+2] = - $cache.Zp[i][k]
                b[2*D*S+3*(D*(i-1)+k-1)+3] = - $cache.Φp[i][k]
                for j in 1:S
                    b[2*D*S+3*(D*(i-1)+k-1)+1] += params.tab.q̃.a[i,j] * $cache.Vi[j][k]
                    b[2*D*S+3*(D*(i-1)+k-1)+2] += params.tab.p̃.a[i,j] * $cache.Fi[j][k]
                end
                for j in 1:R
                    b[2*D*S+3*(D*(i-1)+k-1)+1] += params.tab.q̃.α[i,j] * $cache.Up[j][k]
                    b[2*D*S+3*(D*(i-1)+k-1)+2] += params.tab.p̃.α[i,j] * $cache.Gp[j][k]
                end
            end
        end

        # compute b = - [Λ₁-λ]
        if params.tab.λ.c[1] == 0
            for k in 1:D
                b[2*D*S+3*(k-1)+3] = - $cache.Λp[1][k] + params.λ[k]
            end
        end
    end
end


function update_solution!(int::IntegratorHPARK{DT,TT}, sol::AtomicSolutionPDAE{DT,TT}) where {DT,TT}
    # compute final update
    update_solution!(sol.q, sol.q̃, int.cache.Vi, int.params.tab.q.b, timestep(int))
    update_solution!(sol.p, sol.p̃, int.cache.Fi, int.params.tab.p.b, timestep(int))

    # compute projection
    update_solution!(sol.q, sol.q̃, int.cache.Up, int.params.tab.q.β, timestep(int))
    update_solution!(sol.p, sol.p̃, int.cache.Gp, int.params.tab.p.β, timestep(int))
    # update_multiplier!(sol.λ, int.cache.Λp, int.params.tab.λ.b)
end
