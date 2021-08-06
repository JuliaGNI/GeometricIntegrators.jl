"Holds the tableau of an Specialised Partitioned Additive Runge-Kutta method for Variational systems."
struct TableauVSPARKsecondary{DT <: Number} <: AbstractTableau{DT}
    name::Symbol
    o::Int
    s::Int
    r::Int
    ρ::Int

    q::CoefficientsSPARK{DT}
    p::CoefficientsSPARK{DT}

    q̃::CoefficientsSPARK{DT}
    p̃::CoefficientsSPARK{DT}

    ω::Matrix{DT}
    d::Vector{DT}

    function TableauVSPARKsecondary(name::Symbol, o::Int, s::Int, r::Int,
                        q::CoefficientsSPARK{DT}, p::CoefficientsSPARK{DT},
                        q̃::CoefficientsSPARK{DT}, p̃::CoefficientsSPARK{DT},
                        ω::Matrix{DT}, d=DT[]) where {DT}

        @assert s > 0 "Number of stages s must be > 0"
        @assert r > 0 "Number of stages r must be > 0"

        @assert s==q.s==p.s==q̃.σ==p̃.σ
        @assert r==q.σ==p.σ==q̃.s==p̃.s
        @assert size(ω,1)==r
        @assert size(ω,2)==r+1

        @assert length(d)==0 || length(d)==r

        new{DT}(name, o, s, r, 0, q, p, q̃, p̃, ω, d)
    end
end


"Parameters for right-hand side function of Specialised Partitioned Additive Runge-Kutta methods for Variational systems."
const ParametersVSPARKsecondary = AbstractParametersSPARK{:vspark_secondary}


@doc raw"""
Specialised Partitioned Additive Runge-Kutta integrator for degenerate
variational systems with projection on secondary constraint.

This integrator solves the following system of equations for the internal stages,
```math
\begin{aligned}
Q_{n,i} &= q_{n} + h \sum \limits_{j=1}^{s} a^1_{ij} V_{n,j} + h \sum \limits_{j=1}^{\sigma} a^2_{ij} \tilde{\Lambda}_{n,j} , & i &= 1, ..., s , \\
P_{n,i} &= p_{n} + h \sum \limits_{j=1}^{s} \bar{a}^1_{ij} F^1_{n,j} + h \sum \limits_{j=1}^{s} \bar{a}^2_{ij} F^2_{n,j} + h \sum \limits_{j=1}^{\sigma} \bar{a}^3_{ij} \tilde{F}^3_{n,j} , & i &= 1, ..., s , \\
0 &= \Phi_{n,i} , & i &= 1, ..., s ,
\end{aligned}
```
the projective stages
```math
\begin{aligned}
\tilde{Q}_{n,i} &= q_{n} + h \sum \limits_{j=1}^{s} \alpha^1_{ij} \tilde{V}_{n,j} + h \sum \limits_{j=1}^{\sigma} \alpha^2_{ij} \tilde{\Lambda}_{n,j} , & i &= 1, ..., \sigma , \\
\tilde{P}_{n,i} &= p_{n} + h \sum \limits_{j=1}^{s} \bar{\alpha}^1_{ij} F^1_{n,j} + h \sum \limits_{j=1}^{s} \bar{\alpha}^2_{ij} \tilde{F}^2_{n,j} + h \sum \limits_{j=1}^{\sigma} \bar{\alpha}^3_{ij} \tilde{F}^3_{n,j} , & i &= 1, ..., \sigma , \\
0 &= \tilde{\Phi}_{n,i} , & i &= 1, ..., \sigma , \\
0 &= \sum \limits_{j=1}^{\sigma} \omega_{ij} \tilde{\Psi}_{n,i} , & i &= 1, ..., \sigma-1 ,
\end{aligned}
```
and update rule
```math
\begin{aligned}
q_{n+1} &= q_{n} + h \sum \limits_{i=1}^{s} b^1_{i} V_{n,i} + h \sum \limits_{i=1}^{\sigma} b^2_{i} \tilde{\Lambda}_{n,i} , \\
p_{n+1} &= p_{n} + h \sum \limits_{i=1}^{s} b^1_{i} F^1_{n,i} + h \sum \limits_{i=1}^{s} b^2_{i} F^2_{n,i} + h \sum \limits_{i=1}^{\sigma} b^3_{i} \tilde{F}^3_{n,i} , \\
0 &= \phi (q_{n+1}, p_{n+1}) ,
\end{aligned}
```
with definitions
```math
\begin{aligned}
F^1_{n,i} &= \nabla H (Q_{n,i}) , & i &= 1, ..., s , \\
F^2_{n,i} &= \nabla \vartheta (Q_{n,i}) \cdot V_{n,i} , & i &= 1, ..., s , \\
\tilde{F}^3_{n,i} &= \nabla \phi (\tilde{Q}_{n,i}, \tilde{P}_{n,i}) \cdot \tilde{\Lambda}_{n,i} , & i &= 1, ..., \sigma , \\
\Phi_{n,i} &= \phi (Q_{n,i}, P_{n,i}) = P_{n,i} - \vartheta (Q_{n,i}) , & i &= 1, ..., s , \\
\tilde{\Phi}_{n,i} &= \phi (\tilde{Q}_{n,i}, \tilde{P}_{n,i}) = \tilde{P}_{n,i} - \vartheta (\tilde{Q}_{n,i}) , & i &= 1, ..., \sigma , \\
\tilde{\Psi}_{n,i} &= \psi (\tilde{Q}_{n,i}, \tilde{V}_{n,i}, \tilde{P}_{n,i}, \tilde{F}_{n,i}) = \tilde{F}_{n,i} - \tilde{V}_{n,i} \cdot \nabla \vartheta (\tilde{Q}_{n,i}) , & i &= 1, ..., \sigma ,
\end{aligned}
```
so that
```math
\begin{aligned}
F^1_{n,i} + F^2_{n,i} &= \frac{\partial L}{\partial q} (Q_{n,i}, V_{n,i}) , & i &= 1, ..., s , \\
\phi (Q_{n,i}, P_{n,i}) &= P_{n,i} - \frac{\partial L}{\partial v} (Q_{n,i}, V_{n,i}) , & i &= 1, ..., s , \\
\psi (\tilde{Q}_{n,i}, \tilde{V}_{n,i}, \tilde{P}_{n,i}, \tilde{F}_{n,i}) &= \tilde{F}_{n,i} - \tilde{V}_{n,i} \cdot \nabla \frac{\partial L}{\partial v} (\tilde{Q}_{n,i}, \tilde{V}_{n,i}) , & i &= 1, ..., \sigma .
\end{aligned}
```
"""
struct IntegratorVSPARKsecondary{DT, TT, D, S, R,
                                         PT <: ParametersVSPARKsecondary{DT,TT,D,S,R},
                                         ST <: NonlinearSolver{DT},
                                         IT <: InitialGuessIODE{TT}} <: AbstractIntegratorVSPARK{DT,TT,D,S,R}
    params::PT
    solver::ST
    iguess::IT
    caches::CacheDict{PT}

    function IntegratorVSPARKsecondary(params::ParametersVSPARKsecondary{DT,TT,D,S,R}, solver::ST, iguess::IT, caches) where {DT,TT,D,S,R,ST,IT}
        new{DT, TT, D, S, R, typeof(params), ST, IT}(params, solver, iguess, caches)
    end

    function IntegratorVSPARKsecondary{DT,D}(equations::NamedTuple, tableau::TableauVSPARKsecondary{TT}, Δt::TT) where {DT,TT,D,ST}
        # get number of stages
        S = tableau.s
        R = tableau.r

        N = 4*D*R

        if isdefined(tableau, :d) && length(tableau.d) > 0
            N += D
        end

        # create params
        params = ParametersVSPARKsecondary{DT,D}(equations, tableau, Δt)

        # create cache dict
        caches = CacheDict(params)

        # create solver
        solver = create_nonlinear_solver(DT, N, params, caches)

        # create initial guess
        iguess = InitialGuessIODE(get_config(:ig_extrapolation), equations[:v̄], equations[:f̄], Δt)

        # create integrator
        IntegratorVSPARKsecondary(params, solver, iguess, caches)
    end

    function IntegratorVSPARKsecondary(equation::LDAE{DT}, tableau::TableauVSPARKsecondary, Δt; kwargs...) where {DT}
        IntegratorVSPARKsecondary{DT, ndims(equation)}(get_functions(equation), tableau, Δt; kwargs...)
    end
end


GeometricBase.nconstraints(::IntegratorVSPARKsecondary{DT,TT,D}) where {DT,TT,D} = D


function initial_guess!(int::IntegratorVSPARKsecondary{DT}, sol::AtomicSolutionPDAE{DT},
                        cache::IntegratorCacheSPARK{DT}=int.caches[DT]) where {DT}
    for i in 1:pstages(int)
        evaluate!(int.iguess, sol.q̄, sol.p̄, sol.v̄, sol.f̄,
                              sol.q, sol.p, sol.v, sol.f,
                              cache.q̃, cache.p̃, cache.ṽ, cache.f̃,
                              tableau(int).q̃.c[i], tableau(int).p̃.c[i])

        for k in 1:ndims(int)
            int.solver.x[4*(ndims(int)*(i-1)+k-1)+1] = (cache.q̃[k] - sol.q[k])/timestep(int)
            int.solver.x[4*(ndims(int)*(i-1)+k-1)+2] = (cache.p̃[k] - sol.p[k])/timestep(int)
            int.solver.x[4*(ndims(int)*(i-1)+k-1)+3] = cache.ṽ[k]
            int.solver.x[4*(ndims(int)*(i-1)+k-1)+4] = 0
        end
    end

    if isdefined(tableau(int), :d) && length(tableau(int).d) > 0
        for k in 1:ndims(int)
            int.solver.x[4*ndims(int)*pstages(int)+k] = 0
        end
    end
end


function compute_stages!(x::Vector{ST}, cache::IntegratorCacheSPARK{ST,D,S,R},
                                        params::ParametersVSPARKsecondary{DT,TT,D,S,R}) where {ST,DT,TT,D,S,R}
    local t::TT

    for i in 1:R
        for k in 1:D
            # copy y to Y, Z and Λ
            cache.Yp[i][k] = x[4*(D*(i-1)+k-1)+1]
            cache.Zp[i][k] = x[4*(D*(i-1)+k-1)+2]
            cache.Vp[i][k] = x[4*(D*(i-1)+k-1)+3]
            cache.Λp[i][k] = x[4*(D*(i-1)+k-1)+4]

            # compute Q and P
            cache.Qp[i][k] = params.q[k] + params.Δt * cache.Yp[i][k]
            cache.Pp[i][k] = params.p[k] + params.Δt * cache.Zp[i][k]
        end

        # compute f(X)
        t = params.t + params.Δt * params.tab.p̃.c[i]
        params.equs[:f](t, cache.Qp[i], cache.Vp[i], cache.Fp[i])
        params.equs[:g](t, cache.Qp[i], cache.Vp[i], cache.Gp[i])
        params.equs[:g](t, cache.Qp[i], cache.Λp[i], cache.G̅p[i])

        cache.Hp[i] .= cache.Fp[i] .+ cache.Gp[i]

        params.equs[:ϕ](t, cache.Qp[i], cache.Pp[i], cache.Φp[i])
        params.equs[:ψ](t, cache.Qp[i], cache.Pp[i], cache.Vp[i], cache.Hp[i], cache.Ψp[i])
    end

    if isdefined(params.tab, :d) && length(params.tab.d) > 0
        for k in 1:D
            cache.μ[k] = x[4*D*R+k]
        end
    end

    for i in 1:S
        # compute Q
        for k in 1:D
            cache.Yi[i][k] = 0
            for j in 1:R
                cache.Yi[i][k] += params.tab.q.a[1][i,j] * cache.Vp[j][k]
                cache.Yi[i][k] += params.tab.q.a[2][i,j] * cache.Λp[j][k]
            end
            cache.Qi[i][k] = params.q[k] + params.Δt * cache.Yi[i][k]
        end

        # compute f(X)
        t = params.t + params.Δt * params.tab.p.c[i]
        params.equs[:f](t, cache.Qi[i], cache.Vi[i], cache.Fi[i])
    end

    # compute q and p
    cache.q̃ .= params.q
    cache.p̃ .= params.p
    for i in 1:S
        cache.p̃ .+= params.Δt .* params.tab.p.b[1][i] .* cache.Fi[i]
    end
    for i in 1:R
        cache.q̃ .+= params.Δt .* params.tab.q.b[1][i] .* cache.Vp[i]
        cache.q̃ .+= params.Δt .* params.tab.q.b[2][i] .* cache.Λp[i]
        cache.p̃ .+= params.Δt .* params.tab.p.b[2][i] .* cache.Gp[i]
        cache.p̃ .+= params.Δt .* params.tab.p.b[3][i] .* cache.G̅p[i]
    end

    # compute ϕ(q,p)
    t = params.t + params.Δt
    params.equs[:ϕ](t, cache.q̃, cache.p̃, cache.ϕ̃)
end


"Compute stages of specialised partitioned additive Runge-Kutta methods for variational systems."
function Integrators.function_stages!(y::Vector{ST}, b::Vector{ST}, params::ParametersVSPARKsecondary{DT,TT,D,S,R},
                                      caches::CacheDict) where {ST,DT,TT,D,S,R}

    # get cache for internal stages
    cache = caches[ST]

    compute_stages!(y, cache, params)

    # compute b = - [(Y-AV-AU), (Z-AF-AG), Φ, ωΨ]
    for i in 1:R
        for k in 1:D
            b[4*(D*(i-1)+k-1)+1] = - cache.Yp[i][k]
            b[4*(D*(i-1)+k-1)+2] = - cache.Zp[i][k]
            b[4*(D*(i-1)+k-1)+3] = - cache.Φp[i][k]
            b[4*(D*(i-1)+k-1)+4] = 0
            for j in 1:S
                b[4*(D*(i-1)+k-1)+2] += params.tab.p̃.a[1][i,j] * cache.Fi[j][k]
            end
            for j in 1:R
                b[4*(D*(i-1)+k-1)+1] += params.tab.q̃.a[1][i,j] * cache.Vp[j][k]
                b[4*(D*(i-1)+k-1)+1] += params.tab.q̃.a[2][i,j] * cache.Λp[j][k]
                b[4*(D*(i-1)+k-1)+2] += params.tab.p̃.a[2][i,j] * cache.Gp[j][k]
                b[4*(D*(i-1)+k-1)+2] += params.tab.p̃.a[3][i,j] * cache.G̅p[j][k]
            end
            for j in 1:R
                b[4*(D*(i-1)+k-1)+4] -= params.tab.ω[i,j] * cache.Ψp[j][k]
            end
            b[4*(D*(i-1)+k-1)+4] -= params.tab.ω[i,R+1] * cache.ϕ̃[k]
        end
    end

    if isdefined(params.tab, :d) && length(params.tab.d) > 0
        for i in 1:R
            for k in 1:D
                b[4*(D*(i-1)+k-1)+2] -= cache.μ[k] * params.tab.d[i] / params.tab.p.b[2][i]
            end
        end

        for k in 1:D
            b[4*D*R+k] = 0
            for i in 1:R
                b[4*D*R+k] -= cache.Vp[i][k] * params.tab.d[i]
            end
        end
    end
end


function update_solution!(int::IntegratorVSPARKsecondary{DT,TT}, sol::AtomicSolutionPDAE{DT,TT},
                          cache::IntegratorCacheSPARK{DT}=int.caches[DT]) where {DT,TT}
    # compute final update
    update_solution!(sol.p, sol.p̃, cache.Fi, int.params.tab.p.b[1], timestep(int))

    # compute projection
    update_solution!(sol.q, sol.q̃, cache.Vp, int.params.tab.q.b[1], timestep(int))
    update_solution!(sol.q, sol.q̃, cache.Λp, int.params.tab.q.b[2], timestep(int))
    update_solution!(sol.p, sol.p̃, cache.Gp, int.params.tab.p.b[2], timestep(int))
    update_solution!(sol.p, sol.p̃, cache.G̅p, int.params.tab.p.b[3], timestep(int))
end
