
"Holds the tableau of an Specialised Partitioned Additive Runge-Kutta method."
const TableauSPARK = AbstractTableauSPARK{:spark}
const ParametersSPARK = AbstractParametersSPARK{:spark}

function RungeKutta.check_symplecticity(tab::TableauSPARK{T}; atol=16*eps(T), rtol=16*eps(T)) where {T}
    s_αa_qp̃ = [isapprox(tab.q.b[i] * tab.p.α[i,j] + tab.p.β[j] * tab.q̃.a[j,i], tab.q.b[i] * tab.p.β[j]; atol=atol, rtol=rtol) for i in 1:tab.s, j in 1:tab.r]
    s_α_q̃p̃  = [isapprox(tab.p.β[i] * tab.q̃.α[i,j] + tab.q.β[j] * tab.p̃.α[j,i], tab.p.β[i] * tab.q.β[j]; atol=atol, rtol=rtol) for i in 1:tab.r, j in 1:tab.r]
    s_b_qp  = isapprox.(tab.q.b, tab.p.b; atol=atol, rtol=rtol)
    s_ω     = [isapprox(tab.ω[i,j], (i == j ? 1 : 0); atol=atol, rtol=rtol) for i in 1:tab.s, j in 1:tab.s+1]

    return (s_αa_qp̃, s_α_q̃p̃, s_b_qp, s_ω)
end

function Integrators.symplecticity_conditions(::TableauSPARK)
    (
        """`` b^{1}_{i} b^{3}_{j} = b^{3}_{j} \\tilde{a}^{1}_{ji} + b^{1}_{i} a^{3}_{ij} ``""",
        """`` b^{2}_{i} b^{3}_{j} = b^{3}_{j} \\tilde{a}^{2}_{ji} + b^{2}_{i} \\tilde{a}^{3}_{ij} ``""",
        """`` b^{3}_{i} = b^{2}_{i} ``""",
        """`` \\omega_{ij} = \\begin{cases}  1 & i = j , \\\\  0 & i \\ne j . \\\\  \\end{cases} `` """,
    )
end


@doc raw"""
Specialised Partitioned Additive Runge-Kutta integrator for index-two DAE systems.

This integrator solves the following system of equations for the internal stages,
```math
\begin{aligned}
Q_{n,i} &= q_{n} + h \sum \limits_{j=1}^{s} a_{ij} V_{n,j} + h \sum \limits_{j=1}^{r} \alpha_{ij} U_{n,j} , & i &= 1, ..., s , \\
P_{n,i} &= p_{n} + h \sum \limits_{j=1}^{s} a_{ij} F_{n,j} + h \sum \limits_{j=1}^{r} \alpha_{ij} G_{n,j} , & i &= 1, ..., s , \\
\tilde{Q}_{n,i} &= q_{n} + h \sum \limits_{j=1}^{s} \tilde{a}_{ij} V_{n,j} + h \sum \limits_{j=1}^{r} \tilde{\alpha}_{ij} U_{n,j} , & i &= 1, ..., r , \\
\tilde{P}_{n,i} &= p_{n} + h \sum \limits_{j=1}^{s} \tilde{a}_{ij} F_{n,j} + h \sum \limits_{j=1}^{r} \tilde{\alpha}_{ij} G_{n,j} , & i &= 1, ..., r , \\
0 &= \sum \limits_{j=1}^{r} \omega_{ij} \tilde{\Phi}_{n,j} , & i &= 1, ..., r-1 ,
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
p_{n+1} &= p_{n} + h \sum \limits_{i=1}^{s} b_{i} F_{n,i} + h \sum \limits_{i=1}^{r} \beta_{i} G_{n,i} , \\
0 &= \phi (q_{n+1}, p_{n+1}) .
\end{aligned}
```
"""
struct IntegratorSPARK{DT, TT, D, S, R, PT <: ParametersSPARK{DT,TT,D,S,R},
                                        ST <: NonlinearSolver{DT},
                                        IT <: InitialGuessIODE{TT}} <: AbstractIntegratorSPARK{DT,TT,D,S,R}
    params::PT
    solver::ST
    iguess::IT
    caches::CacheDict{PT}

    function IntegratorSPARK(params::ParametersSPARK{DT,TT,D,S,R}, solver::ST, iguess::IT, caches) where {DT,TT,D,S,R,ST,IT}
        new{DT, TT, D, S, R, typeof(params), ST, IT}(params, solver, iguess, caches)
    end

    function IntegratorSPARK{DT,D}(equations::NamedTuple, tableau::Union{TableauSPARK{TT},TableauSPARK{TT}}, Δt::TT) where {DT,TT,D}
        # @assert tableau.ρ == tableau.r-1

        # get number of stages
        S = tableau.s
        R = tableau.r
        P = tableau.ρ

        N = 2*D*S + 3*D*R

        if isdefined(tableau, :d) && length(tableau.d) > 0
            N += D
        end

        # create params
        params = ParametersSPARK{DT,D}(equations, tableau, Δt)

        # create cache dict
        caches = CacheDict(params)

        # create solver
        solver = create_nonlinear_solver(DT, N, params, caches)

        # create initial guess
        iguess = InitialGuessIODE(get_config(:ig_extrapolation), equations[:v̄], equations[:f̄], Δt)

        # create integrator
        IntegratorSPARK(params, solver, iguess, caches)
    end

    function IntegratorSPARK(equation::Union{IDAEProblem{DT}, LDAEProblem{DT}}, tableau::Union{TableauSPARK,TableauSPARK}, Δt=tstep(equation); kwargs...) where {DT}
        IntegratorSPARK{DT, ndims(equation)}(functions(equation), tableau, Δt; kwargs...)
    end
end


GeometricBase.nconstraints(::IntegratorSPARK{DT,TT,D}) where {DT,TT,D} = D


function Integrators.initialize!(int::IntegratorSPARK, sol::AtomicSolutionPDAE)
    sol.t̄ = sol.t - timestep(int)

    equation(int, :v̄)(sol.t, sol.q, sol.v)
    equation(int, :f̄)(sol.t, sol.q, sol.v, sol.f)

    initialize!(int.iguess, sol.t, sol.q, sol.p, sol.v, sol.f,
                            sol.t̄, sol.q̄, sol.p̄, sol.v̄, sol.f̄)
end

function initial_guess!(int::IntegratorSPARK{DT}, sol::AtomicSolutionPDAE{DT},
                        cache::IntegratorCacheSPARK{DT}=int.caches[DT]) where {DT}
    for i in eachstage(int)
        evaluate!(int.iguess, sol.q̄, sol.p̄, sol.v̄, sol.f̄,
                              sol.q, sol.p, sol.v, sol.f,
                              cache.q̃, cache.p̃, cache.ṽ, cache.f̃,
                              tableau(int).q.c[i], tableau(int).p.c[i])

        for k in eachdim(int)
            int.solver.x[2*(ndims(int)*(i-1)+k-1)+1] = (cache.q̃[k] - sol.q[k])/timestep(int)
            int.solver.x[2*(ndims(int)*(i-1)+k-1)+2] = (cache.p̃[k] - sol.p[k])/timestep(int)
        end
    end

    for i in 1:pstages(int)
        evaluate!(int.iguess, sol.q̄, sol.p̄, sol.v̄, sol.f̄,
                              sol.q, sol.p, sol.v, sol.f,
                              cache.q̃, cache.p̃, cache.ṽ, cache.f̃,
                              tableau(int).q̃.c[i], tableau(int).p̃.c[i])

        for k in eachdim(int)
            int.solver.x[2*ndims(int)*nstages(int)+3*(ndims(int)*(i-1)+k-1)+1] = (cache.q̃[k] - sol.q[k])/timestep(int)
            int.solver.x[2*ndims(int)*nstages(int)+3*(ndims(int)*(i-1)+k-1)+2] = (cache.p̃[k] - sol.p[k])/timestep(int)
            int.solver.x[2*ndims(int)*nstages(int)+3*(ndims(int)*(i-1)+k-1)+3] = cache.ṽ[k]
        end
    end

    # if isdefined(tableau(int), :λ) && tableau(int).λ.c[1] == 0
    #     for k in eachdim(int)
    #         int.solver.x[2*ndims(int)*nstages(int)+3*(k-1)+1] = cache.λ[k]
    #     end
    # end

    if isdefined(tableau(int), :d) && length(tableau(int).d) > 0
        for k in eachdim(int)
            int.solver.x[2*ndims(int)*nstages(int)+3*ndims(int)*pstages(int)+k] = 0
        end
    end
end


function compute_stages!(x::Vector{ST}, cache::IntegratorCacheSPARK{ST,D,S,R},
                                        params::ParametersSPARK{DT,TT,D,S,R}) where {ST,DT,TT,D,S,R}
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
        tpᵢ = params.t + params.Δt * params.tab.p.c[i]
        params.equs[:f](tpᵢ, cache.Qi[i], cache.Vi[i], cache.Fi[i])
        params.equs[:ϑ](tpᵢ, cache.Qi[i], cache.Vi[i], cache.Φi[i])
        cache.Φi[i] .-= cache.Pi[i]
    end

    for i in 1:R
        for k in 1:D
            # copy y to Y, Z and Λ
            cache.Yp[i][k] = x[2*D*S+3*(D*(i-1)+k-1)+1]
            cache.Zp[i][k] = x[2*D*S+3*(D*(i-1)+k-1)+2]
            cache.Vp[i][k] = x[2*D*S+3*(D*(i-1)+k-1)+3]
            cache.Λp[i][k] = cache.Vp[i][k]

            # compute Q and V
            cache.Qp[i][k] = params.q[k] + params.Δt * cache.Yp[i][k]
            cache.Pp[i][k] = params.p[k] + params.Δt * cache.Zp[i][k]
        end

        # compute f(X)
        tλᵢ = params.t + params.Δt * params.tab.λ.c[i]
        params.equs[:u](tλᵢ, cache.Qp[i], cache.Pp[i], cache.Vp[i], cache.Up[i])
        params.equs[:g](tλᵢ, cache.Qp[i], cache.Pp[i], cache.Vp[i], cache.Gp[i])
        params.equs[:ϕ](tλᵢ, cache.Qp[i], cache.Pp[i], cache.Φp[i])
    end

    if isdefined(params.tab, :d) && length(params.tab.d) > 0
        for k in 1:D
            cache.μ[k] = x[2*D*S+3*D*R+k]
        end
    end

    # compute q and p
    cache.q̃ .= params.q
    cache.p̃ .= params.p
    for i in 1:S
        cache.p̃ .+= params.Δt .* params.tab.p.b[i] .* cache.Fi[i]
    end
    for i in 1:R
        cache.q̃ .+= params.Δt .* params.tab.q.β[i] .* cache.Up[i]
        cache.p̃ .+= params.Δt .* params.tab.p.β[i] .* cache.Gp[i]
    end

    # compute ϕ(q,p)
    tλᵢ = params.t + params.Δt
    params.equs[:ϕ](tλᵢ, cache.q̃, cache.p̃, cache.ϕ̃)
end


"Compute stages of specialised partitioned additive Runge-Kutta methods for variational systems."
function Integrators.function_stages!(y::Vector{ST}, b::Vector{ST}, params::ParametersSPARK{DT,TT,D,S,R,P},
                                      caches::CacheDict) where {ST,DT,TT,D,S,R,P}

    # get cache for internal stages
    cache = caches[ST]

    compute_stages!(y, cache, params)

    # compute b = - [(Y-AV-AU), (Z-AF-AG), Φ]
    for i in 1:S
        for k in 1:D
            b[2*(D*(i-1)+k-1)+1] = - cache.Yi[i][k]
            b[2*(D*(i-1)+k-1)+2] = - cache.Zi[i][k]
            for j in 1:S
                b[2*(D*(i-1)+k-1)+2] += params.tab.p.a[i,j] * cache.Fi[j][k]
            end
            for j in 1:R
                b[2*(D*(i-1)+k-1)+1] += params.tab.q.α[i,j] * cache.Up[j][k]
                b[2*(D*(i-1)+k-1)+2] += params.tab.p.α[i,j] * cache.Gp[j][k]
            end
        end
    end

    # compute b = - [(Y-AV-AU), (Z-AF-AG)]
    for i in 1:R
        for k in 1:D
            b[2*D*S+3*(D*(i-1)+k-1)+1] = - cache.Yp[i][k]
            b[2*D*S+3*(D*(i-1)+k-1)+2] = - cache.Zp[i][k]
            b[2*D*S+3*(D*(i-1)+k-1)+3] = 0
            for j in 1:S
                b[2*D*S+3*(D*(i-1)+k-1)+2] += params.tab.p̃.a[i,j] * cache.Fi[j][k]
            end
            for j in 1:R
                b[2*D*S+3*(D*(i-1)+k-1)+1] += params.tab.q̃.α[i,j] * cache.Up[j][k]
                b[2*D*S+3*(D*(i-1)+k-1)+2] += params.tab.p̃.α[i,j] * cache.Gp[j][k]
            end
        end
    end

    # compute b = - ωΦ
    for i in 1:R-P
        for k in 1:D
            for j in 1:R
                b[2*D*S+3*(D*(i-1)+k-1)+3] -= params.tab.ω[i,j] * cache.Φp[j][k]
            end
            b[2*D*S+3*(D*(i-1)+k-1)+3] -= params.tab.ω[i,R+1] * cache.ϕ̃[k]
        end
    end

    # compute b = d_λ ⋅ V
    for i in R-P+1:R
        for k in 1:D
            for j in 1:R
                b[2*D*S+3*(D*(i-1)+k-1)+3] -= params.tab.δ[j] * cache.Vp[j][k]
            end
        end
    end

    if isdefined(params.tab, :d) && length(params.tab.d) > 0
        for i in 1:R
            for k in 1:D
                b[2*(D*(i-1)+k-1)+3] -= cache.μ[k] * params.tab.d[i] / params.tab.p.b[i]
            end
        end

        for k in 1:D
            b[2*D*S+3*D*R+k] = 0
            for i in 1:R
                b[2*D*S+3*D*R+k] -= cache.Vp[i][k] * params.tab.d[i]
            end
        end
    end
end
