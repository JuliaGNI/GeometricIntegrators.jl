
const TableauHSPARK = AbstractTableauSPARK{:hspark}

struct HSPARK{TT <: Union{TableauSPARK,TableauHSPARK}} <: PSPARKMethod
    tableau::TT
end

HSPARK(method::SPARKMethod) = HSPARK(tableau(method))

tableau(method::HSPARK) = method.tableau

solversize(problem::Union{PDAEProblem,HDAEProblem}, method::Union{HSPARK,SPARKMethod}) =
    2 * ndims(problem) * nstages(method) + 3 * ndims(problem) * pstages(method)


@doc raw"""
Specialised Partitioned Additive Runge-Kutta integrator for Hamiltonian systems *EXPERIMENTAL*.

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
V_{n,i} &= \hphantom{-} \frac{\partial H}{\partial p} (Q_{n,i}, P_{n,i}) , & i &= 1, ..., s , \\
F_{n,i} &=           -  \frac{\partial H}{\partial q} (Q_{n,i}, P_{n,i}) , & i &= 1, ..., s , \\
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
const IntegratorHSPARK{DT,TT} = GeometricIntegrator{<:Union{PDAEProblem{DT,TT},HDAEProblem{DT,TT}}, <:HSPARK}

function Base.show(io::IO, int::IntegratorHSPARK)
    print(io, "\nSpecialised Partitioned Additive Runge-Kutta integrator for Hamiltonian systems *EXPERIMENTAL*\n")
    print(io, "   Timestep: $(timestep(problem))\n")
    print(io, "   Tableau:  $(description(method(int)))\n")
    print(io, "   $(string(method(int).q))")
    print(io, "   $(string(method(int).p))")
    # print(io, reference(method(int)))
end


function components!(
    x::Vector{ST},
    solstep::SolutionStepPDAE{DT,TT}, 
    problem::Union{PDAEProblem,HDAEProblem},
    method::HSPARK, 
    caches::CacheDict) where {ST,DT,TT}

    local cache = caches[ST]
    local S = nstages(method)
    local R = pstages(method)
    local D = ndims(problem)

    local tqᵢ::TT
    local tpᵢ::TT
    local tλᵢ::TT

    for i in 1:S
        for k in 1:D
            # copy x to Y, Z
            cache.Yi[i][k] = x[2*(D*(i-1)+k-1)+1]
            cache.Zi[i][k] = x[2*(D*(i-1)+k-1)+2]

            # compute Q and P
            cache.Qi[i][k] = solstep.q̄[k] + timestep(problem) * cache.Yi[i][k]
            cache.Pi[i][k] = solstep.p̄[k] + timestep(problem) * cache.Zi[i][k]
        end

        # compute f(X)
        tqᵢ = solstep.t̄ + timestep(problem) * tableau(method).q.c[i]
        tpᵢ = solstep.t̄ + timestep(problem) * tableau(method).p.c[i]
        functions(problem)[:v](cache.Vi[i], tqᵢ, cache.Qi[i], cache.Pi[i])
        functions(problem)[:f](cache.Fi[i], tpᵢ, cache.Qi[i], cache.Pi[i])
    end

    for i in 1:R
        for k in 1:D
            # copy y to Y, Z and Λ
            cache.Yp[i][k] = x[2*D*S+3*(D*(i-1)+k-1)+1]
            cache.Zp[i][k] = x[2*D*S+3*(D*(i-1)+k-1)+2]
            cache.Λp[i][k] = x[2*D*S+3*(D*(i-1)+k-1)+3]

            # compute Q and V
            cache.Qp[i][k] = solstep.q̄[k] + timestep(problem) * cache.Yp[i][k]
            cache.Pp[i][k] = solstep.p̄[k] + timestep(problem) * cache.Zp[i][k]
        end

        # compute f(X)
        tλᵢ = solstep.t̄ + timestep(problem) * tableau(method).λ.c[i]
        functions(problem)[:u](cache.Up[i], tλᵢ, cache.Qp[i], cache.Pp[i], cache.Λp[i])
        functions(problem)[:g](cache.Gp[i], tλᵢ, cache.Qp[i], cache.Pp[i], cache.Λp[i])
        functions(problem)[:ϕ](cache.Φp[i], tλᵢ, cache.Qp[i], cache.Pp[i])
    end

    # compute q and p
    cache.q̃ .= solstep.q̄
    cache.p̃ .= solstep.p̄
    for i in 1:S
        cache.q̃ .+= timestep(problem) .* tableau(method).q.b[i] .* cache.Vi[i]
        cache.p̃ .+= timestep(problem) .* tableau(method).p.b[i] .* cache.Fi[i]
    end
    for i in 1:R
        cache.q̃ .+= timestep(problem) .* tableau(method).q.β[i] .* cache.Up[i]
        cache.p̃ .+= timestep(problem) .* tableau(method).p.β[i] .* cache.Gp[i]
    end

    # compute ϕ(q,p)
    functions(problem)[:ϕ](cache.ϕ̃, solstep.t, cache.q̃, cache.p̃)
end


# Compute stages of specialised partitioned additive Runge-Kutta methods for variational systems.
function residual!(
    b::Vector{ST},
    x::Vector{ST},
    solstep::SolutionStepPDAE, 
    problem::Union{PDAEProblem,HDAEProblem},
    method::HSPARK, 
    caches::CacheDict) where {ST}

    # get cache for internal stages
    local cache = caches[ST]

    # number of internal stages
    local S = nstages(method)
    local R = pstages(method)
    local P = tableau(method).ρ
    local D = ndims(problem)

    # compute stages from nonlinear solver solution x
    components!(x, solstep, problem, method, caches)

    # compute b = - [(Y-AV-AU), (Z-AF-AG)]
    for i in 1:S
        for k in 1:D
            b[2*(D*(i-1)+k-1)+1] = - cache.Yi[i][k]
            b[2*(D*(i-1)+k-1)+2] = - cache.Zi[i][k]
            for j in 1:S
                b[2*(D*(i-1)+k-1)+1] += tableau(method).q.a[i,j] * cache.Vi[j][k]
                b[2*(D*(i-1)+k-1)+2] += tableau(method).p.a[i,j] * cache.Fi[j][k]
            end
            for j in 1:R
                b[2*(D*(i-1)+k-1)+1] += tableau(method).q.α[i,j] * cache.Up[j][k]
                b[2*(D*(i-1)+k-1)+2] += tableau(method).p.α[i,j] * cache.Gp[j][k]
            end
        end
    end

    # compute b = - [(Y-AV-AU), (Z-AF-AG), ωΦ]
    for i in 1:R
        for k in 1:D
            b[2*D*S+3*(D*(i-1)+k-1)+1] = - cache.Yp[i][k]
            b[2*D*S+3*(D*(i-1)+k-1)+2] = - cache.Zp[i][k]
            for j in 1:S
                b[2*D*S+3*(D*(i-1)+k-1)+1] += tableau(method).q̃.a[i,j] * cache.Vi[j][k]
                b[2*D*S+3*(D*(i-1)+k-1)+2] += tableau(method).p̃.a[i,j] * cache.Fi[j][k]
            end
            for j in 1:R
                b[2*D*S+3*(D*(i-1)+k-1)+1] += tableau(method).q̃.α[i,j] * cache.Up[j][k]
                b[2*D*S+3*(D*(i-1)+k-1)+2] += tableau(method).p̃.α[i,j] * cache.Gp[j][k]
            end
        end
    end
    for i in 1:R-P
        for k in 1:D
            b[2*D*S+3*(D*(i-1)+k-1)+3] = 0
            for j in 1:R
                b[2*D*S+3*(D*(i-1)+k-1)+3] -= tableau(method).ω[i,j] * cache.Φp[j][k]
            end
            b[2*D*S+3*(D*(i-1)+k-1)+3] -= tableau(method).ω[i,R+1] * cache.ϕ̃[k]
        end
    end

    # compute b = d_λ ⋅ Λ
    for i in R-P+1:R
        for k in 1:D
            b[2*D*S+3*(D*(R-1)+k-1)+3] = 0
            for j in 1:R
                b[2*D*S+3*(D*(i-1)+k-1)+3] -= tableau(method).δ[j] * cache.Λp[j][k]
            end
        end
    end
end
