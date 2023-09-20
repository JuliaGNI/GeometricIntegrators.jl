
"Holds the tableau of a Hamiltonian Partitioned Additive Runge-Kutta methods."
const TableauHPARK = AbstractTableauSPARK{:hpark}


struct HPARK{TT <: TableauHPARK} <: PSPARKMethod
    tableau::TT
end

tableau(method::HPARK) = method.tableau

solversize(problem::Union{PDAEProblem,HDAEProblem}, method::HPARK) =
    2 * ndims(problem) * nstages(method) + 3 * ndims(problem) * pstages(method)


@doc raw"""
Partitioned Additive Runge-Kutta integrator for Hamiltonian systems subject
to Dirac constraints *EXPERIMENTAL*.

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
p_{n+1} &= p_{n} + h \sum \limits_{i=1}^{s} b_{i} F_{n,i} + h \sum \limits_{i=1}^{r} \beta_{i} G_{n,i} .
\end{aligned}
```
"""
const IntegratorHPARK{DT,TT} = GeometricIntegrator{<:Union{PDAEProblem{DT,TT},HDAEProblem{DT,TT}}, <:HPARK}

function Base.show(io::IO, int::IntegratorHPARK)
    print(io, "\nPartitioned Additive Runge-Kutta integrator for Hamiltonian systems subject")
    print(io, "\nto Dirac constraints *EXPERIMENTAL*\n")
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
    method::HPARK, 
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
        functions(problem).v(cache.Vi[i], tqᵢ, cache.Qi[i], cache.Pi[i], parameters(solstep))
        functions(problem).f(cache.Fi[i], tpᵢ, cache.Qi[i], cache.Pi[i], parameters(solstep))
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
        functions(problem).u(cache.Up[i], tλᵢ, cache.Qp[i], cache.Pp[i], cache.Λp[i], parameters(solstep))
        functions(problem).g(cache.Gp[i], tλᵢ, cache.Qp[i], cache.Pp[i], cache.Λp[i], parameters(solstep))
        functions(problem).ϕ(cache.Φp[i], tλᵢ, cache.Qp[i], cache.Pp[i], parameters(solstep))
    end
end


# Compute stages of variational partitioned additive Runge-Kutta methods.
function residual!(
    b::Vector{ST},
    x::Vector{ST},
    solstep::SolutionStepPDAE, 
    problem::Union{PDAEProblem,HDAEProblem},
    method::HPARK, 
    caches::CacheDict) where {ST}

    # get cache for internal stages
    local cache = caches[ST]

    # number of internal stages
    local S = nstages(method)
    local R = pstages(method)
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

    # compute b = - [(Y-AV-AU), (Z-AF-AG), Φ]
    for i in 1:R
        for k in 1:D
            b[2*D*S+3*(D*(i-1)+k-1)+1] = - cache.Yp[i][k]
            b[2*D*S+3*(D*(i-1)+k-1)+2] = - cache.Zp[i][k]
            b[2*D*S+3*(D*(i-1)+k-1)+3] = - cache.Φp[i][k]
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

    # compute b = - [Λ₁-λ]
    if tableau(method).λ.c[1] == 0
        for k in 1:D
            b[2*D*S+3*(k-1)+3] = - cache.Λp[1][k] + solstep.λ̄[k]
        end
    end
end
