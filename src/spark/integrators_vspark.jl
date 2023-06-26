
const TableauVSPARK = AbstractTableauSPARK{:vspark}

struct VSPARK{TT <: Union{TableauSPARK,TableauVSPARK}} <: ISPARKMethod
    tableau::TT
end

VSPARK(method::SPARKMethod) = VSPARK(tableau(method))

tableau(method::VSPARK) = method.tableau

solversize(problem::Union{IDAEProblem,LDAEProblem}, method::VSPARK) =
    3 * ndims(problem) * nstages(method) + 3 * ndims(problem) * pstages(method)


@doc raw"""
Specialised Partitioned Additive Runge-Kutta integrator for Variational systems.

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
const IntegratorVSPARK{DT,TT} = Integrator{<:Union{IDAEProblem{DT,TT},LDAEProblem{DT,TT}}, <:VSPARK}

function Base.show(io::IO, int::IntegratorVSPARK)
    print(io, "\nSpecialised Partitioned Additive Runge-Kutta integrator for Variational systems:\n")
    print(io, "   Timestep: $(timestep(problem))\n")
    print(io, "   Tableau:  $(description(method(int)))\n")
    print(io, "   $(string(method(int).q))")
    print(io, "   $(string(method(int).p))")
    # print(io, reference(method(int)))
end


function Integrators.initial_guess!(
    solstep::SolutionStepPDAE{DT}, 
    problem::Union{IDAEProblem,LDAEProblem},
    method::VSPARK, 
    caches::CacheDict, 
    ::NonlinearSolver, 
    iguess::Union{InitialGuess,Extrapolation}) where {DT}

    cache = caches[DT]

    for i in 1:nstages(method)
        # TODO: initialguess! should take two timesteps for c[i] of q and p tableau
        initialguess!(solstep.t̄ + timestep(problem) * tableau(method).q.c[i], cache.Qi[i], cache.Pi[i], cache.Vi[i], cache.Fi[i], solstep, problem, iguess)

        for k in 1:ndims(problem)
            cache.x[3*(ndims(problem)*(i-1)+k-1)+1] = (cache.Qi[i][k] - solstep.q̄[k]) / timestep(problem)
            cache.x[3*(ndims(problem)*(i-1)+k-1)+2] = (cache.Pi[i][k] - solstep.p̄[k]) / timestep(problem)
            cache.x[3*(ndims(problem)*(i-1)+k-1)+3] =  cache.Vi[i][k]
        end

        # Quick fix for dirty implementation of F function
        cache.Vi[i] .= 0
        cache.Fi[i] .= 0
    end

    for i in 1:pstages(method)
        # TODO: initialguess! should take two timesteps for c[i] of q and p tableau
        initialguess!(solstep.t̄ + timestep(problem) * tableau(method).q̃.c[i], cache.Qp[i], cache.Pp[i], cache.Vp[i], cache.Fp[i], solstep, problem, iguess)

        for k in 1:ndims(problem)
            cache.x[3*ndims(problem)*nstages(method)+3*(ndims(problem)*(i-1)+k-1)+1] = (cache.Qp[i][k] - solstep.q̄[k]) / timestep(problem)
            cache.x[3*ndims(problem)*nstages(method)+3*(ndims(problem)*(i-1)+k-1)+2] = (cache.Pp[i][k] - solstep.p̄[k]) / timestep(problem)
            cache.x[3*ndims(problem)*nstages(method)+3*(ndims(problem)*(i-1)+k-1)+3] = 0
        end

        # Quick fix for dirty implementation of F function
        cache.Vp[i] .= 0
        cache.Fp[i] .= 0
    end

    # TODO: Check indices !!!
    # if isdefined(tableau(method), :λ) && tableau(method).λ.c[1] == 0
    #     for k in 1:ndims(problem)
    #         cache.x[3*ndims(problem)*nstages(method)+3*(k-1)+3] = solstep.λ[k]
    #     end
    # end

    if hasnullvector(method)
        for k in 1:ndims(problem)
            cache.x[3*ndims(problem)*nstages(method)+3*ndims(problem)*pstages(method)+k] = 0
        end
    end
end


function components!(
    x::Vector{ST},
    solstep::SolutionStepPDAE{DT,TT}, 
    problem::Union{IDAEProblem,LDAEProblem},
    method::VSPARK, 
    caches::CacheDict) where {ST,DT,TT}

    local cache = caches[ST]
    local S = nstages(method)
    local R = pstages(method)
    local D = ndims(problem)

    local tpᵢ::TT

    for i in 1:S
        for k in 1:D
            # copy x to Y, Z
            cache.Yi[i][k] = x[3*(D*(i-1)+k-1)+1]
            cache.Zi[i][k] = x[3*(D*(i-1)+k-1)+2]
            cache.Vi[i][k] = x[3*(D*(i-1)+k-1)+3]

            # compute Q and P
            cache.Qi[i][k] = solstep.q̄[k] + timestep(problem) * cache.Yi[i][k]
            cache.Pi[i][k] = solstep.p̄[k] + timestep(problem) * cache.Zi[i][k]
        end

        # compute f(X)
        tpᵢ = solstep.t̄ + timestep(problem) * tableau(method).p.c[i]
        functions(problem)[:f](cache.Fi[i], tpᵢ, cache.Qi[i], cache.Vi[i])
        functions(problem)[:ϑ](cache.Φi[i], tpᵢ, cache.Qi[i], cache.Vi[i])
        cache.Φi[i] .-= cache.Pi[i]
    end

    for i in 1:R
        for k in 1:D
            # copy y to Y, Z and Λ
            cache.Yp[i][k] = x[3*D*S+3*(D*(i-1)+k-1)+1]
            cache.Zp[i][k] = x[3*D*S+3*(D*(i-1)+k-1)+2]
            cache.Λp[i][k] = x[3*D*S+3*(D*(i-1)+k-1)+3]

            # compute Q and V
            cache.Qp[i][k] = solstep.q̄[k] + timestep(problem) * cache.Yp[i][k]
            cache.Pp[i][k] = solstep.p̄[k] + timestep(problem) * cache.Zp[i][k]
        end

        # compute f(X)
        tλᵢ = solstep.t̄ + timestep(problem) * tableau(method).λ.c[i]
        functions(problem)[:u](cache.Up[i], tλᵢ, cache.Qp[i], cache.Vp[i], cache.Pp[i], cache.Λp[i])
        functions(problem)[:g](cache.Gp[i], tλᵢ, cache.Qp[i], cache.Vp[i], cache.Pp[i], cache.Λp[i])
        functions(problem)[:ϕ](cache.Φp[i], tλᵢ, cache.Qp[i], cache.Vp[i], cache.Pp[i])
    end

    if hasnullvector(method)
        for k in 1:D
            cache.μ[k] = x[3*D*S+3*D*R+k]
        end
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
    functions(problem)[:ϕ](cache.ϕ̃, solstep.t, cache.q̃, cache.ṽ, cache.p̃)
end


# Compute stages of specialised partitioned additive Runge-Kutta methods for variational systems.
function residual!(
    b::Vector{ST},
    x::Vector{ST},
    solstep::SolutionStepPDAE, 
    problem::Union{IDAEProblem,LDAEProblem},
    method::VSPARK, 
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

    # compute b = - [(Y-AV-AU), (Z-AF-AG), Φ]
    for i in 1:S
        for k in 1:D
            b[3*(D*(i-1)+k-1)+1] = - cache.Yi[i][k]
            b[3*(D*(i-1)+k-1)+2] = - cache.Zi[i][k]
            b[3*(D*(i-1)+k-1)+3] = - cache.Φi[i][k]
            for j in 1:S
                b[3*(D*(i-1)+k-1)+1] += tableau(method).q.a[i,j] * cache.Vi[j][k]
                b[3*(D*(i-1)+k-1)+2] += tableau(method).p.a[i,j] * cache.Fi[j][k]
            end
            for j in 1:R
                b[3*(D*(i-1)+k-1)+1] += tableau(method).q.α[i,j] * cache.Up[j][k]
                b[3*(D*(i-1)+k-1)+2] += tableau(method).p.α[i,j] * cache.Gp[j][k]
            end
        end
    end

    # compute b = - [(Y-AV-AU), (Z-AF-AG)]
    for i in 1:R
        for k in 1:D
            b[3*D*S+3*(D*(i-1)+k-1)+1] = - cache.Yp[i][k]
            b[3*D*S+3*(D*(i-1)+k-1)+2] = - cache.Zp[i][k]
            b[3*D*S+3*(D*(i-1)+k-1)+3] = 0
            for j in 1:S
                b[3*D*S+3*(D*(i-1)+k-1)+1] += tableau(method).q̃.a[i,j] * cache.Vi[j][k]
                b[3*D*S+3*(D*(i-1)+k-1)+2] += tableau(method).p̃.a[i,j] * cache.Fi[j][k]
            end
            for j in 1:R
                b[3*D*S+3*(D*(i-1)+k-1)+1] += tableau(method).q̃.α[i,j] * cache.Up[j][k]
                b[3*D*S+3*(D*(i-1)+k-1)+2] += tableau(method).p̃.α[i,j] * cache.Gp[j][k]
            end
        end
    end

    # compute b = - ωΦ
    for i in 1:R-P
        for k in 1:D
            for j in 1:R
                b[3*D*S+3*(D*(i-1)+k-1)+3] -= tableau(method).ω[i,j] * cache.Φp[j][k]
            end
            b[3*D*S+3*(D*(i-1)+k-1)+3] -= tableau(method).ω[i,R+1] * cache.ϕ̃[k]
        end
    end

    # compute b = d_λ ⋅ Λ
    for i in R-P+1:R
        for k in 1:D
            for j in 1:R
                b[3*D*S+3*(D*(i-1)+k-1)+3] -= tableau(method).δ[j] * cache.Λp[j][k]
            end
        end
    end

    if hasnullvector(method)
        for i in 1:S
            for k in 1:D
                b[3*(D*(i-1)+k-1)+3] -= cache.μ[k] * tableau(method).d[i] / tableau(method).p.b[i]
            end
        end

        for k in 1:D
            b[3*D*S+3*D*R+k] = 0
            for i in 1:S
                b[3*D*S+3*D*R+k] -= cache.Vi[i][k] * tableau(method).d[i]
            end
        end
    end
end
