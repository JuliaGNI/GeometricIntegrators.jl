"Holds the tableau of an Specialised Partitioned Additive Runge-Kutta method for Variational systems."
struct VSPARKsecondary{DT <: Number, DVT} <: LSPARKMethod
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
    d::DVT

    function VSPARKsecondary(name::Symbol, o::Int, s::Int, r::Int,
                        q::CoefficientsSPARK{DT}, p::CoefficientsSPARK{DT},
                        q̃::CoefficientsSPARK{DT}, p̃::CoefficientsSPARK{DT},
                        ω::Matrix{DT}, d::DVT = nothing) where {DT, DVT <: Union{AbstractVector,Nothing}}

        @assert s > 0 "Number of stages s must be > 0"
        @assert r > 0 "Number of stages r must be > 0"

        @assert s==q.s==p.s==q̃.σ==p̃.σ
        @assert r==q.σ==p.σ==q̃.s==p̃.s
        @assert size(ω,1)==r
        @assert size(ω,2)==r+1

        @assert d === nothing || length(d) == r

        new{DT,DVT}(name, o, s, r, 0, q, p, q̃, p̃, ω, d)
    end
end

tableau(method::VSPARKsecondary) = method

nstages(method::VSPARKsecondary) = method.s
pstages(method::VSPARKsecondary) = method.r

hasnullvector(method::VSPARKsecondary{DT,Nothing}) where {DT} = false
hasnullvector(method::VSPARKsecondary{DT,<:AbstractVector}) where {DT} = true

solversize(problem::Union{IDAEProblem,LDAEProblem}, method::VSPARKsecondary) =
    4 * ndims(problem) * pstages(method)


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
const IntegratorVSPARKsecondary{DT,TT} = Integrator{<:LDAEProblem{DT,TT}, <:VSPARKsecondary}

function Base.show(io::IO, int::IntegratorVSPARKsecondary)
    print(io, "\nSpecialised Partitioned Additive Runge-Kutta integrator for degenerate")
    print(io, "\nvariational systems with projection on secondary constraint:\n")
    print(io, "   Timestep: $(timestep(problem))\n")
    print(io, "   Tableau:  $(description(method(int)))\n")
    print(io, "   $(string(method(int).q))")
    print(io, "   $(string(method(int).p))")
    # print(io, reference(method(int)))
end


function Integrators.initial_guess!(
    solstep::SolutionStepPDAE{DT}, 
    problem::LDAEProblem,
    method::VSPARKsecondary, 
    caches::CacheDict, 
    ::NonlinearSolver, 
    iguess::Union{InitialGuess,Extrapolation}) where {DT}

    cache = caches[DT]

    for i in 1:pstages(method)
        # TODO: initialguess! should take two timesteps for c[i] of q and p tableau
        initialguess!(solstep.t̄[1] + timestep(problem) * tableau(method).q̃.c[i], cache.Qp[i], cache.Pp[i], cache.Vp[i], cache.Fp[i], solstep, problem, iguess)

        for k in 1:ndims(problem)
            cache.x[4*(ndims(problem)*(i-1)+k-1)+1] = (cache.Qp[i][k] - solstep.q̄[1][k]) / timestep(problem)
            cache.x[4*(ndims(problem)*(i-1)+k-1)+2] = (cache.Pp[i][k] - solstep.p̄[1][k]) / timestep(problem)
            cache.x[4*(ndims(problem)*(i-1)+k-1)+3] =  cache.Vp[i][k]
            cache.x[4*(ndims(problem)*(i-1)+k-1)+4] = 0
        end
    end

    if hasnullvector(method)
        for k in 1:ndims(problem)
            cache.x[4*ndims(problem)*pstages(method)+k] = 0
        end
    end
end


function compute_stages!(
    x::Vector{ST},
    solstep::SolutionStepPDAE{DT,TT}, 
    problem::Union{IDAEProblem,LDAEProblem},
    method::VSPARKsecondary, 
    caches::CacheDict) where {ST,DT,TT}

    local cache = caches[ST]
    local S = nstages(method)
    local R = pstages(method)
    local D = ndims(problem)

    local t::TT

    for i in 1:R
        for k in 1:D
            # copy y to Y, Z and Λ
            cache.Yp[i][k] = x[4*(D*(i-1)+k-1)+1]
            cache.Zp[i][k] = x[4*(D*(i-1)+k-1)+2]
            cache.Vp[i][k] = x[4*(D*(i-1)+k-1)+3]
            cache.Λp[i][k] = x[4*(D*(i-1)+k-1)+4]

            # compute Q and P
            cache.Qp[i][k] = solstep.q̄[1][k] + timestep(problem) * cache.Yp[i][k]
            cache.Pp[i][k] = solstep.p̄[1][k] + timestep(problem) * cache.Zp[i][k]
        end

        # compute f(X)
        t = solstep.t̄[1] + timestep(problem) * tableau(method).p̃.c[i]
        functions(problem)[:f](cache.Fp[i], t, cache.Qp[i], cache.Vp[i])
        functions(problem)[:g](cache.Gp[i], t, cache.Qp[i], cache.Vp[i], cache.Pp[i], cache.Vp[i])
        functions(problem)[:g](cache.G̅p[i], t, cache.Qp[i], cache.Vp[i], cache.Pp[i], cache.Λp[i])

        cache.Hp[i] .= cache.Fp[i] .+ cache.Gp[i]

        functions(problem)[:ϕ](cache.Φp[i], t, cache.Qp[i], cache.Vp[i], cache.Pp[i])
        functions(problem)[:ψ](cache.Ψp[i], t, cache.Qp[i], cache.Vp[i], cache.Pp[i], cache.Vp[i], cache.Hp[i])
    end

    if hasnullvector(method)
        for k in 1:D
            cache.μ[k] = x[4*D*R+k]
        end
    end

    for i in 1:S
        # compute Q
        for k in 1:D
            cache.Yi[i][k] = 0
            for j in 1:R
                cache.Yi[i][k] += tableau(method).q.a[1][i,j] * cache.Vp[j][k]
                cache.Yi[i][k] += tableau(method).q.a[2][i,j] * cache.Λp[j][k]
            end
            cache.Qi[i][k] = solstep.q̄[1][k] + timestep(problem) * cache.Yi[i][k]
        end

        # compute f(X)
        t = solstep.t̄[1] + timestep(problem) * tableau(method).p.c[i]
        functions(problem)[:f](cache.Fi[i], t, cache.Qi[i], cache.Vi[i])
    end

    # compute q and p
    cache.q̃ .= solstep.q̄[1]
    cache.p̃ .= solstep.p̄[1]
    for i in 1:S
        cache.p̃ .+= timestep(problem) .* tableau(method).p.b[1][i] .* cache.Fi[i]
    end
    for i in 1:R
        cache.q̃ .+= timestep(problem) .* tableau(method).q.b[1][i] .* cache.Vp[i]
        cache.q̃ .+= timestep(problem) .* tableau(method).q.b[2][i] .* cache.Λp[i]
        cache.p̃ .+= timestep(problem) .* tableau(method).p.b[2][i] .* cache.Gp[i]
        cache.p̃ .+= timestep(problem) .* tableau(method).p.b[3][i] .* cache.G̅p[i]
    end

    # compute ϕ(q,p)
    functions(problem)[:ϕ](cache.ϕ̃, solstep.t, cache.q̃, cache.ṽ, cache.p̃)
end


# Compute stages of specialised partitioned additive Runge-Kutta methods for variational systems.
function function_stages!(
    b::Vector{ST},
    x::Vector{ST},
    solstep::SolutionStepPDAE, 
    problem::Union{IDAEProblem,LDAEProblem},
    method::VSPARKsecondary, 
    caches::CacheDict) where {ST}

    # get cache for internal stages
    local cache = caches[ST]

    # number of internal stages
    local S = nstages(method)
    local R = pstages(method)
    local D = ndims(problem)

    # compute stages from nonlinear solver solution x
    compute_stages!(x, solstep, problem, method, caches)

    # compute b = - [(Y-AV-AU), (Z-AF-AG), Φ, ωΨ]
    for i in 1:R
        for k in 1:D
            b[4*(D*(i-1)+k-1)+1] = - cache.Yp[i][k]
            b[4*(D*(i-1)+k-1)+2] = - cache.Zp[i][k]
            b[4*(D*(i-1)+k-1)+3] = - cache.Φp[i][k]
            b[4*(D*(i-1)+k-1)+4] = 0
            for j in 1:S
                b[4*(D*(i-1)+k-1)+2] += tableau(method).p̃.a[1][i,j] * cache.Fi[j][k]
            end
            for j in 1:R
                b[4*(D*(i-1)+k-1)+1] += tableau(method).q̃.a[1][i,j] * cache.Vp[j][k]
                b[4*(D*(i-1)+k-1)+1] += tableau(method).q̃.a[2][i,j] * cache.Λp[j][k]
                b[4*(D*(i-1)+k-1)+2] += tableau(method).p̃.a[2][i,j] * cache.Gp[j][k]
                b[4*(D*(i-1)+k-1)+2] += tableau(method).p̃.a[3][i,j] * cache.G̅p[j][k]
            end
            for j in 1:R
                b[4*(D*(i-1)+k-1)+4] -= tableau(method).ω[i,j] * cache.Ψp[j][k]
            end
            b[4*(D*(i-1)+k-1)+4] -= tableau(method).ω[i,R+1] * cache.ϕ̃[k]
        end
    end

    if hasnullvector(method)
        for i in 1:R
            for k in 1:D
                b[4*(D*(i-1)+k-1)+2] -= cache.μ[k] * tableau(method).d[i] / tableau(method).p.b[2][i]
            end
        end

        for k in 1:D
            b[4*D*R+k] = 0
            for i in 1:R
                b[4*D*R+k] -= cache.Vp[i][k] * tableau(method).d[i]
            end
        end
    end
end


function update_solution!(
    solstep::SolutionStepPDAE{DT,TT}, 
    problem::Union{IDAEProblem,LDAEProblem},
    method::VSPARKsecondary, 
    caches::CacheDict) where {DT,TT}

    # compute final update
    update_solution!(solstep.p, solstep.p̃, caches[DT].Fi, tableau(method).p.b[1], timestep(problem))

    # compute projection
    update_solution!(solstep.q, solstep.q̃, caches[DT].Vp, tableau(method).q.b[1], timestep(problem))
    update_solution!(solstep.q, solstep.q̃, caches[DT].Λp, tableau(method).q.b[2], timestep(problem))
    update_solution!(solstep.p, solstep.p̃, caches[DT].Gp, tableau(method).p.b[2], timestep(problem))
    update_solution!(solstep.p, solstep.p̃, caches[DT].G̅p, tableau(method).p.b[3], timestep(problem))
end
