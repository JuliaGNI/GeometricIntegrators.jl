"Holds all parameters of an Specialised Partitioned Additive Runge-Kutta method for variational systems subject to constraints."
struct SLRK{DT <: Number, DVT} <: LSPARKMethod
    name::Symbol
    o::Int
    s::Int
    r::Int

    q::Tableau{DT}
    p::Tableau{DT}

    q̃::Tableau{DT}
    p̃::Tableau{DT}

    ω::Matrix{DT}
    d::DVT

    function SLRK(name::Symbol, o::Int, s::Int,
                  q::Tableau{DT}, p::Tableau{DT},
                  q̃::Tableau{DT}, p̃::Tableau{DT},
                  ω::Matrix{DT}, d::DVT = nothing) where {DT, DVT <: Union{AbstractVector,Nothing}}

        @assert s > 0 "Number of stages s must be > 0"

        @assert s == q.s == p.s == q̃.s == p̃.s
        @assert size(ω,1) == s
        @assert size(ω,2) == s+1

        @assert d === nothing || length(d) == s

        new{DT,DVT}(name, o, s, s, q, p, q̃, p̃, ω, d)
    end
end

tableau(method::SLRK) = method

nstages(method::SLRK) = method.s
pstages(method::SLRK) = method.r

hasnullvector(method::SLRK{DT,Nothing}) where {DT} = false
hasnullvector(method::SLRK{DT,<:AbstractVector}) where {DT} = true

nonlinearsolversize(problem::LDAEProblem, method::SLRK) =
    4 * ndims(problem) * nstages(method)


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
const IntegratorSLRK{DT,TT} = Integrator{<:LDAEProblem{DT,TT}, <:SLRK}


# function Integrators.initsolver(::Newton, solstep::SolutionStepPDAE{DT}, problem::LDAEProblem, method::SLRK, caches::CacheDict) where {DT}
#     # create wrapper function f!(b,x)
#     f! = (b,x) -> function_stages!(b, x, solstep, problem, method, caches)

#     # create nonlinear solver
#     NewtonSolver(zero(caches[DT].x), zero(caches[DT].x), f!; linesearch = Backtracking(), config = Options(min_iterations = 1, x_abstol = 8eps(), f_abstol = 8eps()))
# end


function Base.show(io::IO, int::IntegratorSLRK)
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
    method::SLRK, 
    caches::CacheDict, 
    ::NonlinearSolver, 
    iguess::Union{InitialGuess,Extrapolation}) where {DT}

    # get caches for nonlinear solver vector and internal stages
    local x = caches[DT].x
    local Q = caches[DT].Qp
    local P = caches[DT].Pp
    local V = caches[DT].Vp
    local F = caches[DT].Fp

    # compute initial guess for internal stages
    for i in 1:pstages(method)
        # TODO: initialguess! should take two timesteps for c[i] of q and p tableau
        initialguess!(solstep.t̄[1] + timestep(problem) * method.p̃.c[i], Q[i], P[i], V[i], F[i], solstep, problem, iguess)
    end

    # assemble initial guess for nonlinear solver solution vector
    for i in 1:pstages(method)
        for k in 1:ndims(problem)
            offset = 4*(ndims(problem)*(i-1)+k-1)
            x[offset+1] = (Q[i][k] - solstep.q̄[1][k]) / timestep(problem)
            x[offset+2] = (P[i][k] - solstep.p̄[1][k]) / timestep(problem)
            x[offset+3] =  V[i][k]
            x[offset+4] = 0
        end
    end

    if hasnullvector(method)
        for k in 1:ndims(problem)
            x[4*ndims(problem)*pstages(method)+k] = 0
        end
    end
end


function compute_stages!(
    x::Vector{ST},
    solstep::SolutionStepPDAE{DT,TT}, 
    problem::LDAEProblem,
    method::SLRK, 
    caches::CacheDict) where {ST,DT,TT}

    local cache = caches[ST]
    local S = nstages(method)
    local D = ndims(problem)
    local t::TT

    for i in 1:nstages(method)
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
        t = solstep.t̄[1] + timestep(problem) * method.p.c[i]
        functions(problem).f(cache.Fp[i], t, cache.Qp[i], cache.Vp[i])
        functions(problem).g(cache.Gp[i], t, cache.Qp[i], cache.Vp[i], cache.Pp[i], cache.Λp[i])
        functions(problem).ϕ(cache.Φp[i], t, cache.Qp[i], cache.Vp[i], cache.Pp[i])
        functions(problem).ψ(cache.Ψp[i], t, cache.Qp[i], cache.Vp[i], cache.Pp[i], cache.Vp[i], cache.Fp[i])
    end

    if hasnullvector(method)
        for k in 1:D
            cache.μ[k] = x[4*D*S+k]
        end
    end

    # compute q and p
    cache.q̃ .= solstep.q̄[1]
    cache.p̃ .= solstep.p̄[1]
    for i in 1:nstages(method)
        cache.q̃ .+= timestep(problem) .* method.q.b[i] .* cache.Vp[i]
        cache.q̃ .+= timestep(problem) .* method.q̃.b[i] .* cache.Λp[i]
        cache.p̃ .+= timestep(problem) .* method.p.b[i] .* cache.Fp[i]
        cache.p̃ .+= timestep(problem) .* method.p̃.b[i] .* cache.Gp[i]
    end

    # compute ϕ(q,p)
    functions(problem).ϕ(cache.ϕ̃, solstep.t, cache.q̃, cache.ṽ, cache.p̃)
end


# Compute stages of specialised partitioned additive Runge-Kutta methods for variational systems.
function function_stages!(
    b::Vector{ST},
    x::Vector{ST},
    solstep::SolutionStepPDAE, 
    problem::LDAEProblem,
    method::SLRK, 
    caches::CacheDict) where {ST}

    # get cache for internal stages
    local cache = caches[ST]

    # number of internal stages
    local S = nstages(method)
    local D = ndims(problem)

    # compute stages from nonlinear solver solution x
    compute_stages!(x, solstep, problem, method, caches)

    # compute b = - [(Y-AV-AU), (Z-AF-AG), Φ, ωΨ]
    for i in 1:nstages(method)
        for k in 1:D
            b[4*(D*(i-1)+k-1)+1] = - cache.Yp[i][k]
            b[4*(D*(i-1)+k-1)+2] = - cache.Zp[i][k]
            b[4*(D*(i-1)+k-1)+3] = - cache.Φp[i][k]
            b[4*(D*(i-1)+k-1)+4] = method.ω[i,S+1] * cache.ϕ̃[k]

            for j in 1:nstages(method)
                b[4*(D*(i-1)+k-1)+1] += method.q.a[i,j] * cache.Vp[j][k]
                b[4*(D*(i-1)+k-1)+1] += method.q̃.a[i,j] * cache.Λp[j][k]
                b[4*(D*(i-1)+k-1)+2] += method.p.a[i,j] * cache.Fp[j][k]
                b[4*(D*(i-1)+k-1)+2] += method.p̃.a[i,j] * cache.Gp[j][k]
                b[4*(D*(i-1)+k-1)+4] += method.ω[i,j]   * cache.Ψp[j][k]
            end
        end
    end

    if hasnullvector(method)
        for i in 1:nstages(method)
            for k in 1:D
                b[4*(D*(i-1)+k-1)+2] += cache.μ[k] * method.d[i] / method.p.b[i]
                b[4*(D*(i-1)+k-1)+3] += cache.μ[k] * method.d[i] / method.p.b[i]
            end
        end

        for k in 1:D
            b[4*D*S+k] = 0
            for i in 1:nstages(method)
                b[4*D*S+k] -= cache.Vp[i][k] * method.d[i]
            end
        end
    end
end


function update_solution!(
    solstep::SolutionStepPDAE{DT,TT}, 
    problem::LDAEProblem,
    method::SLRK, 
    caches::CacheDict) where {DT,TT}

    # compute final update
    update_solution!(solstep.q, caches[DT].Vp, method.q.b, timestep(problem))
    update_solution!(solstep.q, caches[DT].Λp, method.q̃.b, timestep(problem))
    update_solution!(solstep.p, caches[DT].Fp, method.p.b, timestep(problem))
    update_solution!(solstep.p, caches[DT].Gp, method.p̃.b, timestep(problem))
end
