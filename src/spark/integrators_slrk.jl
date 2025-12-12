"Holds all parameters of an Specialised Partitioned Additive Runge-Kutta method for variational systems subject to constraints."
struct SLRK{DT<:Number,DVT} <: LSPARKMethod
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
        ω::Matrix{DT}, d::DVT=nothing) where {DT,DVT<:Union{AbstractVector,Nothing}}

        @assert s > 0 "Number of stages s must be > 0"

        @assert s == q.s == p.s == q̃.s == p̃.s
        @assert size(ω, 1) == s
        @assert size(ω, 2) == s + 1

        @assert d === nothing || length(d) == s

        new{DT,DVT}(name, o, s, s, q, p, q̃, p̃, ω, d)
    end
end

tableau(method::SLRK) = method

nstages(method::SLRK) = method.s
pstages(method::SLRK) = method.r

hasnullvector(method::SLRK{DT,Nothing}) where {DT} = false
hasnullvector(method::SLRK{DT,<:AbstractVector}) where {DT} = true

solversize(problem::LDAEProblem, method::SLRK) =
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
const IntegratorSLRK{DT,TT} = GeometricIntegrator{<:LDAEProblem{DT,TT},<:SLRK}


# function Integrators.initsolver(::Newton, config::Options, solstep::SolutionStepPDAE{DT}, problem::LDAEProblem, method::SLRK, caches::CacheDict) where {DT}
#     # create wrapper function f!(b,x)
#     f! = (b,x) -> residual!(b, x, solstep, problem, method, caches)

#     # create nonlinear solver
#     NewtonSolver(zero(caches[DT].x), zero(caches[DT].x), f!; linesearch = Backtracking(), config = config)
# end


function Base.show(io::IO, int::IntegratorSLRK)
    print(io, "\nSpecialised Partitioned Additive Runge-Kutta integrator for degenerate")
    print(io, "\nvariational systems with projection on secondary constraint:\n")
    print(io, "   Timestep: $(timestep(int))\n")
    print(io, "   Tableau:  $(description(method(int)))\n")
    print(io, "   $(string(method(int).q))")
    print(io, "   $(string(method(int).p))")
    # print(io, reference(method(int)))
end


function initial_guess!(sol, history, params, int::GeometricIntegrator{<:SLRK,<:LDAEProblem})
    # get caches for nonlinear solver vector
    local x = cache(int).x

    # compute initial guess for internal stages
    for i in 1:pstages(method(int))
        # TODO: initialguess! should take two timesteps for c[i] of q and p tableau
        soltmp = (
            t=history.t[1] + timestep(int) * tableau(int).p̃.c[i],
            q=cache(int).Qp[i],
            p=cache(int).Pp[i],
            q̇=cache(int).Vp[i],
            ṗ=cache(int).Fp[i],
        )
        solutionstep!(soltmp, history, problem(int), iguess(int))
    end

    # assemble initial guess for nonlinear solver solution vector
    for i in 1:pstages(method(int))
        for k in 1:ndims(problem(int))
            offset = 4 * (ndims(int) * (i - 1) + k - 1)
            x[offset+1] = (cache(int).Qp[i][k] - sol.q[k]) / timestep(int)
            x[offset+2] = (cache(int).Pp[i][k] - sol.p[k]) / timestep(int)
            x[offset+3] = cache(int).Vp[i][k]
            x[offset+4] = 0
        end
    end

    if hasnullvector(method(int))
        for k in 1:ndims(int)
            x[4*ndims(int)*pstages(method(int))+k] = 0
        end
    end
end


function components!(x::AbstractVector{ST}, sol, params, int::GeometricIntegrator{<:SLRK}) where {ST}
    # get cache for internal stages
    local C = cache(int, ST)

    for i in eachindex(C.Yp, C.Zp, C.Vp, C.Λp)
        for k in eachindex(C.Yp[i], C.Zp[i], C.Vp[i], C.Λp[i])
            # copy y to Y, Z and Λ
            C.Yp[i][k] = x[4*(ndims(int)*(i-1)+k-1)+1]
            C.Zp[i][k] = x[4*(ndims(int)*(i-1)+k-1)+2]
            C.Vp[i][k] = x[4*(ndims(int)*(i-1)+k-1)+3]
            C.Λp[i][k] = x[4*(ndims(int)*(i-1)+k-1)+4]

            # compute Q and P
            C.Qp[i][k] = sol.q[k] + timestep(int) * C.Yp[i][k]
            C.Pp[i][k] = sol.p[k] + timestep(int) * C.Zp[i][k]
        end

        # compute f(X)
        tᵢ = sol.t + timestep(int) * (method(int).p.c[i] - 1)
        equations(int).f(C.Fp[i], tᵢ, C.Qp[i], C.Vp[i], params)
        equations(int).g(C.Gp[i], tᵢ, C.Qp[i], C.Vp[i], C.Pp[i], C.Λp[i], params)
        equations(int).ϕ(C.Φp[i], tᵢ, C.Qp[i], C.Vp[i], C.Pp[i], params)
        equations(int).ψ(C.Ψp[i], tᵢ, C.Qp[i], C.Vp[i], C.Pp[i], C.Vp[i], C.Fp[i], params)
    end

    if hasnullvector(method(int))
        for k in eachindex(C.μ)
            C.μ[k] = x[4*ndims(int)*nstages(int)+k]
        end
    end

    # compute q and p
    C.q̃ .= sol.q
    C.p̃ .= sol.p
    for i in 1:nstages(int)
        C.q̃ .+= timestep(int) .* method(int).q.b[i] .* C.Vp[i]
        C.q̃ .+= timestep(int) .* method(int).q̃.b[i] .* C.Λp[i]
        C.p̃ .+= timestep(int) .* method(int).p.b[i] .* C.Fp[i]
        C.p̃ .+= timestep(int) .* method(int).p̃.b[i] .* C.Gp[i]
    end

    # compute ϕ(q,p)
    equations(int).ϕ(C.ϕ̃, sol.t, C.q̃, C.ṽ, C.p̃, params)
end


# Compute stages of specialised partitioned additive Runge-Kutta methods for variational systems.
function residual!(b::AbstractVector{ST}, x::AbstractVector{ST}, sol, params, int::GeometricIntegrator{<:SLRK}) where {ST}
    # get cache and number of internal stages
    local C = cache(int, ST)
    local S = nstages(int)
    local D = ndims(int)

    # compute stages from nonlinear solver solution x
    components!(x, sol, params, int)

    # compute b = - [(Y-AV-AU), (Z-AF-AG), Φ, ωΨ]
    for i in 1:nstages(int)
        for k in 1:ndims(int)
            b[4*(D*(i-1)+k-1)+1] = -C.Yp[i][k]
            b[4*(D*(i-1)+k-1)+2] = -C.Zp[i][k]
            b[4*(D*(i-1)+k-1)+3] = -C.Φp[i][k]
            b[4*(D*(i-1)+k-1)+4] = method(int).ω[i, S+1] * C.ϕ̃[k]

            for j in 1:nstages(int)
                b[4*(D*(i-1)+k-1)+1] += method(int).q.a[i, j] * C.Vp[j][k]
                b[4*(D*(i-1)+k-1)+1] += method(int).q̃.a[i, j] * C.Λp[j][k]
                b[4*(D*(i-1)+k-1)+2] += method(int).p.a[i, j] * C.Fp[j][k]
                b[4*(D*(i-1)+k-1)+2] += method(int).p̃.a[i, j] * C.Gp[j][k]
                b[4*(D*(i-1)+k-1)+4] += method(int).ω[i, j] * C.Ψp[j][k]
            end
        end
    end

    if hasnullvector(method(int))
        for i in 1:nstages(int)
            for k in 1:ndims(int)
                b[4*(D*(i-1)+k-1)+2] += C.μ[k] * method(int).d[i] / method(int).p.b[i]
                b[4*(D*(i-1)+k-1)+3] += C.μ[k] * method(int).d[i] / method(int).p.b[i]
            end
        end

        for k in 1:ndims(int)
            b[4*D*S+k] = 0
            for i in 1:nstages(int)
                b[4*D*S+k] -= C.Vp[i][k] * method(int).d[i]
            end
        end
    end
end


function update!(sol, params, x::AbstractVector{DT}, int::GeometricIntegrator{<:SLRK,<:LDAEProblem}) where {DT}
    # compute vector field at internal stages
    components!(x, sol, params, int)

    # compute final update
    update!(sol.q, cache(int, DT).Vp, method(int).q, timestep(int))
    update!(sol.q, cache(int, DT).Λp, method(int).q̃, timestep(int))
    update!(sol.p, cache(int, DT).Fp, method(int).p, timestep(int))
    update!(sol.p, cache(int, DT).Gp, method(int).p̃, timestep(int))
end
