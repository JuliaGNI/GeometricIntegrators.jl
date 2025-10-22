"Holds the tableau of an Specialised Partitioned Additive Runge-Kutta method for Variational systems."
struct VSPARKsecondary{DT<:Number,DVT} <: LSPARKMethod
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
        ω::Matrix{DT}, d::DVT=nothing) where {DT,DVT<:Union{AbstractVector,Nothing}}

        @assert s > 0 "Number of stages s must be > 0"
        @assert r > 0 "Number of stages r must be > 0"

        @assert s == q.s == p.s == q̃.σ == p̃.σ
        @assert r == q.σ == p.σ == q̃.s == p̃.s
        @assert size(ω, 1) == r
        @assert size(ω, 2) == r + 1

        @assert d === nothing || length(d) == r

        new{DT,DVT}(name, o, s, r, 0, q, p, q̃, p̃, ω, d)
    end
end

tableau(method::VSPARKsecondary) = method

nstages(method::VSPARKsecondary) = method.s
pstages(method::VSPARKsecondary) = method.r

hasnullvector(method::VSPARKsecondary{DT,Nothing}) where {DT} = false
hasnullvector(method::VSPARKsecondary{DT,<:AbstractVector}) where {DT} = true

solversize(problem::AbstractProblemIDAE, method::VSPARKsecondary) =
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
const IntegratorVSPARKsecondary{DT,TT} = GeometricIntegrator{<:VSPARKsecondary,<:LDAEProblem{DT,TT}}

function Base.show(io::IO, int::IntegratorVSPARKsecondary)
    print(io, "\nSpecialised Partitioned Additive Runge-Kutta integrator for degenerate")
    print(io, "\nvariational systems with projection on secondary constraint:\n")
    print(io, "   Timestep: $(timestep(int))\n")
    print(io, "   Tableau:  $(description(method(int)))\n")
    print(io, "   $(string(method(int).q))")
    print(io, "   $(string(method(int).p))")
    # print(io, reference(method(int)))
end


function initial_guess!(sol, history, params, int::GeometricIntegrator{<:VSPARKsecondary,<:Union{IDAEProblem,LDAEProblem}})
    # get cache for internal stages
    local C = cache(int)

    for i in 1:pstages(method(int))
        # TODO: initialguess! should take two timesteps for c[i] of q and p tableau
        soltmp = (
            t=history.t[1] + timestep(int) * tableau(int).q̃.c[i],
            q=cache(int).Qp[i],
            p=cache(int).Pp[i],
            q̇=cache(int).Vp[i],
            ṗ=cache(int).Fp[i],
        )
        solutionstep!(soltmp, history, problem(int), iguess(int))

        for k in 1:ndims(int)
            C.x[4*(ndims(int)*(i-1)+k-1)+1] = (C.Qp[i][k] - sol.q[k]) / timestep(int)
            C.x[4*(ndims(int)*(i-1)+k-1)+2] = (C.Pp[i][k] - sol.p[k]) / timestep(int)
            C.x[4*(ndims(int)*(i-1)+k-1)+3] = C.Vp[i][k]
            C.x[4*(ndims(int)*(i-1)+k-1)+4] = 0
        end
    end

    if hasnullvector(method(int))
        for k in 1:ndims(int)
            C.x[4*ndims(int)*pstages(method(int))+k] = 0
        end
    end
end


function components!(x::AbstractVector{ST}, sol, params, int::GeometricIntegrator{<:VSPARKsecondary,<:Union{IDAEProblem,LDAEProblem}}) where {ST}
    # get cache and number of internal stages
    local C = cache(int, ST)
    local S = nstages(int)
    local R = pstages(method(int))
    local D = ndims(int)

    for i in 1:R
        for k in 1:D
            # copy y to Y, Z and Λ
            C.Yp[i][k] = x[4*(D*(i-1)+k-1)+1]
            C.Zp[i][k] = x[4*(D*(i-1)+k-1)+2]
            C.Vp[i][k] = x[4*(D*(i-1)+k-1)+3]
            C.Λp[i][k] = x[4*(D*(i-1)+k-1)+4]
        end

        # compute Q and P
        C.Qp[i] .= sol.q .+ timestep(int) .* C.Yp[i]
        C.Pp[i] .= sol.p .+ timestep(int) .* C.Zp[i]

        # compute f(X)
        t = sol.t + timestep(int) * (tableau(int).p̃.c[i] - 1)
        equations(int).f(C.Fp[i], t, C.Qp[i], C.Vp[i], params)
        equations(int).g(C.Gp[i], t, C.Qp[i], C.Vp[i], C.Pp[i], C.Vp[i], params)
        equations(int).g(C.G̅p[i], t, C.Qp[i], C.Vp[i], C.Pp[i], C.Λp[i], params)

        C.Hp[i] .= C.Fp[i] .+ C.Gp[i]

        equations(int).ϕ(C.Φp[i], t, C.Qp[i], C.Vp[i], C.Pp[i], params)
        equations(int).ψ(C.Ψp[i], t, C.Qp[i], C.Vp[i], C.Pp[i], C.Vp[i], C.Hp[i], params)
    end

    if hasnullvector(method(int))
        for k in 1:D
            C.μ[k] = x[4*D*R+k]
        end
    end

    for i in 1:S
        # compute Y
        for k in 1:D
            C.Yi[i][k] = 0
            for j in 1:R
                C.Yi[i][k] += tableau(int).q.a[1][i, j] * C.Vp[j][k]
                C.Yi[i][k] += tableau(int).q.a[2][i, j] * C.Λp[j][k]
            end
        end

        # compute Q
        C.Qi[i] .= sol.q .+ timestep(int) .* C.Yi[i]

        # compute f(X)
        t = sol.t + timestep(int) * (tableau(int).p.c[i] - 1)
        equations(int).f(C.Fi[i], t, C.Qi[i], C.Vi[i], params)
    end

    # compute q and p
    C.q̃ .= sol.q
    C.p̃ .= sol.p
    for i in 1:S
        C.p̃ .+= timestep(int) .* tableau(int).p.b[1][i] .* C.Fi[i]
    end
    for i in 1:R
        C.q̃ .+= timestep(int) .* tableau(int).q.b[1][i] .* C.Vp[i]
        C.q̃ .+= timestep(int) .* tableau(int).q.b[2][i] .* C.Λp[i]
        C.p̃ .+= timestep(int) .* tableau(int).p.b[2][i] .* C.Gp[i]
        C.p̃ .+= timestep(int) .* tableau(int).p.b[3][i] .* C.G̅p[i]
    end

    # compute ϕ(q,p)
    equations(int).ϕ(C.ϕ̃, sol.t, C.q̃, C.ṽ, C.p̃, params)
end


# Compute stages of specialised partitioned additive Runge-Kutta methods for variational systems.
function residual!(b::AbstractVector{ST}, x::AbstractVector{ST}, sol, params, int::GeometricIntegrator{<:VSPARKsecondary,<:Union{IDAEProblem,LDAEProblem}}) where {ST}
    # get cache and number of internal stages
    local C = cache(int, ST)
    local S = nstages(int)
    local R = pstages(method(int))
    local D = ndims(int)

    # compute stages from nonlinear solver solution x
    components!(x, sol, params, int)

    # compute b = - [(Y-AV-AU), (Z-AF-AG), Φ, ωΨ]
    for i in 1:R
        for k in 1:D
            b[4*(D*(i-1)+k-1)+1] = -C.Yp[i][k]
            b[4*(D*(i-1)+k-1)+2] = -C.Zp[i][k]
            b[4*(D*(i-1)+k-1)+3] = -C.Φp[i][k]
            b[4*(D*(i-1)+k-1)+4] = 0
            for j in 1:S
                b[4*(D*(i-1)+k-1)+2] += tableau(int).p̃.a[1][i, j] * C.Fi[j][k]
            end
            for j in 1:R
                b[4*(D*(i-1)+k-1)+1] += tableau(int).q̃.a[1][i, j] * C.Vp[j][k]
                b[4*(D*(i-1)+k-1)+1] += tableau(int).q̃.a[2][i, j] * C.Λp[j][k]
                b[4*(D*(i-1)+k-1)+2] += tableau(int).p̃.a[2][i, j] * C.Gp[j][k]
                b[4*(D*(i-1)+k-1)+2] += tableau(int).p̃.a[3][i, j] * C.G̅p[j][k]
            end
            for j in 1:R
                b[4*(D*(i-1)+k-1)+4] -= tableau(int).ω[i, j] * C.Ψp[j][k]
            end
            b[4*(D*(i-1)+k-1)+4] -= tableau(int).ω[i, R+1] * C.ϕ̃[k]
        end
    end

    if hasnullvector(method(int))
        for i in 1:R
            for k in 1:D
                b[4*(D*(i-1)+k-1)+2] -= C.μ[k] * tableau(int).d[i] / tableau(int).p.b[2][i]
            end
        end

        for k in 1:D
            b[4*D*R+k] = 0
            for i in 1:R
                b[4*D*R+k] -= C.Vp[i][k] * tableau(int).d[i]
            end
        end
    end
end


function update!(sol, params, x::AbstractVector{DT}, int::GeometricIntegrator{<:VSPARKsecondary}) where {DT}
    # compute vector field at internal stages
    components!(x, sol, params, int)

    # compute final update
    update!(sol.p, cache(int, DT).Fi, tableau(int).p.b[1], timestep(int))

    # compute projection
    update!(sol.q, cache(int, DT).Vp, tableau(int).q.b[1], timestep(int))
    update!(sol.q, cache(int, DT).Λp, tableau(int).q.b[2], timestep(int))
    update!(sol.p, cache(int, DT).Gp, tableau(int).p.b[2], timestep(int))
    update!(sol.p, cache(int, DT).G̅p, tableau(int).p.b[3], timestep(int))
end
