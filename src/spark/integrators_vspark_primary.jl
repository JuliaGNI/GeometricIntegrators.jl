
const TableauVSPARKprimary = AbstractTableauSPARK{:vspark_primary}

function RungeKutta.check_symplecticity(tab::TableauVSPARKprimary{T}; atol=16*eps(T), rtol=16*eps(T)) where {T}
    s_a_qp  = [isapprox(tab.p.b[i] * tab.q.a[i,j] + tab.q.b[j] * tab.p.a[j,i], tab.p.b[i] * tab.q.b[j]; atol=atol, rtol=rtol) for i in 1:tab.s, j in 1:tab.s]
    s_α_q̃p̃  = [isapprox(tab.p.β[i] * tab.q̃.α[i,j] + tab.q.β[j] * tab.p̃.α[j,i], tab.p.β[i] * tab.q.β[j]; atol=atol, rtol=rtol) for i in 1:tab.r, j in 1:tab.r]
    s_αa_q̃p = [isapprox(tab.p.b[i] * tab.q.α[i,j] + tab.q.β[j] * tab.p̃.a[j,i], tab.p.b[i] * tab.q.β[j]; atol=atol, rtol=rtol) for i in 1:tab.s, j in 1:tab.r]
    s_αa_qp̃ = [isapprox(tab.q.b[i] * tab.p.α[i,j] + tab.p.β[j] * tab.q̃.a[j,i], tab.q.b[i] * tab.p.β[j]; atol=atol, rtol=rtol) for i in 1:tab.s, j in 1:tab.r]
    s_b_qp  = isapprox.(tab.q.b, tab.p.b; atol=atol, rtol=rtol)
    s_β_qp  = isapprox.(tab.q.β, tab.p.β; atol=atol, rtol=rtol)

    return (s_a_qp, s_α_q̃p̃, s_αa_q̃p, s_αa_qp̃, s_b_qp, s_β_qp)
end


function Integrators.symplecticity_conditions(::TableauVSPARKprimary)
    (
        """`` b^{p}_{i} b^{q}_{j} = b^{p}_{i} a^{q}_{ij} + b^{q}_{j} a^{p}_{ji} ``""",
        """`` \\beta^{p}_{i} \\beta^{q}_{j} = \\beta^{p}_{i} \\tilde{\\alpha}^{q}_{ij} + \\beta^{q}_{j} \\tilde{\\alpha}^{p}_{ji} ``""",
        """`` b^{p}_{i} \\beta^{q}_{j} = b^{p}_{i} \\alpha^{q}_{ij} + \\beta^{q}_{j} \\tilde{a}^{p}_{ji} ``""",
        """`` b^{q}_{i} \\beta^{p}_{j} = b^{q}_{i} \\alpha^{p}_{ij} + \\beta^{p}_{j} \\tilde{a}^{q}_{ji} ``""",
        """`` b^{q}_{i} = b^{p}_{i} ``""",
        """`` \\beta^{q}_{i} = \\beta^{p}_{i} ``""",
    )
end


function compute_conjugate_vspark_primary(a, b, b̄)
    ā = zero(a)
    for i in axes(ā,1)
        for j in axes(ā,2)
            ā[i,j] = b̄[j] / b[i] * ( b[i] - a[j,i] )
        end
    end
    return ā
end

function compute_ã_vspark_primary(α, β, b)
    s = length(b)
    s̃ = length(β)
    ã = zeros(eltype(α), s̃, s)
    for i in 1:s̃
        for j in 1:s
            ã[i,j] = b[j] / β[i] * ( β[i] - α[j,i] )
        end
    end
    return ã
end

function get_ã_vspark_primary(α_q, β_q, b_q, α_p, β_p, b_p)
    ã_q = compute_ã_vspark_primary(α_p, β_p, b_q)
    ã_p = compute_ã_vspark_primary(α_q, β_q, b_p)
    return (ã_q, ã_p)
end

function compute_α_vspark_primary(ã, b, β)
    s = length(b)
    s̃ = length(β)
    α = zeros(eltype(ã), s, s̃)
    for i in 1:s
        for j in 1:s̃
            α[i,j] = β[j] / b[i] * ( b[i] - ã[j,i] )
        end
    end
    return α
end

function get_α_vspark_primary(ã_q, b_q, β_q, ã_p, b_p, β_p)
    α_q = compute_α_vspark_primary(ã_p, b_p, β_q)
    α_p = compute_α_vspark_primary(ã_q, b_q, β_p)
    return (α_q, α_p)
end


struct VSPARKprimary{TT <: TableauVSPARKprimary} <: ISPARKMethod
    tableau::TT
end

tableau(method::VSPARKprimary) = method.tableau

solversize(problem::Union{IDAEProblem,LDAEProblem}, method::VSPARKprimary) =
    2 * ndims(problem) * nstages(method) + 2 * ndims(problem) * pstages(method)


@doc raw"""
Specialised Partitioned Additive Runge-Kutta integrator for Variational systems.

This integrator solves the following system of equations for the internal stages,
```math
\begin{aligned}
Q_{n,i} &= q_{n} + h \sum \limits_{j=1}^{s} a_{ij} V_{n,j} + h \sum \limits_{j=1}^{r} \alpha_{ij} U_{n,j} , & i &= 1, ..., s , \\
P_{n,i} &= p_{n} + h \sum \limits_{j=1}^{s} a_{ij} F_{n,j} + h \sum \limits_{j=1}^{r} \alpha_{ij} G_{n,j} , & i &= 1, ..., s , \\
\tilde{Q}_{n,i} &= q_{n} + h \sum \limits_{j=1}^{s} \tilde{a}_{ij} V_{n,j} + h \sum \limits_{j=1}^{r} \tilde{\alpha}_{ij} U_{n,j} , & i &= 1, ..., r , \\
\tilde{P}_{n,i} &= p_{n} + h \sum \limits_{j=1}^{s} \tilde{a}_{ij} F_{n,j} + h \sum \limits_{j=1}^{r} \tilde{\alpha}_{ij} G_{n,j} , & i &= 1, ..., r , \\
0 &= \sum \limits_{j=1}^{r} \omega_{ij} \tilde{\Phi}_{n,j} , & i &= 1, ..., r-1 , \\
0 &= \sum \limits_{i=1}^{r} \tilde{d}_i \, \Lambda_{n,i} ,
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
p_{n+1} &= p_{n} + h \sum \limits_{i=1}^{s} b_{i} F_{n,i} + h \sum \limits_{i=1}^{r} \beta_{i} G_{n,i} .
\end{aligned}
```
"""
const IntegratorVSPARKprimary{DT,TT} = GeometricIntegrator{<:VSPARKprimary, <:Union{IDAEProblem{DT,TT},LDAEProblem{DT,TT}}}

function Base.show(io::IO, int::IntegratorVSPARKprimary)
    print(io, "\nSpecialised Partitioned Additive Runge-Kutta integrator for Variational systems")
    print(io, "\nwith projection on primary constraint:\n")
    print(io, "   Timestep: $(timestep(int))\n")
    print(io, "   Tableau:  $(description(method(int)))\n")
    print(io, "   $(string(method(int).q))")
    print(io, "   $(string(method(int).p))")
    # print(io, reference(method(int)))
end


function initial_guess!(sol, history, params, int::GeometricIntegrator{<:VSPARKprimary,<:Union{IDAEProblem,LDAEProblem}})
    # get cache for internal stages
    local C = cache(int)

    for i in 1:nstages(int)
        # TODO: initialguess! should take two timesteps for c[i] of q and p tableau
        soltmp = (
            t = history.t[1] + timestep(int) * tableau(int).q.c[i],
            q = cache(int).Qi[i],
            p = cache(int).Pi[i],
            v = cache(int).Vi[i],
            f = cache(int).Fi[i],
        )
        solutionstep!(soltmp, history, problem(int), iguess(int))

        for k in 1:ndims(int)
            C.x[2*(ndims(int)*(i-1)+k-1)+1] =  C.Vi[i][k]
            C.x[2*(ndims(int)*(i-1)+k-1)+2] = (C.Pi[i][k] - sol.p[k]) / timestep(int)
        end
    end

    for i in 1:pstages(method(int))
        # TODO: initialguess! should take two timesteps for c[i] of q and p tableau
        soltmp = (
            t = history.t[1] + timestep(int) * tableau(int).q̃.c[i],
            q = cache(int).Qp[i],
            p = cache(int).Pp[i],
            v = cache(int).Vp[i],
            f = cache(int).Fp[i],
        )
        solutionstep!(soltmp, history, problem(int), iguess(int))

        for k in 1:ndims(int)
            C.x[2*ndims(int)*nstages(int)+2*(ndims(int)*(i-1)+k-1)+1] = 0
            C.x[2*ndims(int)*nstages(int)+2*(ndims(int)*(i-1)+k-1)+2] = (C.Pp[i][k] - sol.p[k]) / timestep(int)
        end
    end

    if isdefined(tableau(int), :λ) && tableau(int).λ.c[1] == 0
        for k in 1:ndims(int)
            C.x[2*ndims(int)*nstages(int)+2*(k-1)+1] = C.λ[k]
        end
    end

    if hasnullvector(method(int))
        for k in 1:ndims(int)
            C.x[2*ndims(int)*nstages(int)+2*ndims(int)*pstages(method(int))+k] = 0
        end
    end
end


function components!(x::AbstractVector{ST}, sol, params, int::GeometricIntegrator{<:VSPARKprimary,<:Union{IDAEProblem,LDAEProblem}}) where {ST}
    # get cache and number of internal stages
    local C = cache(int, ST)
    local S = nstages(int)
    local R = pstages(method(int))
    local D = ndims(int)

    # copy x to Vi and Zi
    for i in 1:S
        for k in 1:D
            C.Vi[i][k] = x[2*(D*(i-1)+k-1)+1]
            C.Zi[i][k] = x[2*(D*(i-1)+k-1)+2]
        end
    end

    # copy x to Λp and Zp
    for i in 1:R
        for k in 1:D
            C.Λp[i][k] = x[2*D*S+2*(D*(i-1)+k-1)+1]
            C.Zp[i][k] = x[2*D*S+2*(D*(i-1)+k-1)+2]
        end
        # tλᵢ = sol.t + timestep(int) * (tableau(int).λ.c[i] - 1)
        # params.f_u(tλᵢ, C.Qp[i], C.Pp[i], C.Λp[i], C.Up[i])
        C.Up[i] .= C.Λp[i]
    end

    for i in 1:S
        # compute Y
        C.Yi[i] .= 0
        for j in 1:S
            C.Yi[i] .+= tableau(int).q.a[i,j] .* C.Vi[j]
        end
        for j in 1:R
            C.Yi[i] .+= tableau(int).q.α[i,j] .* C.Up[j]
        end

        # compute Q and P
        C.Qi[i] .= sol.q .+ timestep(int) .* C.Yi[i]
        C.Pi[i] .= sol.p .+ timestep(int) .* C.Zi[i]

        # compute f(X)
        tpᵢ = sol.t + timestep(int) * (tableau(int).p.c[i] - 1)
        equations(int).f(C.Fi[i], tpᵢ, C.Qi[i], C.Vi[i], params)
        equations(int).ϑ(C.Φi[i], tpᵢ, C.Qi[i], C.Vi[i], params)

        C.Φi[i] .-= C.Pi[i]
    end

    for i in 1:R
        # compute Y
        C.Yp[i] .= 0
        for j in 1:S
            C.Yp[i] .+= tableau(int).q̃.a[i,j] .* C.Vi[j]
        end
        for j in 1:R
            C.Yp[i] .+= tableau(int).q̃.α[i,j] .* C.Up[j]
        end

        # compute Q and P
        C.Qp[i] .= sol.q .+ timestep(int) .* C.Yp[i]
        C.Pp[i] .= sol.p .+ timestep(int) .* C.Zp[i]

        # compute f(X)
        tλᵢ = sol.t + timestep(int) * (tableau(int).λ.c[i] - 1)
        equations(int).g(C.Gp[i], tλᵢ, C.Qp[i], C.Vp[i], C.Pp[i], C.Λp[i], params)
        equations(int).ϕ(C.Φp[i], tλᵢ, C.Qp[i], C.Vp[i], C.Pp[i], params)
    end

    if hasnullvector(method(int))
        for k in 1:D
            C.μ[k] = x[2*D*S+2*D*R+k]
        end
    end

    # compute q and p
    C.q̃ .= sol.q
    C.p̃ .= sol.p
    for i in 1:S
        C.q̃ .+= timestep(int) .* tableau(int).q.b[i] .* C.Vi[i]
        C.p̃ .+= timestep(int) .* tableau(int).p.b[i] .* C.Fi[i]
    end
    for i in 1:R
        C.q̃ .+= timestep(int) .* tableau(int).q.β[i] .* C.Up[i]
        C.p̃ .+= timestep(int) .* tableau(int).p.β[i] .* C.Gp[i]
    end

    # compute ϕ(q,p)
    equations(int).ϕ(C.ϕ̃, sol.t, C.q̃, C.ṽ, C.p̃, params)
end


# Compute stages of specialised partitioned additive Runge-Kutta methods for variational systems.
function residual!(b::AbstractVector{ST}, x::AbstractVector{ST}, sol, params, int::GeometricIntegrator{<:VSPARKprimary,<:Union{IDAEProblem,LDAEProblem}}) where {ST}
    # get cache and number of internal stages
    local C = cache(int, ST)
    local S = nstages(int)
    local R = pstages(method(int))
    local P = tableau(int).ρ
    local D = ndims(int)

    # compute stages from nonlinear solver solution x
    components!(x, sol, params, int)

    # compute b = [Φ, (Z-AF-AG)]
    for i in 1:S
        for k in 1:D
            b[2*(D*(i-1)+k-1)+1] = C.Φi[i][k]
            b[2*(D*(i-1)+k-1)+2] = C.Zi[i][k]
            for j in 1:S
                b[2*(D*(i-1)+k-1)+2] -= tableau(int).p.a[i,j] * C.Fi[j][k]
            end
            for j in 1:R
                b[2*(D*(i-1)+k-1)+2] -= tableau(int).p.α[i,j] * C.Gp[j][k]
            end
        end
    end

    # compute b = Z-AF-AG
    for i in 1:R
        for k in 1:D
            b[2*D*S+2*(D*(i-1)+k-1)+2] = C.Zp[i][k]
            for j in 1:S
                b[2*D*S+2*(D*(i-1)+k-1)+2] -= tableau(int).p̃.a[i,j] * C.Fi[j][k]
            end
            for j in 1:R
                b[2*D*S+2*(D*(i-1)+k-1)+2] -= tableau(int).p̃.α[i,j] * C.Gp[j][k]
            end
        end
    end

    # compute b = ωΦ
    for i in 1:R-P
        for k in 1:D
            b[2*D*S+2*(D*(i-1)+k-1)+1] = 0
            for j in 1:R
                b[2*D*S+2*(D*(i-1)+k-1)+1] += tableau(int).ω[i,j] * C.Φp[j][k]
            end
            b[2*D*S+2*(D*(i-1)+k-1)+1] += tableau(int).ω[i,R+1] * C.ϕ̃[k]
        end
    end

    # compute b = d_λ ⋅ Λ
    for i in 1:P
        for k in 1:D
            b[2*D*S+2*(D*(R-P+i-1)+k-1)+1] = 0
            for j in 1:R
                b[2*D*S+2*(D*(R-P+i-1)+k-1)+1] += tableau(int).δ[i,j] * C.Λp[j][k]
            end
        end
    end

    if hasnullvector(method(int))
        for i in 1:S
            for k in 1:D
                b[2*(D*(i-1)+k-1)+2] += C.μ[k] * tableau(int).d[i] / tableau(int).p.b[i]
            end
        end

        for k in 1:D
            b[2*D*S+2*D*R+k] = 0
            for i in 1:S
                b[2*D*S+2*D*R+k] += C.Vi[i][k] * tableau(int).d[i]
            end
        end
    end
end
