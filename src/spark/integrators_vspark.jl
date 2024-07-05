
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
const IntegratorVSPARK{DT,TT} = GeometricIntegrator{<:VSPARK, <:Union{IDAEProblem{DT,TT},LDAEProblem{DT,TT}}}

function Base.show(io::IO, int::IntegratorVSPARK)
    print(io, "\nSpecialised Partitioned Additive Runge-Kutta integrator for Variational systems:\n")
    print(io, "   Timestep: $(timestep(int))\n")
    print(io, "   Tableau:  $(description(method(int)))\n")
    print(io, "   $(string(method(int).q))")
    print(io, "   $(string(method(int).p))")
    # print(io, reference(method(int)))
end


function initial_guess!(sol, history, params, int::GeometricIntegrator{<:VSPARK,<:Union{IDAEProblem,LDAEProblem}})
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
            C.x[3*(ndims(int)*(i-1)+k-1)+1] = (C.Qi[i][k] - sol.q[k]) / timestep(int)
            C.x[3*(ndims(int)*(i-1)+k-1)+2] = (C.Pi[i][k] - sol.p[k]) / timestep(int)
            C.x[3*(ndims(int)*(i-1)+k-1)+3] =  C.Vi[i][k]
        end

        # Quick fix for dirty implementation of F function
        C.Vi[i] .= 0
        C.Fi[i] .= 0
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
            C.x[3*ndims(int)*nstages(int)+3*(ndims(int)*(i-1)+k-1)+1] = (C.Qp[i][k] - sol.q[k]) / timestep(int)
            C.x[3*ndims(int)*nstages(int)+3*(ndims(int)*(i-1)+k-1)+2] = (C.Pp[i][k] - sol.p[k]) / timestep(int)
            C.x[3*ndims(int)*nstages(int)+3*(ndims(int)*(i-1)+k-1)+3] = 0
        end

        # Quick fix for dirty implementation of F function
        C.Vp[i] .= 0
        C.Fp[i] .= 0
    end

    # TODO: Check indices !!!
    # if isdefined(tableau(int), :λ) && tableau(int).λ.c[1] == 0
    #     for k in 1:ndims(int)
    #         C.x[3*ndims(int)*nstages(int)+3*(k-1)+3] = solstep(int).λ[k]
    #     end
    # end

    if hasnullvector(method(int))
        for k in 1:ndims(int)
            C.x[3*ndims(int)*nstages(int)+3*ndims(int)*pstages(method(int))+k] = 0
        end
    end
end


function components!(x::AbstractVector{ST}, sol, params, int::GeometricIntegrator{<:VSPARK,<:Union{IDAEProblem,LDAEProblem}}) where {ST}
    # get cache and number of internal stages
    local C = cache(int, ST)
    local S = nstages(int)
    local R = pstages(method(int))
    local D = ndims(int)

    for i in 1:S
        for k in 1:D
            # copy x to Y, Z
            C.Yi[i][k] = x[3*(D*(i-1)+k-1)+1]
            C.Zi[i][k] = x[3*(D*(i-1)+k-1)+2]
            C.Vi[i][k] = x[3*(D*(i-1)+k-1)+3]
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
        for k in 1:D
            # copy y to Y, Z and Λ
            C.Yp[i][k] = x[3*D*S+3*(D*(i-1)+k-1)+1]
            C.Zp[i][k] = x[3*D*S+3*(D*(i-1)+k-1)+2]
            C.Λp[i][k] = x[3*D*S+3*(D*(i-1)+k-1)+3]
        end

        # compute Q and V
        C.Qp[i] .= sol.q .+ timestep(int) .* C.Yp[i]
        C.Pp[i] .= sol.p .+ timestep(int) .* C.Zp[i]

        # compute f(X)
        tλᵢ = sol.t + timestep(int) * (tableau(int).λ.c[i] - 1)
        equations(int).u(C.Up[i], tλᵢ, C.Qp[i], C.Vp[i], C.Pp[i], C.Λp[i], params)
        equations(int).g(C.Gp[i], tλᵢ, C.Qp[i], C.Vp[i], C.Pp[i], C.Λp[i], params)
        equations(int).ϕ(C.Φp[i], tλᵢ, C.Qp[i], C.Vp[i], C.Pp[i], params)
    end

    if hasnullvector(method(int))
        for k in 1:D
            C.μ[k] = x[3*D*S+3*D*R+k]
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
function residual!(b::AbstractVector{ST}, x::AbstractVector{ST}, sol, params, int::GeometricIntegrator{<:VSPARK,<:Union{IDAEProblem,LDAEProblem}}) where {ST}
    # get cache and number of internal stages
    local C = cache(int, ST)
    local S = nstages(method(int))
    local R = pstages(method(int))
    local P = tableau(int).ρ
    local D = ndims(int)

    # compute stages from nonlinear solver solution x
    components!(x, sol, params, int)

    # compute b = - [(Y-AV-AU), (Z-AF-AG), Φ]
    for i in 1:S
        for k in 1:D
            b[3*(D*(i-1)+k-1)+1] = - C.Yi[i][k]
            b[3*(D*(i-1)+k-1)+2] = - C.Zi[i][k]
            b[3*(D*(i-1)+k-1)+3] = - C.Φi[i][k]
            for j in 1:S
                b[3*(D*(i-1)+k-1)+1] += tableau(int).q.a[i,j] * C.Vi[j][k]
                b[3*(D*(i-1)+k-1)+2] += tableau(int).p.a[i,j] * C.Fi[j][k]
            end
            for j in 1:R
                b[3*(D*(i-1)+k-1)+1] += tableau(int).q.α[i,j] * C.Up[j][k]
                b[3*(D*(i-1)+k-1)+2] += tableau(int).p.α[i,j] * C.Gp[j][k]
            end
        end
    end

    # compute b = - [(Y-AV-AU), (Z-AF-AG)]
    for i in 1:R
        for k in 1:D
            b[3*D*S+3*(D*(i-1)+k-1)+1] = - C.Yp[i][k]
            b[3*D*S+3*(D*(i-1)+k-1)+2] = - C.Zp[i][k]
            b[3*D*S+3*(D*(i-1)+k-1)+3] = 0
            for j in 1:S
                b[3*D*S+3*(D*(i-1)+k-1)+1] += tableau(int).q̃.a[i,j] * C.Vi[j][k]
                b[3*D*S+3*(D*(i-1)+k-1)+2] += tableau(int).p̃.a[i,j] * C.Fi[j][k]
            end
            for j in 1:R
                b[3*D*S+3*(D*(i-1)+k-1)+1] += tableau(int).q̃.α[i,j] * C.Up[j][k]
                b[3*D*S+3*(D*(i-1)+k-1)+2] += tableau(int).p̃.α[i,j] * C.Gp[j][k]
            end
        end
    end

    # compute b = - ωΦ
    for i in 1:R-P
        for k in 1:D
            for j in 1:R
                b[3*D*S+3*(D*(i-1)+k-1)+3] -= tableau(int).ω[i,j] * C.Φp[j][k]
            end
            b[3*D*S+3*(D*(i-1)+k-1)+3] -= tableau(int).ω[i,R+1] * C.ϕ̃[k]
        end
    end

    # compute b = d_λ ⋅ Λ
    for i in R-P+1:R
        for k in 1:D
            for j in 1:R
                b[3*D*S+3*(D*(i-1)+k-1)+3] -= tableau(int).δ[j] * C.Λp[j][k]
            end
        end
    end

    if hasnullvector(method(int))
        for i in 1:S
            for k in 1:D
                b[3*(D*(i-1)+k-1)+3] -= C.μ[k] * tableau(int).d[i] / tableau(int).p.b[i]
            end
        end

        for k in 1:D
            b[3*D*S+3*D*R+k] = 0
            for i in 1:S
                b[3*D*S+3*D*R+k] -= C.Vi[i][k] * tableau(int).d[i]
            end
        end
    end
end
