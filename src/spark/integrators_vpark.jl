
"Holds the tableau of an Variational Partitioned Additive Runge-Kutta method."
const TableauVPARK = AbstractTableauSPARK{:vpark}

struct VPARK{TT <: TableauVPARK} <: ISPARKMethod
    tableau::TT
end

tableau(method::VPARK) = method.tableau

solversize(problem::AbstractProblemIDAE, method::VPARK) =
    3 * ndims(problem) * nstages(method) + 3 * ndims(problem) * pstages(method)


@doc raw"""
Variational partitioned additive Runge-Kutta integrator.

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
const IntegratorVPARK{DT,TT} = GeometricIntegrator{<:VPARK, <:Union{IDAEProblem{DT,TT},LDAEProblem{DT,TT}}}

function Base.show(io::IO, int::IntegratorVPARK)
    print(io, "\nVariational partitioned additive Runge-Kutta integrator:\n")
    print(io, "   Timestep: $(timestep(int))\n")
    print(io, "   Tableau:  $(description(method(int)))\n")
    print(io, "   $(string(method(int).q))")
    print(io, "   $(string(method(int).p))")
    # print(io, reference(method(int)))
end


function components!(x::AbstractVector{ST}, sol, params, int::GeometricIntegrator{<:VPARK,<:Union{IDAEProblem,LDAEProblem}}) where {ST}
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
        equations(int).ϑ(C.Φi[i], tpᵢ, C.Qi[i], C.Vi[i], params)
        equations(int).f(C.Fi[i], tpᵢ, C.Qi[i], C.Vi[i], params)

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
end


# Compute stages of variational partitioned additive Runge-Kutta methods.
function residual!(b::AbstractVector{ST}, x::AbstractVector{ST}, sol, params, int::GeometricIntegrator{<:VPARK,<:Union{IDAEProblem,LDAEProblem}}) where {ST}
    # get cache and number of internal stages
    local C = cache(int, ST)
    local S = nstages(int)
    local R = pstages(method(int))
    local D = ndims(int)

    # compute stages from nonlinear solver solution x
    components!(x, sol, params, int)

    # compute b = - [(Y-AV-AU), (Z-AF-AG), Φ]
    for i in 1:S
        for k in 1:D
            b[3*(D*(i-1)+k-1)+1] = - C.Yi[i][k]
            b[3*(D*(i-1)+k-1)+2] = - C.Zi[i][k]
            b[3*(D*(i-1)+k-1)+3] = - C.Φi[i][k]
            for j in eachindex(C.Vi,C.Fi)
                b[3*(D*(i-1)+k-1)+1] += tableau(int).q.a[i,j] * C.Vi[j][k]
                b[3*(D*(i-1)+k-1)+2] += tableau(int).p.a[i,j] * C.Fi[j][k]
            end
            for j in eachindex(C.Up,C.Gp)
                b[3*(D*(i-1)+k-1)+1] += tableau(int).q.α[i,j] * C.Up[j][k]
                b[3*(D*(i-1)+k-1)+2] += tableau(int).p.α[i,j] * C.Gp[j][k]
            end
        end
    end

    # compute b = - [(Y-AV-AU), (Z-AF-AG), Φ]
    for i in 1:R
        for k in 1:D
            b[3*D*S+3*(D*(i-1)+k-1)+1] = - C.Yp[i][k]
            b[3*D*S+3*(D*(i-1)+k-1)+2] = - C.Zp[i][k]
            b[3*D*S+3*(D*(i-1)+k-1)+3] = - C.Φp[i][k]
            for j in eachindex(C.Vi,C.Fi)
                b[3*D*S+3*(D*(i-1)+k-1)+1] += tableau(int).q̃.a[i,j] * C.Vi[j][k]
                b[3*D*S+3*(D*(i-1)+k-1)+2] += tableau(int).p̃.a[i,j] * C.Fi[j][k]
            end
            for j in eachindex(C.Up,C.Gp)
                b[3*D*S+3*(D*(i-1)+k-1)+1] += tableau(int).q̃.α[i,j] * C.Up[j][k]
                b[3*D*S+3*(D*(i-1)+k-1)+2] += tableau(int).p̃.α[i,j] * C.Gp[j][k]
            end
        end
    end

    # compute b = - [Λ₁-λ]
    if tableau(int).λ.c[1] == 0
        for k in 1:D
            b[3*D*S+3*(k-1)+3] = - C.Λp[1][k] + sol.λ[k]
        end
    end

    if hasnullvector(method(int))
        for i in 1:S
            for k in eachindex(C.μ)
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
