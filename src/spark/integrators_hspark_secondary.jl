
struct HSPARKsecondary{DT <: Number, DVT} <: HSPARKMethod
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

    function HSPARKsecondary(name::Symbol, o::Int, s::Int, r::Int,
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

tableau(method::HSPARKsecondary) = method

nstages(method::HSPARKsecondary) = method.s
pstages(method::HSPARKsecondary) = method.r

hasnullvector(method::HSPARKsecondary{DT,Nothing}) where {DT} = false
hasnullvector(method::HSPARKsecondary{DT,<:AbstractVector}) where {DT} = true

solversize(problem::HDAEProblem, method::HSPARKsecondary) =
    2 * ndims(problem) * nstages(method) + 4 * ndims(problem) * pstages(method)


@doc raw"""
Specialised Partitioned Additive Runge-Kutta integrator for Hamiltonian systems
subject to Dirac constraints with projection on secondary constraint *EXPERIMENTAL*.

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
const IntegratorHSPARKsecondary{DT,TT} = GeometricIntegrator{<:HSPARKsecondary, <:HDAEProblem{DT,TT}}


function Base.show(io::IO, int::IntegratorHSPARKsecondary)
    print(io, "\nSpecialised Partitioned Additive Runge-Kutta integrator for Hamiltonian systems")
    print(io, "\nsubject to Dirac constraints with projection on secondary constraint *EXPERIMENTAL*\n")
    print(io, "   Timestep: $(timestep(int))\n")
    print(io, "   Tableau:  $(description(method(int)))\n")
    print(io, "   $(string(method(int).q))")
    print(io, "   $(string(method(int).p))")
    # print(io, reference(method(int)))
end


function initial_guess!(int::GeometricIntegrator{<:Union{HSPARKsecondary},<:Union{HDAEProblem,PDAEProblem}})
    # get cache for internal stages
    local C = cache(int)
    local sol = current(solstep(int))

    for i in 1:nstages(int)
        initialguess!(sol.t + timestep(int) * (tableau(int).q.c[i] - 1), C.Qi[i], C.Pi[i], C.Vi[i], C.Fi[i], solstep(int), problem(int), iguess(int))

        for k in 1:ndims(int)
            C.x[2*(ndims(int)*(i-1)+k-1)+1] = (C.Qi[i][k] - sol.q[k]) / timestep(int)
            C.x[2*(ndims(int)*(i-1)+k-1)+2] = (C.Pi[i][k] - sol.p[k]) / timestep(int)
        end
    end

    for i in 1:pstages(method(int))
        # TODO: initialguess! should take two timesteps for c[i] of q and p tableau
        initialguess!(sol.t + timestep(int) * (tableau(int).q̃.c[i] - 1), C.Qp[i], C.Pp[i], C.Vp[i], C.Fp[i], solstep(int), problem(int), iguess(int))

        for k in 1:ndims(int)
            C.x[2*ndims(int)*nstages(int)+4*(ndims(int)*(i-1)+k-1)+1] = (C.Qp[i][k] - sol.q[k]) / timestep(int)
            C.x[2*ndims(int)*nstages(int)+4*(ndims(int)*(i-1)+k-1)+2] = (C.Pp[i][k] - sol.p[k]) / timestep(int)
            C.x[2*ndims(int)*nstages(int)+4*(ndims(int)*(i-1)+k-1)+3] = 0
            C.x[2*ndims(int)*nstages(int)+4*(ndims(int)*(i-1)+k-1)+4] = 0
        end
    end

    if isdefined(tableau(int), :λ) && tableau(int).λ.c[1] == 0
        for k in 1:ndims(int)
            C.x[2*ndims(int)*nstages(int)+4*ndims(int)*pstages(method(int))+k] = 0
        end
    end
end


function components!(x::AbstractVector{ST}, sol, params, int::GeometricIntegrator{<:HSPARKsecondary,<:Union{PDAEProblem,HDAEProblem}}) where {ST}
    # get cache and number of internal stages
    local C = cache(int, ST)
    local S = nstages(method(int))
    local R = pstages(method(int))
    local D = ndims(int)

    for i in 1:S
        for k in 1:D
            # copy x to Y, Z
            C.Yi[i][k] = x[2*(D*(i-1)+k-1)+1]
            C.Zi[i][k] = x[2*(D*(i-1)+k-1)+2]
        end

        # compute Q and P
        C.Qi[i] .= sol.q .+ timestep(int) .* C.Yi[i]
        C.Pi[i] .= sol.p .+ timestep(int) .* C.Zi[i]

        # compute f(X)
        t = sol.t + timestep(int) * (tableau(int).p.c[i] - 1)
        equations(int).v(C.Vi[i], t, C.Qi[i], C.Pi[i], params)
        equations(int).f(C.Fi[i], t, C.Qi[i], C.Pi[i], params)
    end

    for i in 1:R
        for k in 1:D
            # copy y to Y, Z and Λ
            C.Yp[i][k] = x[2*D*S+4*(D*(i-1)+k-1)+1]
            C.Zp[i][k] = x[2*D*S+4*(D*(i-1)+k-1)+2]
            C.Up[i][k] = x[2*D*S+4*(D*(i-1)+k-1)+3]
            C.Λp[i][k] = x[2*D*S+4*(D*(i-1)+k-1)+4]
        end

        # compute Q and V
        C.Qp[i] .= sol.q .+ timestep(int) .* C.Yp[i]
        C.Pp[i] .= sol.p .+ timestep(int) .* C.Zp[i]

        # compute f(X)
        t = sol.t + timestep(int) * (tableau(int).p̃.c[i] - 1)
        equations(int).v(C.Vp[i], t, C.Qp[i], C.Pp[i], params)
        equations(int).f(C.Fp[i], t, C.Qp[i], C.Pp[i], params)

        equations(int).g(C.Gp[i], t, C.Qp[i], C.Pp[i], C.Up[i], params)
        equations(int).g(C.G̅p[i], t, C.Qp[i], C.Pp[i], C.Λp[i], params)

        equations(int).ϕ(C.Φp[i], t, C.Qp[i], C.Pp[i], params)
        equations(int).ψ(C.Ψp[i], t, C.Qp[i], C.Pp[i], C.Vp[i], C.Fp[i], params)
    end

    # if hasnullvector(method(int))
    #     for k in 1:D
    #         C.μ[k] = x[2*D*S+4*D*R+k]
    #     end
    # end

    # compute q and p
    C.q̃ .= sol.q
    C.p̃ .= sol.p
    for i in 1:S
        C.q̃ .+= timestep(int) .* tableau(int).q.b[1][i] .* C.Vi[i]
        C.p̃ .+= timestep(int) .* tableau(int).p.b[1][i] .* C.Fi[i]
    end
    for i in 1:R
        C.q̃ .+= timestep(int) .* tableau(int).q.b[2][i] .* C.Up[i]
        C.q̃ .+= timestep(int) .* tableau(int).q.b[3][i] .* C.Λp[i]
        C.p̃ .+= timestep(int) .* tableau(int).p.b[2][i] .* C.Gp[i]
        C.p̃ .+= timestep(int) .* tableau(int).p.b[3][i] .* C.G̅p[i]
    end

    # compute ϕ(q,p)
    equations(int).ϕ(C.ϕ̃, sol.t, C.q̃, C.p̃, params)
end


# Compute stages of specialised partitioned additive Runge-Kutta methods for variational systems.
function residual!(b::AbstractVector{ST}, x::AbstractVector{ST}, sol, params, int::GeometricIntegrator{<:HSPARKsecondary,<:Union{PDAEProblem,HDAEProblem}}) where {ST}
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
            b[2*(D*(i-1)+k-1)+1] = - C.Yi[i][k]
            b[2*(D*(i-1)+k-1)+2] = - C.Zi[i][k]
            for j in 1:S
                b[2*(D*(i-1)+k-1)+1] += tableau(int).q.a[1][i,j] * C.Vi[j][k]
                b[2*(D*(i-1)+k-1)+2] += tableau(int).p.a[1][i,j] * C.Fi[j][k]
            end
            for j in 1:R
                b[2*(D*(i-1)+k-1)+1] += tableau(int).q.a[2][i,j] * C.Up[j][k]
                b[2*(D*(i-1)+k-1)+1] += tableau(int).q.a[3][i,j] * C.Λp[j][k]
                b[2*(D*(i-1)+k-1)+2] += tableau(int).p.a[2][i,j] * C.Gp[j][k]
                b[2*(D*(i-1)+k-1)+2] += tableau(int).p.a[3][i,j] * C.G̅p[j][k]
            end
        end
    end

    # compute b = - [(Y-AV-AU), (Z-AF-AG), Φ, ωΨ]
    for i in 1:R
        for k in 1:D
            b[2*D*S+4*(D*(i-1)+k-1)+1] = - C.Yp[i][k]
            b[2*D*S+4*(D*(i-1)+k-1)+2] = - C.Zp[i][k]
            b[2*D*S+4*(D*(i-1)+k-1)+3] = - C.Φp[i][k]
            # b[2*D*S+4*(D*(i-1)+k-1)+4] = - C.Ψp[i][k]
            b[2*D*S+4*(D*(i-1)+k-1)+4] = 0
            for j in 1:S
                b[2*D*S+4*(D*(i-1)+k-1)+1] += tableau(int).q̃.a[1][i,j] * C.Vi[j][k]
                b[2*D*S+4*(D*(i-1)+k-1)+2] += tableau(int).p̃.a[1][i,j] * C.Fi[j][k]
            end
            for j in 1:R
                b[2*D*S+4*(D*(i-1)+k-1)+1] += tableau(int).q̃.a[2][i,j] * C.Up[j][k]
                b[2*D*S+4*(D*(i-1)+k-1)+1] += tableau(int).q̃.a[3][i,j] * C.Λp[j][k]
                b[2*D*S+4*(D*(i-1)+k-1)+2] += tableau(int).p̃.a[2][i,j] * C.Gp[j][k]
                b[2*D*S+4*(D*(i-1)+k-1)+2] += tableau(int).p̃.a[3][i,j] * C.G̅p[j][k]
            end
            for j in 1:R
                b[2*D*S+4*(D*(i-1)+k-1)+4] -= tableau(int).ω[i,j] * C.Ψp[j][k]
            end
            b[2*D*S+4*(D*(i-1)+k-1)+4] -= tableau(int).ω[i,R+1] * C.ϕ̃[k]
        end
    end

    # if hasnullvector(method(int))
    #     for i in 1:R
    #         for k in 1:D
    #             b[2*D*S+4*(D*(i-1)+k-1)+2] -= C.μ[k] * tableau(int).d[i] / tableau(int).p.b[2][i]
    #         end
    #     end
    
    #     for k in 1:D
    #         b[2*D*S+4*D*R+k] = 0
    #         for i in 1:R
    #             b[2*D*S+4*D*R+k] -= C.Vp[i][k] * tableau(int).d[i]
    #         end
    #     end
    # end
end


function update!(sol, params, x::AbstractVector{DT}, int::GeometricIntegrator{<:HSPARKsecondary}) where {DT}
    # compute vector field at internal stages
    components!(x, sol, params, int)

    # compute final update
    update!(sol.q, cache(int, DT).Vi, tableau(int).q.b[1], timestep(int))
    update!(sol.p, cache(int, DT).Fi, tableau(int).p.b[1], timestep(int))

    # compute projection
    update!(sol.q, cache(int, DT).Up, tableau(int).q.b[2], timestep(int))
    update!(sol.q, cache(int, DT).Λp, tableau(int).q.b[3], timestep(int))
    update!(sol.p, cache(int, DT).Gp, tableau(int).p.b[2], timestep(int))
    update!(sol.p, cache(int, DT).G̅p, tableau(int).p.b[3], timestep(int))
end
