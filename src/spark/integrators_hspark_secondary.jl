
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
const IntegratorHSPARKsecondary{DT,TT} = GeometricIntegrator{<:HDAEProblem{DT,TT}, <:HSPARKsecondary}


function Base.show(io::IO, int::IntegratorHSPARKsecondary)
    print(io, "\nSpecialised Partitioned Additive Runge-Kutta integrator for Hamiltonian systems")
    print(io, "\nsubject to Dirac constraints with projection on secondary constraint *EXPERIMENTAL*\n")
    print(io, "   Timestep: $(timestep(problem))\n")
    print(io, "   Tableau:  $(description(method(int)))\n")
    print(io, "   $(string(method(int).q))")
    print(io, "   $(string(method(int).p))")
    # print(io, reference(method(int)))
end


function Integrators.initial_guess!(
    solstep::SolutionStepPDAE{DT}, 
    problem::HDAEProblem,
    method::HSPARKsecondary, 
    caches::CacheDict, 
    ::NonlinearSolver, 
    iguess::Union{InitialGuess,Extrapolation}) where {DT}

    cache = caches[DT]

    for i in 1:nstages(method)
        initialguess!(solstep.t̄ + timestep(problem) * tableau(method).q.c[i], cache.Qi[i], cache.Pi[i], cache.Vi[i], cache.Fi[i], solstep, problem, iguess)

        for k in 1:ndims(problem)
            cache.x[2*(ndims(problem)*(i-1)+k-1)+1] = (cache.Qi[i][k] - solstep.q̄[k]) / timestep(problem)
            cache.x[2*(ndims(problem)*(i-1)+k-1)+2] = (cache.Pi[i][k] - solstep.p̄[k]) / timestep(problem)
        end
    end

    for i in 1:pstages(method)
        # TODO: initialguess! should take two timesteps for c[i] of q and p tableau
        initialguess!(solstep.t̄ + timestep(problem) * tableau(method).q̃.c[i], cache.Qp[i], cache.Pp[i], cache.Vp[i], cache.Fp[i], solstep, problem, iguess)

        for k in 1:ndims(problem)
            cache.x[2*ndims(problem)*nstages(method)+4*(ndims(problem)*(i-1)+k-1)+1] = (cache.Qp[i][k] - solstep.q̄[k]) / timestep(problem)
            cache.x[2*ndims(problem)*nstages(method)+4*(ndims(problem)*(i-1)+k-1)+2] = (cache.Pp[i][k] - solstep.p̄[k]) / timestep(problem)
            cache.x[2*ndims(problem)*nstages(method)+4*(ndims(problem)*(i-1)+k-1)+3] = 0
            cache.x[2*ndims(problem)*nstages(method)+4*(ndims(problem)*(i-1)+k-1)+4] = 0
        end
    end

    if isdefined(tableau(method), :λ) && tableau(method).λ.c[1] == 0
        for k in 1:ndims(problem)
            cache.x[2*ndims(problem)*nstages(method)+4*ndims(problem)*pstages(method)+k] = 0
        end
    end
end


function components!(
    x::Vector{ST},
    solstep::SolutionStepPDAE{DT,TT}, 
    problem::HDAEProblem,
    method::HSPARKsecondary, 
    caches::CacheDict) where {ST,DT,TT}

    local cache = caches[ST]
    local S = nstages(method)
    local R = pstages(method)
    local D = ndims(problem)

    local t::TT

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
        t = solstep.t̄ + timestep(problem) * tableau(method).p.c[i]
        functions(problem).v(cache.Vi[i], t, cache.Qi[i], cache.Pi[i], parameters(solstep))
        functions(problem).f(cache.Fi[i], t, cache.Qi[i], cache.Pi[i], parameters(solstep))
    end

    for i in 1:R
        for k in 1:D
            # copy y to Y, Z and Λ
            cache.Yp[i][k] = x[2*D*S+4*(D*(i-1)+k-1)+1]
            cache.Zp[i][k] = x[2*D*S+4*(D*(i-1)+k-1)+2]
            cache.Up[i][k] = x[2*D*S+4*(D*(i-1)+k-1)+3]
            cache.Λp[i][k] = x[2*D*S+4*(D*(i-1)+k-1)+4]

            # compute Q and V
            cache.Qp[i][k] = solstep.q̄[k] + timestep(problem) * cache.Yp[i][k]
            cache.Pp[i][k] = solstep.p̄[k] + timestep(problem) * cache.Zp[i][k]
        end

        # compute f(X)
        t = solstep.t̄ + timestep(problem) * tableau(method).p̃.c[i]
        functions(problem).v(cache.Vp[i], t, cache.Qp[i], cache.Pp[i], parameters(solstep))
        functions(problem).f(cache.Fp[i], t, cache.Qp[i], cache.Pp[i], parameters(solstep))

        functions(problem).g(cache.Gp[i], t, cache.Qp[i], cache.Pp[i], cache.Up[i], parameters(solstep))
        functions(problem).g(cache.G̅p[i], t, cache.Qp[i], cache.Pp[i], cache.Λp[i], parameters(solstep))

        functions(problem).ϕ(cache.Φp[i], t, cache.Qp[i], cache.Pp[i], parameters(solstep))
        functions(problem).ψ(cache.Ψp[i], t, cache.Qp[i], cache.Pp[i], cache.Vp[i], cache.Fp[i], parameters(solstep))
    end

    # if hasnullvector(method)
    #     for k in 1:D
    #         cache.μ[k] = x[2*D*S+4*D*R+k]
    #     end
    # end

    # compute q and p
    cache.q̃ .= solstep.q̄
    cache.p̃ .= solstep.p̄
    for i in 1:S
        cache.q̃ .+= timestep(problem) .* tableau(method).q.b[1][i] .* cache.Vi[i]
        cache.p̃ .+= timestep(problem) .* tableau(method).p.b[1][i] .* cache.Fi[i]
    end
    for i in 1:R
        cache.q̃ .+= timestep(problem) .* tableau(method).q.b[2][i] .* cache.Up[i]
        cache.q̃ .+= timestep(problem) .* tableau(method).q.b[3][i] .* cache.Λp[i]
        cache.p̃ .+= timestep(problem) .* tableau(method).p.b[2][i] .* cache.Gp[i]
        cache.p̃ .+= timestep(problem) .* tableau(method).p.b[3][i] .* cache.G̅p[i]
    end

    # compute ϕ(q,p)
    functions(problem).ϕ(cache.ϕ̃, solstep.t, cache.q̃, cache.p̃, parameters(solstep))
end


# Compute stages of specialised partitioned additive Runge-Kutta methods for variational systems.
function residual!(
    b::Vector{ST},
    x::Vector{ST},
    solstep::SolutionStepPDAE, 
    problem::HDAEProblem,
    method::HSPARKsecondary, 
    caches::CacheDict) where {ST}

    # get cache for internal stages
    local cache = caches[ST]

    # number of internal stages
    local S = nstages(method)
    local R = pstages(method)
    local D = ndims(problem)

    # compute stages from nonlinear solver solution x
    components!(x, solstep, problem, method, caches)

    # compute b = - [(Y-AV-AU), (Z-AF-AG), Φ]
    for i in 1:S
        for k in 1:D
            b[2*(D*(i-1)+k-1)+1] = - cache.Yi[i][k]
            b[2*(D*(i-1)+k-1)+2] = - cache.Zi[i][k]
            for j in 1:S
                b[2*(D*(i-1)+k-1)+1] += tableau(method).q.a[1][i,j] * cache.Vi[j][k]
                b[2*(D*(i-1)+k-1)+2] += tableau(method).p.a[1][i,j] * cache.Fi[j][k]
            end
            for j in 1:R
                b[2*(D*(i-1)+k-1)+1] += tableau(method).q.a[2][i,j] * cache.Up[j][k]
                b[2*(D*(i-1)+k-1)+1] += tableau(method).q.a[3][i,j] * cache.Λp[j][k]
                b[2*(D*(i-1)+k-1)+2] += tableau(method).p.a[2][i,j] * cache.Gp[j][k]
                b[2*(D*(i-1)+k-1)+2] += tableau(method).p.a[3][i,j] * cache.G̅p[j][k]
            end
        end
    end

    # compute b = - [(Y-AV-AU), (Z-AF-AG), Φ, ωΨ]
    for i in 1:R
        for k in 1:D
            b[2*D*S+4*(D*(i-1)+k-1)+1] = - cache.Yp[i][k]
            b[2*D*S+4*(D*(i-1)+k-1)+2] = - cache.Zp[i][k]
            b[2*D*S+4*(D*(i-1)+k-1)+3] = - cache.Φp[i][k]
            # b[2*D*S+4*(D*(i-1)+k-1)+4] = - cache.Ψp[i][k]
            b[2*D*S+4*(D*(i-1)+k-1)+4] = 0
            for j in 1:S
                b[2*D*S+4*(D*(i-1)+k-1)+1] += tableau(method).q̃.a[1][i,j] * cache.Vi[j][k]
                b[2*D*S+4*(D*(i-1)+k-1)+2] += tableau(method).p̃.a[1][i,j] * cache.Fi[j][k]
            end
            for j in 1:R
                b[2*D*S+4*(D*(i-1)+k-1)+1] += tableau(method).q̃.a[2][i,j] * cache.Up[j][k]
                b[2*D*S+4*(D*(i-1)+k-1)+1] += tableau(method).q̃.a[3][i,j] * cache.Λp[j][k]
                b[2*D*S+4*(D*(i-1)+k-1)+2] += tableau(method).p̃.a[2][i,j] * cache.Gp[j][k]
                b[2*D*S+4*(D*(i-1)+k-1)+2] += tableau(method).p̃.a[3][i,j] * cache.G̅p[j][k]
            end
            for j in 1:R
                b[2*D*S+4*(D*(i-1)+k-1)+4] -= tableau(method).ω[i,j] * cache.Ψp[j][k]
            end
            b[2*D*S+4*(D*(i-1)+k-1)+4] -= tableau(method).ω[i,R+1] * cache.ϕ̃[k]
        end
    end

    # if hasnullvector(method)
    #     for i in 1:R
    #         for k in 1:D
    #             b[2*D*S+4*(D*(i-1)+k-1)+2] -= cache.μ[k] * tableau(method).d[i] / tableau(method).p.b[2][i]
    #         end
    #     end
    
    #     for k in 1:D
    #         b[2*D*S+4*D*R+k] = 0
    #         for i in 1:R
    #             b[2*D*S+4*D*R+k] -= cache.Vp[i][k] * tableau(method).d[i]
    #         end
    #     end
    # end
end


function update_solution!(
    solstep::SolutionStepPDAE{DT,TT}, 
    problem::HDAEProblem,
    method::HSPARKsecondary, 
    caches::CacheDict) where {DT,TT}

    # compute final update
    update!(solstep.q, solstep.q̃, caches[DT].Vi, tableau(method).q.b[1], zero(tableau(method).q.b[1]), timestep(problem))
    update!(solstep.p, solstep.p̃, caches[DT].Fi, tableau(method).p.b[1], zero(tableau(method).p.b[1]), timestep(problem))

    # compute projection
    update!(solstep.q, solstep.q̃, caches[DT].Up, tableau(method).q.b[2], zero(tableau(method).q.b[2]), timestep(problem))
    update!(solstep.q, solstep.q̃, caches[DT].Λp, tableau(method).q.b[3], zero(tableau(method).q.b[3]), timestep(problem))
    update!(solstep.p, solstep.p̃, caches[DT].Gp, tableau(method).p.b[2], zero(tableau(method).p.b[2]), timestep(problem))
    update!(solstep.p, solstep.p̃, caches[DT].G̅p, tableau(method).p.b[3], zero(tableau(method).p.b[3]), timestep(problem))
end
