
"Explicit Runge-Kutta integrator cache."
struct IntegratorCacheERK{DT,D,S} <: ODEIntegratorCache{DT,D}
    Q::Vector{Vector{DT}}
    V::Vector{Vector{DT}}

    function IntegratorCacheERK{DT,D,S}() where {DT,D,S}
        Q = create_internal_stage_vector(DT, D, S)
        V = create_internal_stage_vector(DT, D, S)
        new(Q, V)
    end
end

function Cache{ST}(problem::GeometricProblem, method::ERK; kwargs...) where {ST}
    S = nstages(tableau(method))
    D = ndims(problem)
    IntegratorCacheERK{ST,D,S}(; kwargs...)
end

@inline CacheType(ST, problem::GeometricProblem, method::ERK) = IntegratorCacheERK{ST, ndims(problem), nstages(tableau(method))}


@doc raw"""
Explicit Runge-Kutta integrator solving the system
```math
\begin{aligned}
V_{n,i} &= v (t_i, Q_{n,i}) , &
Q_{n,i} &= q_{n} + h \sum \limits_{j=1}^{s} a_{ij} \, V_{n,j} , &
q_{n+1} &= q_{n} + h \sum \limits_{i=1}^{s} b_{i} \, V_{n,i} .
\end{aligned}
```
"""
const IntegratorERK{DT,TT} = Integrator{<:ODEProblem{DT,TT}, <:ERK}

initmethod(method::ERK) = method
initmethod(method::ERKMethod) = ERK(method)

function Base.show(io::IO, int::IntegratorERK)
    print(io, "\nExplicit Runge-Kutta Integrator with:\n")
    print(io, "   Timestep: $(timestep(int))\n")
    print(io, "   Tableau:  $(description(tableau(int)))\n")
    print(io, "   $(string(tableau(int)))")
    # print(io, reference(tableau(int)))
end


function integrate_step!(solstep::SolutionStepODE{DT,TT}, problem::ODEProblem{DT,TT}, method::ERK, caches::CacheDict, ::NoSolver) where {DT,TT}
    # obtain cache
    local Q::Vector{Vector{DT}} = caches[DT].Q
    local V::Vector{Vector{DT}} = caches[DT].V

    # temporary variables
    local tᵢ::TT
    local yᵢ::DT

    # compute internal stages
    for i in eachstage(method)
        tᵢ = solstep.t̄[1] + timestep(problem) * tableau(method).c[i]
        for k in eachindex(Q[i], V[i])
            yᵢ = 0
            for j in 1:i-1
                yᵢ += tableau(method).a[i,j] * V[j][k]
            end
            Q[i][k] = solstep.q̄[1][k] + timestep(problem) * yᵢ
        end
        functions(problem).v(V[i], tᵢ, Q[i])
    end

    # compute final update
    update!(solstep, V, tableau(method), timestep(problem))
end
