
@doc raw"""
Variational Partitioned Runge-Kutta Integrator.

```math
\begin{aligned}
P_{n,i} &= \dfrac{\partial L}{\partial v} (Q_{n,i}, V_{n,i}) , &
Q_{n,i} &= q_{n} + h \sum \limits_{j=1}^{s} a_{ij} \, V_{n,j} , &
q_{n+1} &= q_{n} + h \sum \limits_{i=1}^{s} b_{i} \, V_{n,i} , \\
F_{n,i} &= \dfrac{\partial L}{\partial q} (Q_{n,i}, V_{n,i}) , &
P_{n,i} &= p_{n} + h \sum \limits_{i=1}^{s} \bar{a}_{ij} \, F_{n,j} - d_i \lambda , &
p_{n+1} &= p_{n} + h \sum \limits_{i=1}^{s} \bar{b}_{i} \, F_{n,i} , \\
&&
0 &= \sum \limits_{i=1}^{s} d_i V_i , &&
\end{aligned}
```
satisfying the symplecticity conditions
```math
\begin{aligned}
b_{i} \bar{a}_{ij} + b_{j} a_{ji} &= b_{i} b_{j} , &
\bar{b}_i &= b_i .
\end{aligned}
```
"""
const IntegratorVPRK{DT,TT} = Integrator{<:Union{IODEProblem{DT,TT},LODEProblem{DT,TT}}, <:VPRKMethod}

description(::IntegratorVPRK) = "Variational Partitioned Runge-Kutta Integrator"

default_solver(::VPRKMethod) = Newton()
default_iguess(::VPRKMethod) = HermiteExtrapolation()

initmethod(method::VPRKMethod) = Methods.VPRK(method)

solversize(problem::VPRKProblem, method::VPRKMethod) =
    ndims(problem) * nstages(method)

function initsolver(::NewtonMethod, method::VPRKMethod, caches::CacheDict{<:VPRKProblem})
    # create nonlinear solver
    NewtonSolver(zero(nlsolution(caches)), zero(nlsolution(caches)); linesearch = Backtracking(), config = Options(min_iterations = 1, x_abstol = 8eps(), f_abstol = 8eps()))
end


# function Integrators.get_internal_variables(int::IntegratorVPRK{DT,TT,D,S}) where {DT, TT, D, S}
#     Q = create_internal_stage_vector(DT, D, S)
#     P = create_internal_stage_vector(DT, D, S)
#     V = create_internal_stage_vector(DT, D, S)
#     F = create_internal_stage_vector(DT, D, S)

#     solver = get_solver_status(int.solver)

#     (Q=Q, P=P, V=V, F=F, solver=solver)
# end


function initial_guess!(
    solstep::SolutionStepPODE{DT}, 
    problem::VPRKProblem,
    method::VPRKMethod, 
    caches::CacheDict, 
    ::NonlinearSolver, 
    iguess::Union{InitialGuess,Extrapolation}) where {DT}

    cache = caches[DT]

    for i in eachstage(method)
        # TODO: Check initial guess !
        initialguess!(solstep.t̄[1] + timestep(problem) * tableau(method).q.c[i], cache.Q[i], cache.V[i], solstep, problem, iguess)
        for k in 1:ndims(problem)
            cache.x[ndims(problem)*(i-1)+k] = 0#solstep.v̄[1][k]#cache.V[i][k]
        end
        # println("  t = $(solstep.t̄[1] + timestep(problem) * tableau(method).q.c[i]),",
        #         "  q̄ = $(solstep.q̄[1]), v̄ = $(solstep.v̄[1]), ",
        #         "  q = $(cache.Q[i]), v = $(cache.V[i])")
    end
end


function components_v!(
    x::AbstractVector{ST},
    solstep::SolutionStepPODE, 
    problem::VPRKProblem,
    method::VPRKMethod,
    caches::CacheDict) where {ST}

    local S = nstages(method)
    local D = ndims(problem)
    local V = caches[ST].V

    # copy x to V
    for i in 1:S
        @assert D == length(V[i])
        for k in 1:D
            V[i][k] = x[D*(i-1)+k]
        end
    end
end

function components_q!(
    x::AbstractVector{ST},
    solstep::SolutionStepPODE, 
    problem::VPRKProblem,
    method::VPRKMethod,
    caches::CacheDict) where {ST}

    local S = nstages(method)
    local D = ndims(problem)
    local q̄ = caches[ST].q̄
    local Q = caches[ST].Q
    local V = caches[ST].V

    local y1::ST
    local y2::ST

    # compute Q
    for i in 1:S
        @assert D == length(Q[i]) == length(V[i])
        for k in 1:D
            y1 = 0
            y2 = 0
            for j in 1:S
                y1 += tableau(method).q.a[i,j] * V[j][k]
                y2 += tableau(method).q.â[i,j] * V[j][k]
            end
            Q[i][k] = q̄[k] + timestep(problem) * (y1 + y2)
        end
        # if ST == Float64
        #     println("    Q[$i] = $(Q[i]),  V[$i] = $(V[i])")
        # end
    end
end

function components_p!(
    x::AbstractVector{ST},
    solstep::SolutionStepPODE, 
    problem::VPRKProblem,
    method::VPRKMethod,
    caches::CacheDict) where {ST}
                            
    local S = nstages(method)
    local Q = caches[ST].Q
    local P = caches[ST].P
    local V = caches[ST].V
    local F = caches[ST].F
    
    local tᵢ::timetype(problem)

    # compute P=ϑ(Q,V) and F=f(Q,V)
    for i in 1:S
        tᵢ = solstep.t̄[1] + timestep(problem) * tableau(method).p.c[i]
        functions(problem).ϑ(P[i], tᵢ, Q[i], V[i])
        functions(problem).f(F[i], tᵢ, Q[i], V[i])
        # if ST == Float64
        #     println("    P[$i] = $(P[i]),  F[$i] = $(F[i])")
        # end
    end
end

function components!(x::AbstractVector{ST}, int::IntegratorVPRK) where {ST}
    # copy x to V
    components_v!(x, solstep(int), problem(int), method(int), caches(int))

    # compute Q
    components_q!(x, solstep(int), problem(int), method(int), caches(int))

    # compute P and F
    components_p!(x, solstep(int), problem(int), method(int), caches(int))
end


function residual!(
    b::AbstractVector{ST},
    solstep::SolutionStepPODE, 
    problem::VPRKProblem,
    method::VPRKMethod,
    caches::CacheDict) where {ST}

    local S = nstages(method)
    local D = ndims(problem)
    local p̄ = caches[ST].p̄
    local P = caches[ST].P
    local F = caches[ST].F

    local z1::ST
    local z2::ST

    # compute b = - [(P-p-AF)]
    for i in 1:S
        for k in 1:D
            z1 = 0
            z2 = 0
            for j in 1:S
                z1 += tableau(method).p.a[i,j] * F[j][k]
                z2 += tableau(method).p.â[i,j] * F[j][k]
            end
            b[D*(i-1)+k] = - ( P[i][k] - p̄[k] ) + timestep(problem) * (z1 + z2)
        end
    end
end


function residual_correction!(
    b::AbstractVector{ST},
    problem::VPRKProblem,
    method::VPRKMethod, 
    caches::CacheDict) where {ST}

    local S = nstages(method)
    local D = ndims(problem)
    local V = caches[ST].V
    local μ = caches[ST].μ

    local sl::Int = div(S+1, 2)

    if hasnullvector(method)
        # compute μ
        for k in 1:D
            μ[k] = tableau(method).p.b[sl] / nullvector(method)[sl] * b[D*(sl-1)+k]
        end

        # replace equation for Pₗ with constraint on V
        for k in 1:D
            b[D*(sl-1)+k] = 0
            for i in 1:S
                b[D*(sl-1)+k] += V[i][k] * nullvector(method)[i]
            end
        end

        # modify P₁, ..., Pₛ except for Pₗ
        for i in 1:S
            if i ≠ sl
                z = nullvector(method)[i] / tableau(method).p.b[i]
                for k in 1:D
                    b[D*(i-1)+k] -= z * μ[k]
                end
            end
        end
    end
end


# Compute stages of variational partitioned Runge-Kutta methods.
function residual!(b::AbstractVector, int::IntegratorVPRK)
    residual!(b, solstep(int), problem(int), method(int), caches(int))
    residual_correction!(b, problem(int), method(int), caches(int))
end


# Compute stages of Variational Partitioned Runge-Kutta methods.
function residual!(b::AbstractVector{ST}, x::AbstractVector{ST}, int::IntegratorVPRK) where {ST}
    @assert axes(x) == axes(b)

    # copy previous solution from solstep to cache
    reset!(cache(int, ST), current(solstep(int))...)

    # compute stages from nonlinear solver solution x
    components!(x, int)

    # compute residual vector
    residual!(b, int)
end


function update!(x::AbstractVector{DT}, int::IntegratorVPRK) where {DT}
    # copy previous solution from solstep to cache
    reset!(cache(int, DT), current(solstep(int))...)

    # compute vector field at internal stages
    components!(x, int)

    # compute final update
    update!(solstep(int), cache(int, DT).V, cache(int, DT).F, tableau(int), timestep(int))
end


function integrate_step!(int::IntegratorVPRK)
    # call nonlinear solver
    solve!(nlsolution(int), (b,x) -> residual!(b, x, int), solver(int))

    # check_jacobian(solver(int))
    # print_jacobian(solver(int))

    # print solver status
    # println(status(solver(int)))

    # check if solution contains NaNs or error bounds are violated
    # println(meets_stopping_criteria(status(solver(int))))

    # compute final update
    update!(nlsolution(int), int)

    # update vector field for initial guess
    update_vector_fields!(solstep(int), problem(int))

    # copy internal stage variables
    # solstep(int).internal[:Q] .= cache(int).Q
    # solstep(int).internal[:P] .= cache(int).P
    # solstep(int).internal[:V] .= cache(int).V
    # solstep(int).internal[:F] .= cache(int).F

    # copy solver status
    # get_solver_status!(solver(int), solstep(int).internal[:solver])
end
