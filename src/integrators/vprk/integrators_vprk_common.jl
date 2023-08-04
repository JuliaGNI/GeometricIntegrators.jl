
function update_solution!(int::GeometricIntegratorVPRK{DT,TT}, sol::Union{SolutionStepPODE{DT,TT}, SolutionStepPDAE{DT,TT}},
                          cache::IntegratorCacheVPRK{DT}) where {DT,TT}
    update!(sol, cache.V, cache.F, tableau(int), timestep(int))
end

function project_solution!(int::GeometricIntegratorVPRK{DT,TT}, sol::Union{SolutionStepPODE{DT,TT}, SolutionStepPDAE{DT,TT}}, R::Vector{TT},
                           cache::IntegratorCacheVPRK{DT}) where {DT,TT}
    update!(sol.q, sol.q̃, cache.U, R, timestep(int))
    update!(sol.p, sol.p̃, cache.G, R, timestep(int))
end

function project_solution!(int::GeometricIntegratorVPRK{DT,TT}, sol::Union{SolutionStepPODE{DT,TT}, SolutionStepPDAE{DT,TT}}, RU::Vector{TT}, RG::Vector{TT},
                           cache::IntegratorCacheVPRK{DT}) where {DT,TT}
    update!(sol.q, sol.q̃, cache.U, RU, timestep(int))
    update!(sol.p, sol.p̃, cache.G, RG, timestep(int))
end




function components!(
    x::AbstractVector{ST},
    solstep::SolutionStepPODE,
    problem::VPRKProblem,
    method::ProjectedVPRK,
    caches::CacheDict) where {ST}

    # get cache
    local cache = caches[ST]

    # copy x to V
    compute_stages_v!(x, cache.V, solstep, problem, method)

    # compute U, G and p̄
    compute_projection_vprk!(x, cache.q, cache.p, cache.v, cache.λ, cache.Q, cache.V, cache.U, cache.G, solstep, problem, method)

    # compute Q
    compute_stages_q!(x, cache.Q, cache.V, cache.U, solstep, problem, method)

    # compute P and F
    compute_stages_p_vprk!(cache.Q, cache.V, cache.P, cache.F, solstep, problem, method)
end

function compute_stages_λ_vprk!(
    x::AbstractVector{ST}, 
    Λ::Vector{Vector{ST}}, 
    solstep::SolutionStepPODE, 
    problem::VPRKProblem,
    method::VPRKMethod) where {ST}

    local S = nstages(method)
    local D = ndims(problem)

    @assert S == length(Λ)

    # copy x to Λ
    for i in 1:S
        @assert D == length(Λ[i])
        for k in 1:D
            Λ[i][k] = x[D*(S+i-1)+k]
        end
    end
end


function compute_stages_q!(
    x::AbstractVector{ST},
    solstep::SolutionStepPODE, 
    problem::VPRKProblem,
    method::ProjectedVPRK,
    caches::CacheDict) where {ST}

    local S = nstages(method)
    local D = ndims(problem)
    local Q = caches[ST].Q
    local U = caches[ST].U

    # compute Q
    compute_stages_q!(x, solstep, problem, parent(method), caches)

    # add projection to Q
    for i in 1:S
        for k in 1:D
            Q[i][k] += timestep(problem) * params.pparams[:R][1] * U[1][k]
        end
    end
end





function compute_rhs_vprk!(b::Vector{ST}, P::Vector{Vector{ST}}, F::Vector{Vector{ST}}, G::Vector{Vector{ST}}, 
    solstep::SolutionStepPODE, 
    problem::VPRKProblem,
    method::VPRKMethod) where {ST}

    local S = nstages(method)
    local D = ndims(problem)

    # compute b = - [(P-p-AF)]
    compute_rhs_vprk!(b, P, F, solstep, problem, method)

    # compute b += G
    for i in 1:S
        for k in 1:D
            b[D*(i-1)+k] += timestep(problem) * params.pparams[:R][1] * G[1][k]
        end
    end
end


function compute_rhs_vprk!(b::Vector{ST}, p̄::Vector{ST},
                                          P::Vector{Vector{ST}}, F::Vector{Vector{ST}},
                                          R::Vector{Vector{ST}}, G::Vector{Vector{ST}}, 
                                          solstep::SolutionStepPODE, 
                                          problem::VPRKProblem,
                                          method::VPRKMethod) where {ST}
                                      
    local S = nstages(method)
    local D = ndims(problem)

    @assert S == length(P) == length(F)
    @assert S == length(R) == length(G)

    local z1::ST
    local z2::ST

    # compute b = - [(P-G-AF)]
    for i in 1:S
        @assert D == length(P[i]) == length(F[i])
        for k in 1:D
            z1 = 0
            z2 = 0
            for j in 1:S
                z1 += tableau(method).p.a[i,j] * (F[j][k] + R[j][k])
                z2 += tableau(method).p.â[i,j] * (F[j][k] + R[j][k])
            end
            b[D*(i-1)+k] = - ( P[i][k] - p̄[k] ) + timestep(problem) * (z1 + z2) + timestep(problem) * params.pparams[:R][1] * G[1][k]
        end
    end
end


function compute_rhs_vprk_projection_q!(b::Vector{ST}, q::Vector{ST}, V::Vector{Vector{ST}}, U::Vector{Vector{ST}}, offset::Int, 
    solstep::SolutionStepPODE, 
    problem::VPRKProblem,
    method::VPRKMethod) where {ST}

    local S = nstages(method)
    local D = ndims(problem)

    @assert S == length(V)

    local y1::ST
    local y2::ST

    for k in 1:D
        y1 = 0
        y2 = 0
        for j in 1:S
            y1 += tableau(method).q.b[j] * V[j][k]
            y2 += tableau(method).q.b̂[j] * V[j][k]
        end
        b[offset+k] = - ( q[k] - solstep.q̄[k] ) + timestep(problem) * (y1 + y2) + timestep(problem) * (params.pparams[:R][1] * U[1][k] + params.pparams[:R][2] * U[2][k])
    end
end


function compute_rhs_vprk_projection_p!(b::Vector{ST}, p::Vector{ST}, F::Vector{Vector{ST}}, G::Vector{Vector{ST}}, offset::Int, 
    solstep::SolutionStepPODE, 
    problem::VPRKProblem,
    method::VPRKMethod) where {ST}

    local S = nstages(method)
    local D = ndims(problem)

    @assert S == length(F)

    local z1::ST
    local z2::ST

    for k in 1:D
        z1 = 0
        z2 = 0
        for j in 1:S
            z1 += tableau(method).p.b[j] * F[j][k]
            z2 += tableau(method).p.b̂[j] * F[j][k]
        end
        b[offset+k] = - ( p[k] - solstep.p̄[k] ) + timestep(problem) * (z1 + z2) + timestep(problem) * (params.pparams[:R][1] * G[1][k] + params.pparams[:R][2] * G[2][k])
    end
end


function compute_rhs_vprk_projection_p!(b::Vector{ST}, offset::Int, 
    solstep::SolutionStepPODE, 
    problem::VPRKProblem,
    method::VPRKMethod) where {ST}

    local S = nstages(method)
    local D = ndims(problem)

    @assert S == length(F)

    local z1::ST
    local z2::ST

    for k in 1:D
        z1 = 0
        z2 = 0
        for j in 1:S
            z1 += tableau(method).p.b[j] * (F[j][k] + R[j][k])
            z2 += tableau(method).p.b̂[j] * (F[j][k] + R[j][k])
        end
        b[offset+k] = - ( p[k] - solstep.p̄[k] ) + timestep(problem) * (z1 + z2) + timestep(problem) * (params.pparams[:R][1] * G[1][k] + params.pparams[:R][2] * G[2][k])
    end
end

