
function Integrators.update_solution!(int::AbstractIntegratorVPRK{DT,TT}, sol::Union{SolutionStepPODE{DT,TT}, SolutionStepPDAE{DT,TT}},
                          cache::IntegratorCacheVPRK{DT}) where {DT,TT}
    update_solution!(sol.q, sol.q̄[1], sol.q̃, cache.V, tableau(int).q.b, tableau(int).q.b̂, timestep(int))
    update_solution!(sol.p, sol.p̄[1], sol.p̃, cache.F, tableau(int).p.b, tableau(int).p.b̂, timestep(int))
end

function project_solution!(int::AbstractIntegratorVPRK{DT,TT}, sol::Union{SolutionStepPODE{DT,TT}, SolutionStepPDAE{DT,TT}}, R::Vector{TT},
                           cache::IntegratorCacheVPRK{DT}) where {DT,TT}
    update_solution!(sol.q, sol.q̃, cache.U, R, timestep(int))
    update_solution!(sol.p, sol.p̃, cache.G, R, timestep(int))
end

function project_solution!(int::AbstractIntegratorVPRK{DT,TT}, sol::Union{SolutionStepPODE{DT,TT}, SolutionStepPDAE{DT,TT}}, RU::Vector{TT}, RG::Vector{TT},
                           cache::IntegratorCacheVPRK{DT}) where {DT,TT}
    update_solution!(sol.q, sol.q̃, cache.U, RU, timestep(int))
    update_solution!(sol.p, sol.p̃, cache.G, RG, timestep(int))
end


function Integrators.initialize!(int::AbstractIntegratorVPRK{DT}, sol::Union{SolutionStepPODE{DT,TT}, SolutionStepPDAE{DT,TT}}) where {DT,TT}
    initialize!(int.iguess, sol.t̄[0], sol.q̄[0], sol.p̄[0], sol.v̄[0], sol.f̄[0],
                            sol.t̄[1], sol.q̄[1], sol.p̄[1], sol.v̄[1], sol.f̄[1])
end


function compute_stages!(x, Q, V, P, F, params::AbstractParametersVPRK)
    # copy x to V
    compute_stages_v_vprk!(x, V, params)

    # compute Q
    compute_stages_q_vprk!(Q, V, params)

    # compute P and F
    compute_stages_p_vprk!(Q, V, P, F, params)
end

function compute_stages!(x, q, p, v, λ, Q, V, U, P, F, G, params::AbstractParametersVPRK)
    # copy x to V
    compute_stages_v_vprk!(x, V, params)

    # compute U, G and p̄
    compute_projection_vprk!(x, q, p, v, λ, Q, V, U, G, params)

    # compute Q
    compute_stages_q_vprk!(Q, V, U, params)

    # compute P and F
    compute_stages_p_vprk!(Q, V, P, F, params)
end

function compute_stages_v_vprk!(x::Vector{ST}, V::Vector{Vector{ST}}, params::AbstractParametersVPRK{IT,DT,TT,D,S}) where {IT,ST,DT,TT,D,S}
    @assert S == length(V)

    # copy x to V
    for i in 1:S
        @assert D == length(V[i])
        for k in 1:D
            V[i][k] = x[D*(i-1)+k]
        end
    end
end

function compute_stages_λ_vprk!(x::Vector{ST}, Λ::Vector{Vector{ST}}, params::AbstractParametersVPRK{IT,DT,TT,D,S}) where {IT,ST,DT,TT,D,S}
    @assert S == length(Λ)

    # copy x to Λ
    for i in 1:S
        @assert D == length(Λ[i])
        for k in 1:D
            Λ[i][k] = x[D*(S+i-1)+k]
        end
    end
end

function compute_stages_q_vprk!(Q::Vector{Vector{ST}}, V::Vector{Vector{ST}}, params::AbstractParametersVPRK{IT,DT,TT,D,S}) where {IT,ST,DT,TT,D,S}
    @assert S == length(Q) == length(V)

    local y1::ST
    local y2::ST

    # compute Q
    for i in 1:S
        @assert D == length(Q[i]) == length(V[i])
        for k in 1:D
            y1 = 0
            y2 = 0
            for j in 1:S
                y1 += params.tab.q.a[i,j] * V[j][k]
                y2 += params.tab.q.â[i,j] * V[j][k]
            end
            Q[i][k] = params.q̄[k] + params.Δt * (y1 + y2)
        end
    end
end


function compute_stages_q_vprk!(Q::Vector{Vector{ST}}, V::Vector{Vector{ST}}, U::Vector{Vector{ST}},
                                params::AbstractParametersVPRK{IT,DT,TT,D,S}) where {IT,ST,DT,TT,D,S}
    @assert S == length(Q) == length(V)
    @assert D == length(U[1]) == length(U[2])

    local y1::ST
    local y2::ST

    # compute Q
    for i in 1:S
        @assert D == length(Q[i]) == length(V[i])
        for k in 1:D
            y1 = 0
            y2 = 0
            for j in 1:S
                y1 += params.tab.q.a[i,j] * V[j][k]
                y2 += params.tab.q.â[i,j] * V[j][k]
            end
            Q[i][k] = params.q̄[k] + params.Δt * (y1 + y2) + params.Δt * params.pparams[:R][1] * U[1][k]
        end
    end
end


function compute_stages_p_vprk!(Q::Vector{Vector{ST}}, V::Vector{Vector{ST}},
                                P::Vector{Vector{ST}}, F::Vector{Vector{ST}},
                                params::AbstractParametersVPRK{IT,DT,TT,D,S}) where {IT,ST,DT,TT,D,S}
    @assert S == length(Q) == length(V) == length(P) == length(F)

    local tᵢ::TT

    # compute P=ϑ(Q) and F=f(Q,V)
    for i in 1:S
        @assert D == length(Q[i]) == length(V[i]) == length(P[i]) == length(F[i])
        tᵢ = params.t̄ + params.Δt * params.tab.q.c[i]
        params.equ.ϑ(P[i], tᵢ, Q[i], V[i])
        params.equ.f(F[i], tᵢ, Q[i], V[i])
    end
end


function compute_rhs_vprk!(b::Vector{ST}, P::Vector{Vector{ST}}, F::Vector{Vector{ST}},
                                params::AbstractParametersVPRK{IT,DT,TT,D,S}) where {IT,ST,DT,TT,D,S}
    @assert S == length(P) == length(F)

    local z1::ST
    local z2::ST

    # compute b = - [(P-p-AF)]
    for i in 1:S
        @assert D == length(P[i]) == length(F[i])
        for k in 1:D
            z1 = 0
            z2 = 0
            for j in 1:S
                z1 += params.tab.p.a[i,j] * F[j][k]
                z2 += params.tab.p.â[i,j] * F[j][k]
            end
            b[D*(i-1)+k] = - (P[i][k] - params.p̄[k]) + params.Δt * (z1 + z2)
        end
    end
end


function compute_rhs_vprk!(b::Vector{ST}, P::Vector{Vector{ST}}, F::Vector{Vector{ST}}, G::Vector{Vector{ST}},
                                    params::AbstractParametersVPRK{IT,DT,TT,D,S}) where {IT,ST,DT,TT,D,S}

    # compute b = - [(P-p-AF)]
    compute_rhs_vprk!(b, P, F, params)

    # compute b += G
    for i in 1:S
        for k in 1:D
            b[D*(i-1)+k] += params.Δt * params.pparams[:R][1] * G[1][k]
        end
    end
end


function compute_rhs_vprk!(b::Vector{ST}, P::Vector{Vector{ST}}, F::Vector{Vector{ST}},
                                          R::Vector{Vector{ST}}, G::Vector{Vector{ST}},
                                          params::AbstractParametersVPRK{IT,DT,TT,D,S}) where {IT,ST,DT,TT,D,S}
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
                z1 += params.tab.p.a[i,j] * (F[j][k] + R[j][k])
                z2 += params.tab.p.â[i,j] * (F[j][k] + R[j][k])
            end
            b[D*(i-1)+k] = - (P[i][k] - params.p̄[k]) + params.Δt * (z1 + z2) + params.Δt * params.pparams[:R][1] * G[1][k]
        end
    end
end


function compute_rhs_vprk_projection_q!(b::Vector{ST}, q::Vector{ST}, V::Vector{Vector{ST}}, U::Vector{Vector{ST}}, offset::Int,
                                    params::AbstractParametersVPRK{IT,DT,TT,D,S}) where {IT,ST,DT,TT,D,S}
    @assert S == length(V)

    local y1::ST
    local y2::ST

    for k in 1:D
        y1 = 0
        y2 = 0
        for j in 1:S
            y1 += params.tab.q.b[j] * V[j][k]
            y2 += params.tab.q.b̂[j] * V[j][k]
        end
        b[offset+k] = - (q[k] - params.q̄[k]) + params.Δt * (y1 + y2) + params.Δt * (params.pparams[:R][1] * U[1][k] + params.pparams[:R][2] * U[2][k])
    end
end


function compute_rhs_vprk_projection_p!(b::Vector{ST}, p::Vector{ST}, F::Vector{Vector{ST}}, G::Vector{Vector{ST}}, offset::Int,
                                    params::AbstractParametersVPRK{IT,DT,TT,D,S}) where {IT,ST,DT,TT,D,S}
    @assert S == length(F)

    local z1::ST
    local z2::ST

    for k in 1:D
        z1 = 0
        z2 = 0
        for j in 1:S
            z1 += params.tab.p.b[j] * F[j][k]
            z2 += params.tab.p.b̂[j] * F[j][k]
        end
        b[offset+k] = - (p[k] - params.p̄[k]) + params.Δt * (z1 + z2) + params.Δt * (params.pparams[:R][1] * G[1][k] + params.pparams[:R][2] * G[2][k])
    end
end


function compute_rhs_vprk_projection_p!(b::Vector{ST}, p::Vector{ST},
                                    F::Vector{Vector{ST}}, R::Vector{Vector{ST}}, G::Vector{Vector{ST}}, offset::Int,
                                    params::AbstractParametersVPRK{IT,DT,TT,D,S}) where {IT,ST,DT,TT,D,S}
    @assert S == length(F)

    local z1::ST
    local z2::ST

    for k in 1:D
        z1 = 0
        z2 = 0
        for j in 1:S
            z1 += params.tab.p.b[j] * (F[j][k] + R[j][k])
            z2 += params.tab.p.b̂[j] * (F[j][k] + R[j][k])
        end
        b[offset+k] = - (p[k] - params.p̄[k]) + params.Δt * (z1 + z2) + params.Δt * (params.pparams[:R][1] * G[1][k] + params.pparams[:R][2] * G[2][k])
    end
end


function compute_rhs_vprk_correction!(b::Vector{ST}, V::Vector{Vector{ST}},
                params::AbstractParametersVPRK{IT,DT,TT,D,S}) where {IT,ST,DT,TT,D,S}

    @assert S == length(V)

    local sl::Int = div(S+1, 2)
    local μ = zeros(ST,D)

    if !isnothing(params.nullvec)
        # compute μ
        for k in 1:D
            μ[k] = params.tab.p.b[sl] / params.nullvec[sl] * b[D*(sl-1)+k]
        end

        # replace equation for Pₗ with constraint on V
        for k in 1:D
            b[D*(sl-1)+k] = 0
            for i in 1:S
                b[D*(sl-1)+k] += V[i][k] * params.nullvec[i]
            end
        end

        # modify P₁, ..., Pₛ except for Pₗ
        for i in 1:S
            if i ≠ sl
                z = params.nullvec[i] / params.tab.p.b[i]
                for k in 1:D
                    b[D*(i-1)+k] -= z * μ[k]
                end
            end
        end
    end
end
