
abstract type AbstractParametersVPRK{DT,TT,ET,D,S} <: Parameters{DT,TT} end
abstract type AbstractIntegratorVPRK{DT,TT} <: DeterministicIntegrator{DT,TT} end

equation(integrator::AbstractIntegratorVPRK) = integrator.params.equ
timestep(integrator::AbstractIntegratorVPRK) = integrator.params.Δt
tableau(integrator::AbstractIntegratorVPRK) = integrator.params.tab
dims(integrator::AbstractIntegratorVPRK) = integrator.params.equ.d


struct NonlinearFunctionCacheVPRK{ST}
    Q::Vector{Vector{ST}}
    P::Vector{Vector{ST}}
    V::Vector{Vector{ST}}
    F::Vector{Vector{ST}}

    v::Vector{ST}
    y::Vector{ST}
    z::Vector{ST}

    q::Vector{ST}
    p::Vector{ST}

    function NonlinearFunctionCacheVPRK{ST}(D,S) where {ST}
        # create internal stage vectors
        Q = create_internal_stage_vector(ST, D, S)
        P = create_internal_stage_vector(ST, D, S)
        V = create_internal_stage_vector(ST, D, S)
        F = create_internal_stage_vector(ST, D, S)

        # create update vectors
        v = zeros(ST,D)
        y = zeros(ST,D)
        z = zeros(ST,D)

        # create solution vector
        q = zeros(ST,D)
        p = zeros(ST,D)

        new(Q, V, P, F, v, y, z, q, p)
    end
end

struct NonlinearFunctionCacheVPRKprojection{ST}
    q::Vector{ST}
    p::Vector{ST}
    λ::Vector{ST}

    Λ::Vector{Vector{ST}}
    Φ::Vector{Vector{ST}}

    U::Vector{Vector{ST}}
    G::Vector{Vector{ST}}
    R::Vector{Vector{ST}}

    u::Vector{ST}
    g::Vector{ST}

    function NonlinearFunctionCacheVPRKprojection{ST}(D,S) where {ST}
        # create projected solution vectors
        q = zeros(ST,D)
        p = zeros(ST,D)
        λ = zeros(ST,D)

        Λ = create_internal_stage_vector(ST, D, S)
        Φ = create_internal_stage_vector(ST, D, S)

        # create projection vectors
        U = create_internal_stage_vector(ST, D, 2)
        G = create_internal_stage_vector(ST, D, 2)
        R = create_internal_stage_vector(ST, D, S)

        # create update vectors
        u = zeros(ST,D)
        g = zeros(ST,D)

        new(q, p, λ, Λ, Φ, U, G, R, u, g)
    end
end


function compute_stages!(x, Q, V, P, F, params::AbstractParametersVPRK)
    # copy x to V
    compute_stages_v_vprk!(x, V, params)

    # compute Q
    compute_stages_q_vprk!(Q, V, params)

    # compute P and F
    compute_stages_p_vprk!(Q, V, P, F, params)
end

function compute_stages!(x, q, p, λ, Q, V, U, P, F, G, params::AbstractParametersVPRK)
    # copy x to V
    compute_stages_v_vprk!(x, V, params)

    # compute U, G and p̅
    compute_projection_vprk!(x, q, p, λ, V, U, G, params)

    # compute Q
    compute_stages_q_vprk!(Q, V, U, params)

    # compute P and F
    compute_stages_p_vprk!(Q, V, P, F, params)
end

function compute_stages_v_vprk!(x::Vector{ST}, V::Vector{Vector{ST}}, params::AbstractParametersVPRK{DT,TT,ET,D,S}) where {ST,DT,TT,ET,D,S}
    @assert S == length(V)

    # copy x to V
    for i in 1:S
        @assert D == length(V[i])
        for k in 1:D
            V[i][k] = x[D*(i-1)+k]
        end
    end
end

function compute_stages_λ_vprk!(x::Vector{ST}, Λ::Vector{Vector{ST}}, params::AbstractParametersVPRK{DT,TT,ET,D,S}) where {ST,DT,TT,ET,D,S}
    @assert S == length(Λ)

    # copy x to Λ
    for i in 1:S
        @assert D == length(Λ[i])
        for k in 1:D
            Λ[i][k] = x[D*(S+i-1)+k]
        end
    end
end

function compute_stages_q_vprk!(Q::Vector{Vector{ST}}, V::Vector{Vector{ST}}, params::AbstractParametersVPRK{DT,TT,ET,D,S}) where {ST,DT,TT,ET,D,S}
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
            Q[i][k] = params.q[k] + params.Δt * (y1 + y2)
        end
    end
end


function compute_stages_q_vprk!(Q::Matrix{ST}, V::Matrix{ST}, U::Matrix{ST}, params::AbstractParametersVPRK{DT,TT,ET,D,S}) where {ST,DT,TT,ET,D,S}
    @assert D == size(Q,1) == size(V,1) == size(U,1)
    @assert S == size(Q,2) == size(V,2)

    local y1::ST
    local y2::ST

    # compute Q
    for i in 1:S
        for k in 1:D
            y1 = 0
            y2 = 0
            for j in 1:S
                y1 += params.tab.q.a[i,j] * V[k,j]
                y2 += params.tab.q.â[i,j] * V[k,j]
            end
            Q[k,i] = params.q[k] + params.Δt * (y1 + y2 + params.R[1] * U[k,1])
        end
    end
end


function compute_stages_p_vprk!(Q::Vector{Vector{ST}}, V::Vector{Vector{ST}},
                                P::Vector{Vector{ST}}, F::Vector{Vector{ST}},
                                params::AbstractParametersVPRK{DT,TT,ET,D,S}) where {ST,DT,TT,ET,D,S}
    @assert S == length(Q) == length(V) == length(P) == length(F)

    local tᵢ::TT

    # compute P=ϑ(Q) and F=f(Q,V)
    for i in 1:S
        @assert D == length(Q[i]) == length(V[i]) == length(P[i]) == length(F[i])
        tᵢ = params.t + params.Δt * params.tab.q.c[i]
        params.equ.ϑ(tᵢ, Q[i], V[i], P[i])
        params.equ.f(tᵢ, Q[i], V[i], F[i])
    end
end


function compute_rhs_vprk!(b::Vector{ST}, P::Vector{Vector{ST}}, F::Vector{Vector{ST}},
                                params::AbstractParametersVPRK{DT,TT,ET,D,S}) where {ST,DT,TT,ET,D,S}
    @assert S == length(P) == length(F)

    local z1::ST
    local z2::ST

    # compute b = - [(P-AF)]
    for i in 1:S
        @assert D == length(P[i]) == length(F[i])
        for k in 1:D
            z1 = 0
            z2 = 0
            for j in 1:S
                z1 += params.tab.p.a[i,j] * F[j][k]
                z2 += params.tab.p.â[i,j] * F[j][k]
            end
            b[D*(i-1)+k] = - (P[i][k] - params.p[k]) + params.Δt * (z1 + z2)
        end
    end
end


function compute_rhs_vprk!(b::Vector{ST}, P::Vector{Vector{ST}}, F::Vector{Vector{ST}}, G::Vector{Vector{ST}},
                                    params::AbstractParametersVPRK{DT,TT,ET,D,S}) where {ST,DT,TT,ET,D,S}
    @assert S == length(P) == length(F)

    local z1::ST
    local z2::ST

    # compute b = - [(P-G-AF)]
    for i in 1:S
        @assert D == length(P[i]) == length(F[i])
        for k in 1:D
            z1 = 0
            z2 = 0
            for j in 1:S
                z1 += params.tab.p.a[i,j] * F[j][k]
                z2 += params.tab.p.â[i,j] * F[j][k]
            end
            b[D*(i-1)+k] = - (P[i][k] - params.p[k]) + params.Δt * (z1 + z2 + params.R[1] * G[1][k])
        end
    end
end


function compute_rhs_vprk!(b::Vector{ST}, P::Vector{Vector{ST}}, F::Vector{Vector{ST}},
                                          R::Vector{Vector{ST}}, G::Vector{Vector{ST}},
                                          params::AbstractParametersVPRK{DT,TT,ET,D,S}) where {ST,DT,TT,ET,D,S}
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
            b[D*(i-1)+k] = - (P[i][k] - params.p[k]) + params.Δt * (z1 + z2 + params.R[1] * G[1][k])
        end
    end
end


function compute_rhs_vprk_projection_q!(b::Vector{ST}, q::Vector{ST}, V::Vector{Vector{ST}}, U::Vector{Vector{ST}}, offset::Int,
                                    params::AbstractParametersVPRK{DT,TT,ET,D,S}) where {ST,DT,TT,ET,D,S}
    @assert S == length(V)

    local y1::ST
    local y2::ST

    for k in 1:D
        @assert D == length(V[i])
        y1 = 0
        y2 = 0
        for j in 1:S
            y1 += params.tab.q.b[j] * V[j][k]
            y2 += params.tab.q.b̂[j] * V[j][k]
        end
        b[offset+k] = - (q[k] - params.q[k]) + params.Δt * (y1 + y2) + params.Δt * (params.R[1] * U[1][k] + params.R[2] * U[2][k])
    end
end


function compute_rhs_vprk_projection_p!(b::Vector{ST}, p::Vector{ST}, F::Vector{Vector{ST}}, G::Vector{Vector{ST}}, offset::Int,
                                    params::AbstractParametersVPRK{DT,TT,ET,D,S}) where {ST,DT,TT,ET,D,S}
    @assert S == length(F)

    local z1::ST
    local z2::ST

    for k in 1:D
        @assert D == length(F[i])
        z1 = 0
        z2 = 0
        for j in 1:S
            z1 += params.tab.p.b[j] * F[j][k]
            z2 += params.tab.p.b̂[j] * F[j][k]
        end
        b[offset+k] = - (p[k] - params.p[k]) + params.Δt * (z1 + z2) + params.Δt * (params.R[1] * G[1][k] + params.R[2] * G[2][k])
    end
end


function compute_rhs_vprk_projection_p!(b::Vector{ST}, p::Vector{ST},
                                    F::Vector{Vector{ST}}, R::Vector{Vector{ST}}, G::Vector{Vector{ST}}, offset::Int,
                                    params::AbstractParametersVPRK{DT,TT,ET,D,S}) where {ST,DT,TT,ET,D,S}
    @assert S == length(F)

    local z1::ST
    local z2::ST

    for k in 1:D
        @assert D == length(F[i])
        z1 = 0
        z2 = 0
        for j in 1:S
            z1 += params.tab.p.b[j] * (F[j][k] + R[j][k])
            z2 += params.tab.p.b̂[j] * (F[j][k] + R[j][k])
        end
        b[offset+k] = - (p[k] - params.p[k]) + params.Δt * (z1 + z2) + params.Δt * (params.R[1] * G[1][k] + params.R[2] * G[2][k])
    end
end


@generated function compute_rhs_vprk_correction!(b::Vector{ST}, V::Vector{Vector{ST}},
                                    params::AbstractParametersVPRK{DT,TT,ET,D,S}) where {ST,DT,TT,ET,D,S}
    μ = zeros(ST,D)

    quote
        @assert S == length(V)

        local sl::Int = div(S+1, 2)

        if isdefined(params.tab, :d)
            # compute μ
            for k in 1:D
                $μ[k] = params.tab.p.b[sl] / params.tab.d[sl] * b[D*(sl-1)+k]
            end

            # replace equation for Pₗ with constraint on V
            for k in 1:D
                b[D*(sl-1)+k] = 0
                for i in 1:S
                    b[D*(sl-1)+k] += V[i][k] * params.tab.d[i]
                end
            end

            # modify P₁, ..., Pₛ except for Pₗ
            for i in 1:S
                if i ≠ sl
                    z = params.tab.d[i] / params.tab.p.b[i]
                    for k in 1:D
                        b[D*(i-1)+k] -= z * $μ[k]
                    end
                end
            end
        end
    end
end


function update_solution!(int::AbstractIntegratorVPRK{DT,TT}, cache::NonlinearFunctionCacheVPRK{DT}) where {DT,TT}
    update_solution!(int.q, cache.V, params.tab.q.b, params.tab.q.b̂, int.params.Δt)
    update_solution!(int.p, cache.F, params.tab.p.b, params.tab.p.b̂, int.params.Δt)
end


function project_solution!(int::AbstractIntegratorVPRK{DT,TT}, cache::NonlinearFunctionCacheVPRKprojection{DT}, R::Vector{TT}) where {DT,TT}
    update_solution!(int.q, cache.U, R, int.params.Δt)
    update_solution!(int.p, cache.G, R, int.params.Δt)
end


function cut_periodic_solution!(int::AbstractIntegratorVPRK)
    cut_periodic_solution!(int.q, int.params.equ.periodicity)
end
