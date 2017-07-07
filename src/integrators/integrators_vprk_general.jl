
abstract type AbstractNonlinearFunctionParametersVPRK{DT,TT,ΑT,FT,D,S} <: NonlinearFunctionParameters{DT} end
abstract type AbstractIntegratorVPRK{DT,TT} <: Integrator{DT,TT} end

struct NonlinearFunctionCacheVPRK{ST}
    Q::Matrix{ST}
    V::Matrix{ST}
    P::Matrix{ST}
    F::Matrix{ST}

    v::Array{ST,1}
    y::Array{ST,1}
    z::Array{ST,1}

    function NonlinearFunctionCacheVPRK{ST}(D,S) where {ST}
        # create internal stage vectors
        Q = zeros(ST,D,S)
        V = zeros(ST,D,S)
        P = zeros(ST,D,S)
        F = zeros(ST,D,S)

        # create update vectors
        v = zeros(ST,D)
        y = zeros(ST,D)
        z = zeros(ST,D)

        new(Q, V, P, F, v, y, z)
    end
end

struct NonlinearFunctionCacheVPRKprojection{ST}
    q̅::Vector{ST}
    p̅::Vector{ST}
    λ::Vector{ST}

    Λ::Matrix{ST}
    Φ::Matrix{ST}

    U::Array{ST,2}
    G::Array{ST,2}
    R::Array{ST,2}

    u::Array{ST,1}
    g::Array{ST,1}

    function NonlinearFunctionCacheVPRKprojection{ST}(D,S) where {ST}
        # create projected solution vectors
        q̅ = zeros(ST,D)
        p̅ = zeros(ST,D)
        λ = zeros(ST,D)

        Λ = zeros(ST,D,S)
        Φ = zeros(ST,D,S)

        # create projection vectors
        U = zeros(ST,D,2)
        G = zeros(ST,D,2)
        R = zeros(ST,D,S)

        # create update vectors
        u = zeros(ST,D)
        g = zeros(ST,D)

        new(q̅, p̅, λ, Λ, Φ, U, G, R, u, g)
    end
end


function compute_stages_vprk!(x, Q, V, P, F, params)
    # copy x to V
    compute_stages_v_vprk!(x, V, params)

    # compute Q
    compute_stages_q_vprk!(Q, V, params)

    # compute P and F
    compute_stages_p_vprk!(Q, V, P, F, params)
end

function compute_stages_vprk!(x, q̅, p̅, λ, Q, V, U, P, F, G, params)
    # copy x to V
    compute_stages_v_vprk!(x, V, params)

    # compute U, G and p̅
    compute_projection_vprk!(x, q̅, p̅, λ, V, U, G, params)

    # compute Q
    compute_stages_q_vprk!(Q, V, U, params)

    # compute P and F
    compute_stages_p_vprk!(Q, V, P, F, params)
end

function compute_stages_v_vprk!(x::Vector{ST}, V::Matrix{ST}, params::AbstractNonlinearFunctionParametersVPRK{DT,TT,AT,FT,D,S}) where {ST,DT,TT,AT,FT,D,S}
    @assert D == size(V,1)
    @assert S == size(V,2)

    # copy x to V
    for i in 1:S
        for k in 1:D
            V[k,i] = x[D*(i-1)+k]
        end
    end
end

function compute_stages_λ_vprk!(x::Vector{ST}, Λ::Matrix{ST}, params::AbstractNonlinearFunctionParametersVPRK{DT,TT,AT,FT,D,S}) where {ST,DT,TT,AT,FT,D,S}
    @assert D == size(Λ,1)
    @assert S == size(Λ,2)

    # copy x to V
    for i in 1:S
        for k in 1:D
            Λ[k,i] = x[D*(S+i-1)+k]
        end
    end
end

function compute_stages_q_vprk!(Q::Matrix{ST}, V::Matrix{ST}, params::AbstractNonlinearFunctionParametersVPRK{DT,TT,AT,FT,D,S}) where {ST,DT,TT,AT,FT,D,S}
    @assert D == size(Q,1) == size(V,1)
    @assert S == size(Q,2) == size(V,2)

    local y1::ST
    local y2::ST

    # compute Q
    for i in 1:S
        for k in 1:D
            y1 = 0
            y2 = 0
            for j in 1:S
                y1 += params.t_q.a[i,j] * V[k,j]
                y2 += params.t_q.â[i,j] * V[k,j]
            end
            Q[k,i] = params.q[k] + params.Δt * (y1 + y2)
        end
    end
end


function compute_stages_q_vprk!(Q::Matrix{ST}, V::Matrix{ST}, U::Matrix{ST}, params::AbstractNonlinearFunctionParametersVPRK{DT,TT,AT,FT,D,S}) where {ST,DT,TT,AT,FT,D,S}
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
                y1 += params.t_q.a[i,j] * V[k,j]
                y2 += params.t_q.â[i,j] * V[k,j]
            end
            Q[k,i] = params.q[k] + params.Δt * (y1 + y2 + params.R[1] * U[k,1])
        end
    end
end


@generated function compute_stages_p_vprk!(Q::Matrix{ST}, V::Matrix{ST}, P::Matrix{ST}, F::Matrix{ST}, params::AbstractNonlinearFunctionParametersVPRK{DT,TT,AT,FT,D,S}) where {ST,DT,TT,AT,FT,D,S}
    # create temporary vectors
    tQ = zeros(ST,D)
    tV = zeros(ST,D)
    tP = zeros(ST,D)
    tF = zeros(ST,D)

    compute_stages_vprk = quote
        @assert D == size(Q,1) == size(V,1) == size(P,1) == size(F,1)
        @assert S == size(Q,2) == size(V,2) == size(P,2) == size(F,2)

        local tᵢ::TT

        # compute P=α(Q) and F=f(Q)
        for i in 1:S
            tᵢ = params.t + params.Δt * params.t_q.c[i]
            simd_copy_xy_first!($tQ, Q, i)
            simd_copy_xy_first!($tV, V, i)
            params.α(tᵢ, $tQ, $tV, $tP)
            params.f(tᵢ, $tQ, $tV, $tF)
            simd_copy_yx_first!($tP, P, i)
            simd_copy_yx_first!($tF, F, i)
        end
    end

    return compute_stages_vprk
end


@generated function compute_rhs_vprk!(b::Vector{ST}, P::Matrix{ST}, F::Matrix{ST}, params::AbstractNonlinearFunctionParametersVPRK{DT,TT,AT,FT,D,S}) where {ST,DT,TT,AT,FT,D,S}
    compute_stages_vprk = quote
        local z1::ST
        local z2::ST

        # compute b = - [(P-AF)]
        for i in 1:S
            for k in 1:D
                z1 = 0
                z2 = 0
                for j in 1:S
                    z1 += params.t_p.a[i,j] * F[k,j]
                    z2 += params.t_p.â[i,j] * F[k,j]
                end
                b[D*(i-1)+k] = - (P[k,i] - params.p[k]) + params.Δt * (z1 + z2)
            end
        end
    end

    return compute_stages_vprk
end


@generated function compute_rhs_vprk!(b::Vector{ST}, P::Matrix{ST}, F::Matrix{ST}, G::Matrix{ST}, params::AbstractNonlinearFunctionParametersVPRK{DT,TT,AT,FT,D,S}) where {ST,DT,TT,AT,FT,D,S}
    compute_stages_vprk = quote
        local z1::ST
        local z2::ST

        # compute b = - [(P-G-AF)]
        for i in 1:S
            for k in 1:D
                z1 = 0
                z2 = 0
                for j in 1:S
                    z1 += params.t_p.a[i,j] * F[k,j]
                    z2 += params.t_p.â[i,j] * F[k,j]
                end
                b[D*(i-1)+k] = - (P[k,i] - params.p[k]) + params.Δt * (z1 + z2 + params.R[1] * G[k,1])
            end
        end
    end

    return compute_stages_vprk
end


@generated function compute_rhs_vprk!(b::Vector{ST}, P::Matrix{ST}, F::Matrix{ST}, R::Matrix{ST}, G::Matrix{ST}, params::AbstractNonlinearFunctionParametersVPRK{DT,TT,AT,FT,D,S}) where {ST,DT,TT,AT,FT,D,S}
    compute_stages_vprk = quote
        local z1::ST
        local z2::ST

        # compute b = - [(P-G-AF)]
        for i in 1:S
            for k in 1:D
                z1 = 0
                z2 = 0
                for j in 1:S
                    z1 += params.t_p.a[i,j] * (F[k,j] + R[k,j])
                    z2 += params.t_p.â[i,j] * (F[k,j] + R[k,j])
                end
                b[D*(i-1)+k] = - (P[k,i] - params.p[k]) + params.Δt * (z1 + z2 + params.R[1] * G[k,1])
            end
        end
    end

    return compute_stages_vprk
end


function compute_rhs_vprk_projection_q!(b::Vector{ST}, q̅::Vector{ST}, V::Matrix{ST}, U::Matrix{ST}, offset::Int, params::AbstractNonlinearFunctionParametersVPRK{DT,TT,AT,FT,D,S}) where {ST,DT,TT,AT,FT,D,S}
    local y1::ST
    local y2::ST

    for k in 1:D
        y1 = 0
        y2 = 0
        for j in 1:S
            y1 += params.t_q.b[j] * V[k,j]
            y2 += params.t_q.b̂[j] * V[k,j]
        end
        b[offset+k] = - (q̅[k] - params.q[k]) + params.Δt * (y1 + y2) + params.Δt * (params.R[1] * U[k,1] + params.R[2] * U[k,2])
    end
end


function compute_rhs_vprk_projection_p!(b::Vector{ST}, p̅::Vector{ST}, F::Matrix{ST}, G::Matrix{ST}, offset::Int, params::AbstractNonlinearFunctionParametersVPRK{DT,TT,AT,FT,D,S}) where {ST,DT,TT,AT,FT,D,S}
    local z1::ST
    local z2::ST

    for k in 1:D
        z1 = 0
        z2 = 0
        for j in 1:S
            z1 += params.t_p.b[j] * F[k,j]
            z2 += params.t_p.b̂[j] * F[k,j]
        end
        b[offset+k] = - (p̅[k] - params.p[k]) + params.Δt * (z1 + z2) + params.Δt * (params.R[1] * G[k,1] + params.R[2] * G[k,2])
    end
end


function compute_rhs_vprk_projection_p!(b::Vector{ST}, p̅::Vector{ST}, F::Matrix{ST}, R::Matrix{ST}, G::Matrix{ST}, offset::Int, params::AbstractNonlinearFunctionParametersVPRK{DT,TT,AT,FT,D,S}) where {ST,DT,TT,AT,FT,D,S}
    local z1::ST
    local z2::ST

    for k in 1:D
        z1 = 0
        z2 = 0
        for j in 1:S
            z1 += params.t_p.b[j] * (F[k,j] + R[k,j])
            z2 += params.t_p.b̂[j] * (F[k,j] + R[k,j])
        end
        b[offset+k] = - (p̅[k] - params.p[k]) + params.Δt * (z1 + z2) + params.Δt * (params.R[1] * G[k,1] + params.R[2] * G[k,2])
    end
end


@generated function compute_rhs_vprk_correction!(b::Vector{ST}, V::Matrix{ST}, params::AbstractNonlinearFunctionParametersVPRK{DT,TT,AT,FT,D,S}) where {ST,DT,TT,AT,FT,D,S}
    μ = zeros(ST,D)

    compute_stages_vprk = quote
        local sl::Int = div(S+1, 2)

        if length(params.d_v) > 0
            # compute μ
            for k in 1:D
                $μ[k] = params.t_p.b[sl] / params.d_v[sl] * b[D*(sl-1)+k]
            end

            # replace equation for Pₗ with constraint on V
            for k in 1:D
                b[D*(sl-1)+k] = 0
                for i in 1:S
                    b[D*(sl-1)+k] += V[k,i] * params.d_v[i]
                end
            end

            # modify P₁, ..., Pₛ except for Pₗ
            for i in 1:S
                if i ≠ sl
                    z = params.d_v[i] / params.t_p.b[i]
                    for k in 1:D
                        b[D*(i-1)+k] -= z * $μ[k]
                    end
                end
            end
        end
    end

    return compute_stages_vprk
end


function update_solution!(int::AbstractIntegratorVPRK{DT,TT}, cache::NonlinearFunctionCacheVPRK{DT}) where {DT,TT}
    update_solution!(int.q, int.qₑᵣᵣ, cache.V, int.tableau.q.b, int.tableau.q.b̂, int.Δt)
    update_solution!(int.p, int.pₑᵣᵣ, cache.F, int.tableau.p.b, int.tableau.p.b̂, int.Δt)
end


function project_solution!(int::AbstractIntegratorVPRK{DT,TT}, cache::NonlinearFunctionCacheVPRKprojection{DT}, R::Vector{TT}) where {DT,TT}
    update_solution!(int.q, int.qₑᵣᵣ, cache.U, R, int.Δt)
    update_solution!(int.p, int.pₑᵣᵣ, cache.G, R, int.Δt)
end


function cut_periodic_solution!(int::AbstractIntegratorVPRK)
    cut_periodic_solution!(int.q, int.qₑᵣᵣ, int.equation.periodicity)
end
