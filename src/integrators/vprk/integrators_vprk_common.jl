
abstract type AbstractParametersVPRK{DT,TT,ET,D,S} <: Parameters{DT,TT} end
abstract type AbstractIntegratorVPRK{DT,TT} <: Integrator{DT,TT} end

equation(integrator::AbstractIntegratorVPRK) = integrator.params.equ
timestep(integrator::AbstractIntegratorVPRK) = integrator.params.Δt
tableau(integrator::AbstractIntegratorVPRK) = integrator.params.tab
dims(integrator::AbstractIntegratorVPRK) = integrator.params.equ.d


struct NonlinearFunctionCacheVPRK{ST}
    Q::Matrix{ST}
    V::Matrix{ST}
    P::Matrix{ST}
    F::Matrix{ST}

    v::Vector{ST}
    y::Vector{ST}
    z::Vector{ST}

    q̅::Vector{ST}
    p̅::Vector{ST}

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

        # create solution vector
        q̅ = zeros(ST,D)
        p̅ = zeros(ST,D)

        new(Q, V, P, F, v, y, z, q̅, p̅)
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


function compute_stages!(x, Q, V, P, F, params::AbstractParametersVPRK)
    # copy x to V
    compute_stages_v_vprk!(x, V, params)

    # compute Q
    compute_stages_q_vprk!(Q, V, params)

    # compute P and F
    compute_stages_p_vprk!(Q, V, P, F, params)
end

function compute_stages!(x, q̅, p̅, λ, Q, V, U, P, F, G, params::AbstractParametersVPRK)
    # copy x to V
    compute_stages_v_vprk!(x, V, params)

    # compute U, G and p̅
    compute_projection_vprk!(x, q̅, p̅, λ, V, U, G, params)

    # compute Q
    compute_stages_q_vprk!(Q, V, U, params)

    # compute P and F
    compute_stages_p_vprk!(Q, V, P, F, params)
end

function compute_stages_v_vprk!(x::Vector{ST}, V::Matrix{ST}, params::AbstractParametersVPRK{DT,TT,ET,D,S}) where {ST,DT,TT,ET,D,S}
    @assert D == size(V,1)
    @assert S == size(V,2)

    # copy x to V
    for i in 1:S
        for k in 1:D
            V[k,i] = x[D*(i-1)+k]
        end
    end
end

function compute_stages_λ_vprk!(x::Vector{ST}, Λ::Matrix{ST}, params::AbstractParametersVPRK{DT,TT,ET,D,S}) where {ST,DT,TT,ET,D,S}
    @assert D == size(Λ,1)
    @assert S == size(Λ,2)

    # copy x to Λ
    for i in 1:S
        for k in 1:D
            Λ[k,i] = x[D*(S+i-1)+k]
        end
    end
end

function compute_stages_q_vprk!(Q::Matrix{ST}, V::Matrix{ST}, params::AbstractParametersVPRK{DT,TT,ET,D,S}) where {ST,DT,TT,ET,D,S}
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
                y1 += params.tab.q.a[i,j] * V[k,j]
                y2 += params.tab.q.â[i,j] * V[k,j]
            end
            Q[k,i] = params.q[k] + params.Δt * (y1 + y2)
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


@generated function compute_stages_p_vprk!(Q::Matrix{ST}, V::Matrix{ST}, P::Matrix{ST}, F::Matrix{ST}, params::AbstractParametersVPRK{DT,TT,ET,D,S}) where {ST,DT,TT,ET,D,S}
    # create temporary vectors
    tQ = zeros(ST,D)
    tV = zeros(ST,D)
    tP = zeros(ST,D)
    tF = zeros(ST,D)

    compute_stages_vprk = quote
        @assert D == size(Q,1) == size(V,1) == size(P,1) == size(F,1)
        @assert S == size(Q,2) == size(V,2) == size(P,2) == size(F,2)

        local tᵢ::TT

        # compute P=α(Q) and F=f(Q,V)
        for i in 1:S
            tᵢ = params.t + params.Δt * params.tab.q.c[i]
            simd_copy_xy_first!($tQ, Q, i)
            simd_copy_xy_first!($tV, V, i)
            params.equ.α(tᵢ, $tQ, $tV, $tP)
            params.equ.f(tᵢ, $tQ, $tV, $tF)
            simd_copy_yx_first!($tP, P, i)
            simd_copy_yx_first!($tF, F, i)
        end
    end

    return compute_stages_vprk
end


@generated function compute_rhs_vprk!(b::Vector{ST}, P::Matrix{ST}, F::Matrix{ST}, params::AbstractParametersVPRK{DT,TT,ET,D,S}) where {ST,DT,TT,ET,D,S}
    compute_stages_vprk = quote
        local z1::ST
        local z2::ST

        # compute b = - [(P-AF)]
        for i in 1:S
            for k in 1:D
                z1 = 0
                z2 = 0
                for j in 1:S
                    z1 += params.tab.p.a[i,j] * F[k,j]
                    z2 += params.tab.p.â[i,j] * F[k,j]
                end
                b[D*(i-1)+k] = - (P[k,i] - params.p[k]) + params.Δt * (z1 + z2)
                # b[D*(i-1)+k]  = - (P[k,i] - params.p[k])
                # b[D*(i-1)+k] += params.Δt * z1
                # b[D*(i-1)+k] += params.Δt * z2
            end
        end
    end

    return compute_stages_vprk
end


@generated function compute_rhs_vprk!(b::Vector{ST}, P::Matrix{ST}, F::Matrix{ST}, G::Matrix{ST}, params::AbstractParametersVPRK{DT,TT,ET,D,S}) where {ST,DT,TT,ET,D,S}
    compute_stages_vprk = quote
        local z1::ST
        local z2::ST

        # compute b = - [(P-G-AF)]
        for i in 1:S
            for k in 1:D
                z1 = 0
                z2 = 0
                for j in 1:S
                    z1 += params.tab.p.a[i,j] * F[k,j]
                    z2 += params.tab.p.â[i,j] * F[k,j]
                end
                b[D*(i-1)+k] = - (P[k,i] - params.p[k]) + params.Δt * (z1 + z2 + params.R[1] * G[k,1])
                # b[D*(i-1)+k]  = - (P[k,i] - params.p[k])
                # b[D*(i-1)+k] += params.Δt * z1
                # b[D*(i-1)+k] += params.Δt * params.R[1] * G[k,1]
                # b[D*(i-1)+k] += params.Δt * z2
            end
        end
    end

    return compute_stages_vprk
end


function compute_rhs_vprk!(b::Vector{ST}, P::Matrix{ST}, F::Matrix{ST}, R::Matrix{ST}, G::Matrix{ST}, params::AbstractParametersVPRK{DT,TT,ET,D,S}) where {ST,DT,TT,ET,D,S}
    local z1::ST
    local z2::ST

    # compute b = - [(P-G-AF)]
    for i in 1:S
        for k in 1:D
            z1 = 0
            z2 = 0
            for j in 1:S
                z1 += params.tab.p.a[i,j] * (F[k,j] + R[k,j])
                z2 += params.tab.p.â[i,j] * (F[k,j] + R[k,j])
            end
            b[D*(i-1)+k] = - (P[k,i] - params.p[k]) + params.Δt * (z1 + z2 + params.R[1] * G[k,1])
            # b[D*(i-1)+k]  = - (P[k,i] - params.p[k])
            # b[D*(i-1)+k] += params.Δt * z1
            # b[D*(i-1)+k] += params.Δt * params.R[1] * G[k,1]
            # b[D*(i-1)+k] += params.Δt * z2
        end
    end
end


function compute_rhs_vprk_projection_q!(b::Vector{ST}, q̅::Vector{ST}, V::Matrix{ST}, U::Matrix{ST}, offset::Int, params::AbstractParametersVPRK{DT,TT,ET,D,S}) where {ST,DT,TT,ET,D,S}
    local y1::ST
    local y2::ST

    for k in 1:D
        y1 = 0
        y2 = 0
        for j in 1:S
            y1 += params.tab.q.b[j] * V[k,j]
            y2 += params.tab.q.b̂[j] * V[k,j]
        end
        b[offset+k] = - (q̅[k] - params.q[k]) + params.Δt * (y1 + y2) + params.Δt * (params.R[1] * U[k,1] + params.R[2] * U[k,2])
        # b[offset+k]  = - (q̅[k] - params.q[k])
        # b[offset+k] += params.Δt * y1
        # b[offset+k] += params.Δt * (params.R[1] * U[k,1] + params.R[2] * U[k,2])
        # b[offset+k] += params.Δt * y2
    end
end


function compute_rhs_vprk_projection_p!(b::Vector{ST}, p̅::Vector{ST}, F::Matrix{ST}, G::Matrix{ST}, offset::Int, params::AbstractParametersVPRK{DT,TT,ET,D,S}) where {ST,DT,TT,ET,D,S}
    local z1::ST
    local z2::ST

    for k in 1:D
        z1 = 0
        z2 = 0
        for j in 1:S
            z1 += params.tab.p.b[j] * F[k,j]
            z2 += params.tab.p.b̂[j] * F[k,j]
        end
        b[offset+k] = - (p̅[k] - params.p[k]) + params.Δt * (z1 + z2) + params.Δt * (params.R[1] * G[k,1] + params.R[2] * G[k,2])
        # b[offset+k]  = - (p̅[k] - params.p[k])
        # b[offset+k] += params.Δt * z1
        # b[offset+k] += params.Δt * (params.R[1] * G[k,1] + params.R[2] * G[k,2])
        # b[offset+k] += params.Δt * z2
    end
end


function compute_rhs_vprk_projection_p!(b::Vector{ST}, p̅::Vector{ST}, F::Matrix{ST}, R::Matrix{ST}, G::Matrix{ST}, offset::Int, params::AbstractParametersVPRK{DT,TT,ET,D,S}) where {ST,DT,TT,ET,D,S}
    local z1::ST
    local z2::ST

    for k in 1:D
        z1 = 0
        z2 = 0
        for j in 1:S
            z1 += params.tab.p.b[j] * (F[k,j] + R[k,j])
            z2 += params.tab.p.b̂[j] * (F[k,j] + R[k,j])
        end
        b[offset+k] = - (p̅[k] - params.p[k]) + params.Δt * (z1 + z2) + params.Δt * (params.R[1] * G[k,1] + params.R[2] * G[k,2])
        # b[offset+k]  = - (p̅[k] - params.p[k])
        # b[offset+k] += params.Δt * z1
        # b[offset+k] += params.Δt * (params.R[1] * G[k,1] + params.R[2] * G[k,2])
        # b[offset+k] += params.Δt * z2
    end
end


@generated function compute_rhs_vprk_correction!(b::Vector{ST}, V::Matrix{ST}, params::AbstractParametersVPRK{DT,TT,ET,D,S}) where {ST,DT,TT,ET,D,S}
    μ = zeros(ST,D)

    quote
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
                    b[D*(sl-1)+k] += V[k,i] * params.tab.d[i]
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
