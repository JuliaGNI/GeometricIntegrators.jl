
"""
Compute p(x) where p is the unique polynomial of degree length(xi),
such that p(x[i]) = y[i]) for all i.

    ti: interpolation nodes
    xi: interpolation values
    t:  evaluation point
    x:  evaluation value
"""
function aitken_neville(ti::Vector{TT}, xi::Matrix{DT}, t::TT, x::Vector{DT}) where {DT,TT}
    @assert length(ti) == size(xi,2)
    @assert length(x)  == size(xi,1)

    for j in indices(t)
        for i in 1:(length(t)-j)
            for k in indices(x,1)
                xi[k,i] = xi[k,i+1] + (xi[k,i] - xi[k,i+1]) * (ti[i+j] - t) / (ti[i+j] - ti[i])
            end
        end
    end
    for k in eachindex(x)
        x[k] = xi[k,1]
    end
end


"""
Euler extrapolation method with arbitrary order p.

    v:  function to compute vector field
    t₀: initial time
    t₁: final   time
    x₀: initial value
    x₁: final   value
    s:  number of interpolations (order p=s+1)

# TODO This is probably broken!
"""
function euler_extrapolation(v::Function, t₀::TT, t₁::TT, x₀::Vector{DT}, x₁::Vector{DT}, s::Int) where {DT,TT}
    @assert size(x₀) == size(x₁)

    F   = collect(1:(s+1))
    Δt  = t₁ - t₀
    σ   = Δt ./ F
    pts = repmat(x₀, 1, s+1)

    local xᵢ = zeros(x₀)
    local vᵢ = zeros(x₀)

    for i in 1:(s+1)
        for j in 1:(F[i]-1)
            tᵢ  = t₀ + σ[i]
            for k in indices(pts,1)
                xᵢ[k] = pts[k,i]
            end
            v(tᵢ, xᵢ, vᵢ)
            for k in indices(pts,1)
                pts[k,i] += σ[i] * vᵢ[k]
            end
        end
    end

    aitken_neville(σ, pts, zero(TT), x₁)
end


"""
Midpoint extrapolation method with arbitrary order p.

    v:  function to compute vector field
    t₀: initial time
    t₁: final   time
    x₀: initial value
    x₁: final   value
    s:  number of interpolations (order p=2s+2)
"""
function midpoint_extrapolation(v::Function, t₀::TT, t₁::TT, x₀::Vector{DT}, x₁::Vector{DT}, s::Int) where {DT,TT}
    @assert size(x₀) == size(x₁)

    local F   = 2*collect(1:(s+1))
    local Δt  = t₁ - t₀
    local σ   = Δt ./ F
    local σ²  = σ.^2
    local pts = zeros(eltype(x₀), length(x₀), s+1)

    local xᵢ₁= zeros(x₀)
    local xᵢ₂= zeros(x₀)
    local xᵢₜ= zeros(x₀)
    local vᵢ = zeros(x₀)
    local v₀ = zeros(x₀)

    v(t₀, x₀, v₀)

    for i in 1:s+1
        tᵢ   = t₀ + σ[i]
        xᵢ₁ .= x₀
        xᵢ₂ .= x₀ + σ[i] .* v₀
        for j in 1:(F[i]-1)
            v(tᵢ, xᵢ₂, vᵢ)
            xᵢₜ .= xᵢ₁ + 2σ[i] .* vᵢ
            xᵢ₁ .= xᵢ₂
            xᵢ₂ .= xᵢₜ
        end
        for k in indices(pts,1)
            pts[k,i] += xᵢ₂[k]
        end
    end

    aitken_neville(σ², pts, zero(TT), x₁)
end
