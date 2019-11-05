
"""
Compute p(x) where p is the unique polynomial of degree length(xi),
such that p(x[i]) = y[i]) for all i.

    ti: interpolation nodes
    xi: interpolation values
    t:  evaluation point
    x:  evaluation value
"""
function aitken_neville(ti::Vector{TT}, xi::Matrix, t::TT, x::Vector) where {TT}
    @assert length(ti) == size(xi,2)
    @assert length(x)  == size(xi,1)

    for j in axes(t)
        for i in 1:(length(t)-j)
            for k in axes(x,1)
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
function euler_extrapolation(v::Function, t₀::TT, t₁::TT, x₀::Vector, x₁::Vector, s::Int) where {TT}
    @assert size(x₀) == size(x₁)

    F   = collect(1:(s+1))
    Δt  = t₁ - t₀
    σ   = Δt ./ F
    pts = repmat(x₀, 1, s+1)

    local xᵢ = zero(x₀)
    local vᵢ = zero(x₀)

    for i in 1:(s+1)
        for j in 1:(F[i]-1)
            tᵢ  = t₀ + σ[i]
            for k in axes(pts,1)
                xᵢ[k] = pts[k,i]
            end
            v(tᵢ, xᵢ, vᵢ)
            for k in axes(pts,1)
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
function midpoint_extrapolation(v::Function, t₀::TT, t₁::TT, x₀::Vector, x₁::Vector, s::Int) where {TT}
    @assert size(x₀) == size(x₁)

    local F   = 2*collect(1:(s+1))
    local Δt  = t₁ - t₀
    local σ   = Δt ./ F
    local σ²  = σ.^2
    local pts = zeros(eltype(x₀), length(x₀), s+1)

    local xᵢ₁= zero(x₀)
    local xᵢ₂= zero(x₀)
    local xᵢₜ= zero(x₀)
    local vᵢ = zero(x₀)
    local v₀ = zero(x₀)

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
        for k in axes(pts,1)
            pts[k,i] += xᵢ₂[k]
        end
    end

    aitken_neville(σ², pts, zero(TT), x₁)
end


"""
Midpoint extrapolation method with arbitrary order p.

    v:  function to compute vector field
    f:  function to compute force  field
    t₀: initial time
    t₁: final   time
    q₀: initial positions
    p₀: initial momenta
    q₁: final   positions
    p₁: final   momenta
    s:  number of interpolations (order p=2s+2)
"""
function midpoint_extrapolation(v::Function, f::Function, t₀::TT, t₁::TT, q₀::Vector, q₁::Vector, p₀::Vector, p₁::Vector, s::Int) where {TT}
    @assert size(q₀) == size(q₁) == size(p₀) == size(p₁)

    local F   = 2*collect(1:(s+1))
    local Δt  = t₁ - t₀
    local σ   = Δt ./ F
    local σ2  = σ.^2

    local qts = zeros(eltype(q₀), length(q₀), s+1)
    local pts = zeros(eltype(p₀), length(p₀), s+1)

    local qᵢ₁= zero(q₀)
    local qᵢ₂= zero(q₀)
    local qᵢₜ= zero(q₀)

    local pᵢ₁= zero(p₀)
    local pᵢ₂= zero(p₀)
    local pᵢₜ= zero(p₀)

    local v₀ = zero(q₀)
    local vᵢ = zero(q₀)

    local f₀ = zero(p₀)
    local fᵢ = zero(p₀)

    v(t₀, q₀, p₀, v₀)
    f(t₀, q₀, v₀, f₀)

    for i in 1:(s+1)
        tᵢ   = t₀ + σ[i]
        qᵢ₁ .= q₀
        qᵢ₂ .= q₀ + σ[i] .* v₀
        pᵢ₁ .= p₀
        pᵢ₂ .= p₀ + σ[i] .* f₀
        for j in 1:(F[i]-1)
            v(tᵢ, qᵢ₂, pᵢ₂, vᵢ)
            f(tᵢ, qᵢ₂, vᵢ,  fᵢ)
            qᵢₜ .= qᵢ₁ + 2σ[i] .* vᵢ
            qᵢ₁ .= qᵢ₂
            qᵢ₂ .= qᵢₜ
            pᵢₜ .= pᵢ₁ + 2σ[i] .* fᵢ
            pᵢ₁ .= pᵢ₂
            pᵢ₂ .= pᵢₜ
        end
        for k in axes(qts,1)
            qts[k,i] += qᵢ₂[k]
        end
        for k in axes(pts,1)
            pts[k,i] += pᵢ₂[k]
        end
    end

    aitken_neville(σ2, qts, zero(TT), q₁)
    aitken_neville(σ2, pts, zero(TT), p₁)
end
