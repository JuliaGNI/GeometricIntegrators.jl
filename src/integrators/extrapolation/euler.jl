
"""
Euler extrapolation method with arbitrary order p.

    v:  function to compute vector field
    t₀: initial time
    t₁: final   time
    x₀: initial value
    x₁: final   value
    s:  number of interpolations (order p=s+1)

"""
function euler_extrapolation(v::Function, t₀::TT, t₁::TT, x₀::Vector, x₁::Vector, s::Int) where {TT}
    @assert size(x₀) == size(x₁)

    F   = collect(1:(s+1))
    Δt  = t₁ - t₀
    σ   = Δt ./ F
    pts = repeat(x₀, outer = [1, s+1])

    local xᵢ = zero(x₀)
    local vᵢ = zero(x₀)

    for i in F
        for _ in 1:(F[i]-1)
            tᵢ = t₀ + σ[i]
            for k in axes(pts,1)
                xᵢ[k] = pts[k,i]
            end
            v(tᵢ, xᵢ, vᵢ)
            for k in axes(pts,1)
                pts[k,i] += σ[i] * vᵢ[k]
            end
        end
    end

    aitken_neville!(σ, pts, zero(TT), x₁)
end
