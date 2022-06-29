
@doc raw"""
Midpoint extrapolation method with arbitrary order p.

Solves the ordinary differential equation
```math
\begin{aligned}
\dot{x} &= v(t, x) , &
x(t_0) &= x_0 ,
\end{aligned}
```
for $x_1 = x(t_1)$.

Call with
```julia
_midpoint_extrapolation_ode!(v, t₀, t₁, x₀, x₁, s)
```

where

* `v`:  function to compute vector field with signature `v(t,x,ẋ)`
* `t₀`: initial time
* `t₁`: final   time
* `x₀`: initial value $x_0 = x(t_0)$
* `x₁`: final   value $x_1 = x(t_1)$
* `s`:  number of interpolations (order $p=2s+2$)

"""
function _midpoint_extrapolation_ode!(v::Function, t₀::TT, t₁::TT,
                                      x₀::AbstractVector{DT}, x₁::AbstractVector{DT}, s::Int) where {DT,TT}
    @assert size(x₀) == size(x₁)

    local F   = [2i*one(TT) for i in 1:(s+1)]
    local Δt  = t₁ - t₀
    local σ   = Δt ./ F
    local σ²  = σ.^2
    local pts = zeros(DT, length(x₀), s+1)

    local xᵢ₁ = zero(x₀)
    local xᵢ₂ = zero(x₀)
    local xᵢₜ = zero(x₀)
    local vᵢ  = zero(x₀)
    local v₀  = zero(x₀)

    v(t₀, x₀, v₀)

    for i in 1:s+1
        tᵢ   = t₀ + σ[i]
        xᵢ₁ .= x₀
        xᵢ₂ .= x₀ .+ σ[i] .* v₀
        for _ in 1:(F[i]-1)
            v(tᵢ, xᵢ₂, vᵢ)
            xᵢₜ .= xᵢ₁ .+ 2σ[i] .* vᵢ
            xᵢ₁ .= xᵢ₂
            xᵢ₂ .= xᵢₜ
        end
        for k in axes(pts,1)
            pts[k,i] += xᵢ₂[k]
        end
    end

    aitken_neville!(σ², pts, zero(TT), x₁)
    return x₁
end


struct MidpointExtrapolationODE{VT} <: Extrapolation
    s::Int
    v::VT

    function MidpointExtrapolationODE{VT}(v, s) where {VT}
        new(s,v)
    end
end

function MidpointExtrapolationODE(v::VT, s::Int) where {VT}
    MidpointExtrapolationODE{VT}(v, s)
end

function MidpointExtrapolation(equ::ODEProblem, s::Int)
    MidpointExtrapolationODE(functions(equ).v, s)
end

function GeometricBase.evaluate!(extrap::MidpointExtrapolationODE, t₀::TT, t₁::TT,
                          x₀::AbstractVector{DT}, x₁::AbstractVector{DT}) where {DT,TT}
    _midpoint_extrapolation_ode!(extrap.v, t₀, t₁, x₀, x₁, extrap.s)
end


@doc raw"""
Midpoint extrapolation method with arbitrary order p.

Solves the implicit differential equation
```math
\begin{aligned}
\dot{q} &= v(t, q, p) , &
q(t_0) &= q_0 , \\
\dot{p} &= f(t, q, v) , &
p(t_0) &= p_0 , 
\end{aligned}
```
for $q_1 = q(t_1)$ and $p_1 = p(t_1)$.

Call with
```julia
_midpoint_extrapolation_iode!(v, t₀, t₁, q₀, q₁, p₀, p₁, s)
```

where

* `v`:  function to compute vector field with signature `v(t,x,p,ẋ)`
* `f`:  function to compute force  field with signature `f(t,x,v,ṗ)`
* `t₀`: initial time
* `t₁`: final   time
* `q₀`: initial position $q_0 = q(t_0)$
* `q₁`: final   position $q_1 = q(t_1)$
* `p₀`: initial momentum $p_0 = p(t_0)$
* `p₁`: final   momentum $p_1 = p(t_1)$
* `s`:  number of interpolations (order $p=2s+2$)

"""
function _midpoint_extrapolation_iode!(v::Function, f::Function, t₀::TT, t₁::TT,
                                       q₀::AbstractVector{DT}, q₁::AbstractVector{DT},
                                       p₀::AbstractVector{DT}, p₁::AbstractVector{DT}, s::Int) where {DT,TT}
    @assert size(q₀) == size(q₁) == size(p₀) == size(p₁)

    local F   = [2i*one(TT) for i in 1:(s+1)]
    local Δt  = t₁ - t₀
    local σ   = Δt ./ F
    local σ2  = σ.^2

    local qts = zeros(DT, length(q₀), s+1)
    local pts = zeros(DT, length(p₀), s+1)

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

    v(t₀, q₀, v₀)
    f(t₀, q₀, v₀, f₀)

    for i in 1:(s+1)
        tᵢ   = t₀ + σ[i]
        qᵢ₁ .= q₀
        qᵢ₂ .= q₀ .+ σ[i] .* v₀
        pᵢ₁ .= p₀
        pᵢ₂ .= p₀ .+ σ[i] .* f₀
        for _ in 1:(F[i]-1)
            v(tᵢ, qᵢ₂, vᵢ)
            f(tᵢ, qᵢ₂, vᵢ,  fᵢ)
            qᵢₜ .= qᵢ₁ .+ 2σ[i] .* vᵢ
            qᵢ₁ .= qᵢ₂
            qᵢ₂ .= qᵢₜ
            pᵢₜ .= pᵢ₁ .+ 2σ[i] .* fᵢ
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

    aitken_neville!(σ2, qts, zero(TT), q₁)
    aitken_neville!(σ2, pts, zero(TT), p₁)
    return (q₁, p₁)
end


struct MidpointExtrapolationIODE{VT,FT} <: Extrapolation
    s::Int
    v::VT
    f::FT
    
    function MidpointExtrapolationIODE{VT,FT}(v, f, s) where {VT,FT}
        new(s,v,f)
    end
end

function MidpointExtrapolationIODE(v::VT, f::FT, s::Int) where {VT,FT}
    MidpointExtrapolationIODE{VT,FT}(v, f, s)
end

function MidpointExtrapolation(equ::IODEProblem, s::Int)
    MidpointExtrapolationIODE(functions(equ).v̄, functions(equ).f̄, s)
end

function GeometricBase.evaluate!(extrap::MidpointExtrapolationIODE, t₀::TT, t₁::TT,
                          q₀::AbstractVector{DT}, q₁::AbstractVector{DT},
                          p₀::AbstractVector{DT}, p₁::AbstractVector{DT}) where {DT,TT}
    _midpoint_extrapolation_iode!(extrap.v, extrap.f, t₀, t₁, q₀, q₁, p₀, p₁, extrap.s)
end


@doc raw"""
Midpoint extrapolation method with arbitrary order p.

Solves the partitioned ordinary differential equation
```math
\begin{aligned}
\dot{q} &= v(t, q, p) , &
q(t_0) &= q_0 , \\
\dot{p} &= f(t, q, p) , &
p(t_0) &= p_0 , 
\end{aligned}
```
for $q_1 = q(t_1)$ and $p_1 = p(t_1)$.

Call with
```julia
_midpoint_extrapolation_pode!(v, t₀, t₁, q₀, q₁, p₀, p₁, s)
```

where

* `v`:  function to compute velocity field with signature `v(t,x,p,ẋ)`
* `f`:  function to compute force    field with signature `f(t,x,p,ṗ)`
* `t₀`: initial time
* `t₁`: final   time
* `q₀`: initial position $q_0 = q(t_0)$
* `q₁`: final   position $q_1 = q(t_1)$
* `p₀`: initial momentum $p_0 = p(t_0)$
* `p₁`: final   momentum $p_1 = p(t_1)$
* `s`:  number of interpolations (order $p=2s+2$)

"""
function _midpoint_extrapolation_pode!(v::Function, f::Function, t₀::TT, t₁::TT,
                                       q₀::AbstractVector{DT}, q₁::AbstractVector{DT},
                                       p₀::AbstractVector{DT}, p₁::AbstractVector{DT}, s::Int) where {DT,TT}
    @assert size(q₀) == size(q₁) == size(p₀) == size(p₁)

    local F   = [2i*one(TT) for i in 1:(s+1)]
    local Δt  = t₁ - t₀
    local σ   = Δt ./ F
    local σ2  = σ.^2

    local qts = zeros(DT, length(q₀), s+1)
    local pts = zeros(DT, length(p₀), s+1)

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
    f(t₀, q₀, p₀, f₀)

    for i in 1:(s+1)
        tᵢ   = t₀ + σ[i]
        qᵢ₁ .= q₀
        qᵢ₂ .= q₀ .+ σ[i] .* v₀
        pᵢ₁ .= p₀
        pᵢ₂ .= p₀ .+ σ[i] .* f₀
        for _ in 1:(F[i]-1)
            v(tᵢ, qᵢ₂, pᵢ₂, vᵢ)
            f(tᵢ, qᵢ₂, pᵢ₂, fᵢ)
            qᵢₜ .= qᵢ₁ .+ 2σ[i] .* vᵢ
            qᵢ₁ .= qᵢ₂
            qᵢ₂ .= qᵢₜ
            pᵢₜ .= pᵢ₁ .+ 2σ[i] .* fᵢ
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

    aitken_neville!(σ2, qts, zero(TT), q₁)
    aitken_neville!(σ2, pts, zero(TT), p₁)
    return (q₁, p₁)
end


struct MidpointExtrapolationPODE{VT,FT} <: Extrapolation
    s::Int
    v::VT
    f::FT
    
    function MidpointExtrapolationPODE{VT,FT}(v, f, s) where {VT,FT}
        new(s,v,f)
    end
end

function MidpointExtrapolationPODE(v::VT, f::FT, s::Int) where {VT,FT}
    MidpointExtrapolationPODE{VT,FT}(v, f, s)
end

function MidpointExtrapolation(equ::PODEProblem, s::Int)
    MidpointExtrapolationPODE(functions(equ).v, functions(equ).f, s)
end

function GeometricBase.evaluate!(extrap::MidpointExtrapolationPODE, t₀::TT, t₁::TT,
                          q₀::AbstractVector{DT}, q₁::AbstractVector{DT},
                          p₀::AbstractVector{DT}, p₁::AbstractVector{DT}) where {DT,TT}
    _midpoint_extrapolation_pode!(extrap.v, extrap.f, t₀, t₁, q₀, q₁, p₀, p₁, extrap.s)
end
