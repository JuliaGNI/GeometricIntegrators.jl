@doc raw"""
Midpoint extrapolation method with arbitrary order p.

For an [`ODEProblem`](@ref), this solves the ordinary differential equation

```math
\begin{aligned}
\dot{x} &= v(t, x) , &
x(t_0) &= x_0 ,
\end{aligned}
```

for $x_1 = x(t_1)$, and is called with

```julia
extrapolate!(t₀, x₀, t₁, x₁, ::ODEProblem, MidpointExtrapolation(s))
```

where

* `t₀`: initial time
* `x₀`: initial value $x_0 = x(t_0)$
* `t₁`: final   time
* `x₁`: final   value $x_1 = x(t_1)$
* `s`:  number of interpolations (order $p=2s+2$)


For a [`PODEProblem`](@ref) or [`HODEProblem`](@ref),
this solves the partitioned ordinary differential equation

```math
\begin{aligned}
\dot{q} &= v(t, q, p) , &
q(t_0) &= q_0 , \\
\dot{p} &= f(t, q, p) , &
p(t_0) &= p_0 , 
\end{aligned}
```

for $q_1 = q(t_1)$ and $p_1 = p(t_1)$m and is called with

```julia
extrapolate!(t₀, q₀, p₀, t₁, q₁, p₁, ::PODEProblem, MidpointExtrapolation(s))
extrapolate!(t₀, q₀, p₀, t₁, q₁, p₁, ::HODEProblem, MidpointExtrapolation(s))
```
    
where

* `t₀`: initial time
* `q₀`: initial position $q_0 = q(t_0)$
* `p₀`: initial momentum $p_0 = p(t_0)$
* `t₁`: final   time
* `q₁`: final   position $q_1 = q(t_1)$
* `p₁`: final   momentum $p_1 = p(t_1)$
* `s`:  number of interpolations (order $p=2s+2$)


Similarly, for a [`IODEProblem`](@ref) or [`LODEProblem`](@ref),
this solves the explicit dynamical equation

```math
\begin{aligned}
\dot{q} &= v(t, q) , &
q(t_0) &= q_0 , \\
\dot{p} &= f(t, q, v) , &
p(t_0) &= p_0 , 
\end{aligned}
```

corresponding to the implicit problem, for $q_1 = q(t_1)$ and $p_1 = p(t_1)$, and is called with

    ```julia
extrapolate!(t₀, q₀, p₀, t₁, q₁, p₁, ::IODEProblem, MidpointExtrapolation(s))
extrapolate!(t₀, q₀, p₀, t₁, q₁, p₁, ::LODEProblem, MidpointExtrapolation(s))
```

where

* `t₀`: initial time
* `q₀`: initial position $q_0 = q(t_0)$
* `p₀`: initial momentum $p_0 = p(t_0)$
* `t₁`: final   time
* `q₁`: final   position $q_1 = q(t_1)$
* `p₁`: final   momentum $p_1 = p(t_1)$
* `s`:  number of interpolations (order $p=2s+2$)

"""
struct MidpointExtrapolation <: Extrapolation
    s::Int
    MidpointExtrapolation(s=default_extrapolation_stages) = new(s)
end


function extrapolate_ode!(t₀::TT, x₀::AbstractVector{DT},
                      t₁::TT, x₁::AbstractVector{DT},
                      v::Callable, params::OptionalParameters,
                      extrap::MidpointExtrapolation) where {DT,TT}
    @assert size(x₀) == size(x₁)

    local F   = [2i*one(TT) for i in 1:extrap.s+1]
    local σ   = (t₁ - t₀) ./ F
    local σ²  = σ.^2
    local pts = zeros(DT, length(x₀), extrap.s+1)

    local xᵢ₁ = zero(x₀)
    local xᵢ₂ = zero(x₀)
    local xᵢₜ = zero(x₀)
    local vᵢ  = zero(x₀)
    local v₀  = zero(x₀)

    v(v₀, t₀, x₀, params)

    for i in 1:extrap.s+1
        tᵢ   = t₀ + σ[i]
        xᵢ₁ .= x₀
        xᵢ₂ .= x₀ .+ σ[i] .* v₀
        for _ in 1:(F[i]-1)
            v(vᵢ, tᵢ, xᵢ₂, params)
            xᵢₜ .= xᵢ₁ .+ 2σ[i] .* vᵢ
            xᵢ₁ .= xᵢ₂
            xᵢ₂ .= xᵢₜ
        end
        for k in axes(pts,1)
            pts[k,i] += xᵢ₂[k]
        end
    end

    aitken_neville!(x₁, zero(TT), σ², pts)
    
    return x₁
end

function extrapolate!(sol, history, problem::AbstractProblemODE, extrap::MidpointExtrapolation)
    extrapolate_ode!(history.t[1], history.q[1], sol.t, sol.q, functions(problem).v, parameters(problem), extrap)
end


function extrapolate_pode!(t₀::TT, q₀::AbstractVector{DT}, p₀::AbstractVector{DT}, 
                      t₁::TT, q₁::AbstractVector{DT}, p₁::AbstractVector{DT}, 
                      v::Callable, f::Callable, params::OptionalParameters,
                      extrap::MidpointExtrapolation) where {DT,TT}
    @assert size(q₀) == size(q₁) == size(p₀) == size(p₁)

    local F   = [2i*one(TT) for i in 1:extrap.s+1]
    local σ   = (t₁ - t₀) ./ F
    local σ2  = σ.^2

    local qts = zeros(DT, length(q₀), extrap.s+1)
    local pts = zeros(DT, length(p₀), extrap.s+1)

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

    v(v₀, t₀, q₀, p₀, params)
    f(f₀, t₀, q₀, p₀, params)

    for i in 1:extrap.s+1
        tᵢ   = t₀ + σ[i]
        qᵢ₁ .= q₀
        qᵢ₂ .= q₀ .+ σ[i] .* v₀
        pᵢ₁ .= p₀
        pᵢ₂ .= p₀ .+ σ[i] .* f₀
        for _ in 1:(F[i]-1)
            v(vᵢ, tᵢ, qᵢ₂, pᵢ₂, params)
            f(fᵢ, tᵢ, qᵢ₂, pᵢ₂, params)
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

    aitken_neville!(q₁, zero(TT), σ2, qts)
    aitken_neville!(p₁, zero(TT), σ2, pts)

    return (q₁, p₁)
end

function extrapolate!(sol, history, problem::AbstractProblemPODE, extrap::MidpointExtrapolation)
    extrapolate_pode!(history.t[1], history.q[1], history.p[1], sol.t, sol.q, sol.p, functions(problem).v, functions(problem).f, parameters(problem), extrap)
end


function extrapolate_iode!(t₀::TT, q₀::AbstractVector{DT}, p₀::AbstractVector{DT}, 
                      t₁::TT, q₁::AbstractVector{DT}, p₁::AbstractVector{DT}, 
                      v::Callable, f::Callable, params::OptionalParameters,
                      extrap::MidpointExtrapolation) where {DT,TT}
    @assert size(q₀) == size(q₁) == size(p₀) == size(p₁)

    local F   = [2i*one(TT) for i in 1:extrap.s+1]
    local σ   = (t₁ - t₀) ./ F
    local σ2  = σ.^2

    local qts = zeros(DT, length(q₀), extrap.s+1)
    local pts = zeros(DT, length(p₀), extrap.s+1)

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

    v(v₀, t₀, q₀, p₀, params)
    f(f₀, t₀, q₀, v₀, params)

    for i in 1:extrap.s+1
        tᵢ   = t₀ + σ[i]
        qᵢ₁ .= q₀
        qᵢ₂ .= q₀ .+ σ[i] .* v₀
        pᵢ₁ .= p₀
        pᵢ₂ .= p₀ .+ σ[i] .* f₀
        for _ in 1:(F[i]-1)
            v(vᵢ, tᵢ, qᵢ₂, pᵢ₂, params)
            f(fᵢ, tᵢ, qᵢ₂, vᵢ, params)
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

    aitken_neville!(q₁, zero(TT), σ2, qts)
    aitken_neville!(p₁, zero(TT), σ2, pts)

    return (q₁, p₁)
end

function extrapolate!(sol, history, problem::AbstractProblemIODE, extrap::MidpointExtrapolation)
    extrapolate_iode!(history.t[1], history.q[1], history.p[1], sol.t, sol.q, sol.p, functions(problem).v̄, functions(problem).f̄, parameters(problem), extrap)
end
