@doc raw"""
Euler extrapolation method with arbitrary order p.

Solves the ordinary differential equation
```math
\begin{aligned}
\dot{x} &= v(t, x) , &
x(t_0) &= x_0 ,
\end{aligned}
```

for $x_1 = x(t_1)$, and is called with

```julia
extrapolate!(t₀, x₀, t₁, x₁, problem, EulerExtrapolation(s))
```

where

* `t₀`: initial time
* `t₁`: final   time
* `x₀`: initial value $x_0 = x(t_0)$
* `x₁`: final   value $x_1 = x(t_1)$
* `problem`: [`ODEProblem`](@ref) whose solution to extrapolate
* `s`:  number of interpolations (order $p=s+1$)

"""
struct EulerExtrapolation <: Extrapolation
    s::Int
    function EulerExtrapolation(s=default_extrapolation_stages)
        @assert s ≥ 0
        new(s)
    end
end


function extrapolate!(t₀::TT, x₀::AbstractArray{DT},
                      t₁::TT, x₁::AbstractArray{DT},
                      problem::AbstractProblemODE,
                      extrap::EulerExtrapolation) where {DT,TT}

    @assert axes(x₀) == axes(x₁)

    local F   = collect(1:(extrap.s+1))
    local σ   = (t₁ - t₀) ./ F
    local pts = [copy(x₀) for _ in F]

    local vᵢ = zero(x₀)

    for i in F
        for _ in 1:(F[i]-1)
            tᵢ  = t₀ + σ[i]
            initialguess(problem).v(vᵢ, tᵢ, pts[i], parameters(problem))
            pts[i] .+= σ[i] * vᵢ
        end
    end

    aitken_neville!(x₁, zero(TT), σ, pts)

    return x₁
end

function solutionstep!(sol, history, problem::Union{AbstractProblemODE, SODEProblem}, extrap::EulerExtrapolation)
    extrapolate!(history.t[1], history.q[1], sol.t, sol.q, problem, extrap)
    update_vectorfields!(sol, problem)
    return sol
end
