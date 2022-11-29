
abstract type SplittingMethod <: SODEMethod end

issodemethod(::SplittingMethod) = true

Integrators.Integrator(problem::SODEProblem, method::SplittingMethod; kwargs...) = Integrator(problem, tableau(method); kwargs...)


@doc raw"""
Lie-Trotter Splitting A

For a vector field $\dot{x} = f_1 (t,x) +  f_2 (t,x)$, the splitting reads
```math
\Phi_{\Delta t} = \varphi^{2}_{\Delta t} \circ \varphi^{1}_{\Delta t}
```

Reference:

    H. F. Trotter.
    On the product of semi-groups of operators.
    Proceedings of the American Mathematical Society, Volume 10, Pages 545-551, 1959.
    doi: 10.1090/S0002-9939-1959-0108732-6.

"""
struct LieA <: SplittingMethod end

function tableau(::LieA, ::Type{T}=Float64) where {T}
    a = Array{T}([ 1 ])
    b = Array{T}([ 0 ])
    TableauSplittingNS(:LieTrotterSplittingA, 1, a, b)
end


@doc raw"""
Lie-Trotter Splitting B

For a vector field $\dot{x} = f_1 (t,x) +  f_2 (t,x)$, the splitting reads
```math
\Phi_{\Delta t} = \varphi^{1}_{\Delta t} \circ \varphi^{2}_{\Delta t}
```

Reference:

    H. F. Trotter.
    On the product of semi-groups of operators.
    Proceedings of the American Mathematical Society, Volume 10, Pages 545-551, 1959.
    doi: 10.1090/S0002-9939-1959-0108732-6.

"""
struct LieB <: SplittingMethod end

function tableau(::LieB, ::Type{T}=Float64) where {T}
    a = Array{T}([ 0 ])
    b = Array{T}([ 1 ])
    TableauSplittingNS(:LieTrotterSplittingB, 1, a, b)
end


@doc raw"""
Strang Splitting

For a vector field $\dot{x} = f_1 (t,x) + f_2 (t,x)$, the splitting reads
```math
\Phi_{\Delta t} = \varphi^{1}_{\Delta t / 2} \circ \varphi^{2}_{\Delta t / 2} \circ \varphi^{2}_{\Delta t / 2} \circ \varphi^{1}_{\Delta t / 2}
```

For vector fields with two components, this is not the most efficient implementation.
For such cases [`StrangA`](@ref) or [`StrangB`](@ref) should be used instead.


References:

    Gilbert Strang.
    On the construction and comparison of difference schemes.
    SIAM Journal on Numerical Analysis, Volume 5, Pages 506-517, 1968.
    doi: 10.1137/0705041.

    Gurij Ivanovich Marchuk.
    Some applications of splitting-up methods to the solution of mathematical physics problems.
    Aplikace Matematiky, Volume 13, Pages 103-132, 1968.
    doi: 10.21136/AM.1968.103142.

"""
struct Strang <: SplittingMethod end

"Alias for [`Strang`](@ref)"
struct Marchuk <: SplittingMethod end

function tableau(::Union{Strang,Marchuk}, ::Type{T}=Float64) where {T}
    a = Array{T}([ 1//2 ])
    b = Array{T}([ 1//2 ])
    TableauSplittingNS(:StrangSplitting, 2, a, b)
end


@doc raw"""
Strang Splitting A for a vector field $\dot{x} = f_1 (t,x) + f_2 (t,x)$.

The splitting reads
```math
\Phi_{\Delta t} = \varphi^{1}_{\Delta t / 2} \circ \varphi^{2}_{\Delta t} \circ \varphi^{1}_{\Delta t / 2}
```

References:

    Gilbert Strang.
    On the construction and comparison of difference schemes.
    SIAM Journal on Numerical Analysis, Volume 5, Pages 506-517, 1968.
    doi: 10.1137/0705041.

    Gurij Ivanovich Marchuk.
    Some applications of splitting-up methods to the solution of mathematical physics problems.
    Aplikace Matematiky, Volume 13, Pages 103-132, 1968.
    doi: 10.21136/AM.1968.103142.

"""
struct StrangA <: SplittingMethod end

function tableau(::StrangA, ::Type{T}=Float64) where {T}
    a = Array{T}([ 1//2, 1//2 ])
    b = Array{T}([ 1//1, 0//1 ])
    TableauSplitting(:StrangSplittingA, 2, a, b)
end


@doc raw"""
Strang Splitting B for a vector field $\dot{x} = f_1 (t,x) + f_2 (t,x)$

The splitting reads
```math
\Phi_{\Delta t} = \varphi^{2}_{\Delta t / 2} \circ \varphi^{1}_{\Delta t} \circ \varphi^{2}_{\Delta t / 2}
```

References:

    Gilbert Strang.
    On the construction and comparison of difference schemes.
    SIAM Journal on Numerical Analysis, Volume 5, Pages 506-517, 1968.
    doi: 10.1137/0705041.

    Gurij Ivanovich Marchuk.
    Some applications of splitting-up methods to the solution of mathematical physics problems.
    Aplikace Matematiky, Volume 13, Pages 103-132, 1968.
    doi: 10.21136/AM.1968.103142.

"""
struct StrangB <: SplittingMethod end

function tableau(::StrangB, ::Type{T}=Float64) where {T}
    a = Array{T}([ 0//1, 1//1 ])
    b = Array{T}([ 1//2, 1//2 ])
    TableauSplitting(:StrangSplittingB, 2, a, b)
end


@doc raw"""
McLachlan's 2nd order symmetric, minimum error composition method

The composition reads
```math
\Phi_{\Delta t} = \varphi_{\alpha \Delta t} \circ \varphi^{*}_{(1/2 - \alpha) \Delta t} \circ \varphi_{(1/2 - \alpha) \Delta t} \circ \varphi^{*}_{\alpha \Delta t} ,
```
where the parameter $\alpha$ can be optimized, e.g., to minimize the solution error.
McLachlan arrives at $\alpha  = 0.1932$ as a generally useful value.

Reference:

    Robert I. McLachlan.
    On the Numerical Integration of Ordinary Differential Equations by Symmetric Composition Methods
    SIAM Journal on Scientific Computing, Volume 16, Pages 151-168, 1995.
    doi: 10.1137/0916010.

"""
struct McLachlan2 <: SplittingMethod end

function tableau(::McLachlan2, ::Type{T}=Float64; α=0.1932) where {T}
    a = Array{T}([ α, 0.5 - α ])
    b = Array{T}([ 0.5 - α, α ])
    TableauSplittingNS(:McLachlanSplitting, 2, a, b)
end


@doc raw"""
McLachlan's 4th order symmetric, minimum error composition method

The composition reads
```math
\Phi_{\Delta t} = \varphi_{\alpha_5 \Delta t} \circ \varphi^{*}_{\beta_5 \Delta t} \circ \dotsc \circ \varphi_{\alpha_2 \Delta t} \circ \varphi^{*}_{\beta_2 \Delta t} \circ \varphi_{\alpha_1 \Delta t} \circ \varphi^{*}_{\beta_1 \Delta t} ,
```
with
```math
\begin{aligned}
\beta_1 &= \alpha_5 = \frac{14 - \sqrt{19}}{108} , &
\alpha_1 &= \beta_5 = \frac{146 + 5 \sqrt{19}}{540} , & & \\
\beta_2 &= \alpha_4 = \frac{- 23 - 20 \sqrt{19}}{270} , &
\alpha_2 &= \beta_4 = \frac{-2 + 10 \sqrt{19}}{135} , & 
\beta_3 &= \alpha_3 = \frac{1}{5} .
\end{aligned}
```
The coefficients are optimised to provide an integrator with minimal solution error.

Reference:

    Robert I. McLachlan.
    On the Numerical Integration of Ordinary Differential Equations by Symmetric Composition Methods
    SIAM Journal on Scientific Computing, Volume 16, Pages 151-168, 1995.
    doi: 10.1137/0916010.

"""
struct McLachlan4 <: SplittingMethod end

function tableau(::McLachlan4, ::Type{T}=Float64) where {T}
    a = Array{T}(@big [ (146 +  5*√19) / 540,
                        ( -2 + 10*√19) / 135,
                                     1 / 5,
                        (-23 - 20*√19) / 270,
                        ( 14 -    √19) / 108])
    TableauSplittingNS(:McLachlanSplitting, 4, a, a[end:-1:1])
end

@doc raw"""
4th order "Triple Jump" composition method.

The composition reads
```math
\Phi_{\Delta t} = \varphi_{\gamma_3 \Delta t} \circ \varphi_{\gamma_2 \Delta t} \circ \varphi_{\gamma_1 \Delta t}
```
with
```math
\gamma_1 = \gamma_3 = \frac{1}{2 - 2^{1/(p+1)}} , \qquad
\gamma_2 = - \frac{2^{1/(p+1)}}{2 - 2^{1/(p+1)}} .
```

References:

    Michael Creutz and Andreas Gocksch.
    Higher-order hybrid Monte Carlo algorithms.
    Physical Review Letters, Volume 63, Pages 9-12, 1989.
    doi: 10.1103/PhysRevLett.63.9.

    Etienne Forest.
    Canonical integrators as tracking codes (or how to integrate perturbation theory with tracking).
    AIP Conference Proceedings, Volume 184, Pages 1106-1136, 1989.
    doi: 10.1063/1.38062.

    Masuo Suzuki
    Fractal decomposition of exponential operators with applications to many-body theories and Monte Carlo simulations.
    Physics Letters A, Volume 146, Pages 319-323, 1990.
    doi: 10.1016/0375-9601(90)90962-N

    Haruo Yoshida.
    Construction of higher order symplectic integrators.
    Physics Letters A, Volume 150, Pages 262-268, 1990.
    doi: 10.1016/0375-9601(90)90092-3

"""
struct TripleJump <: SplittingMethod end

function tableau(::TripleJump, ::Type{T}=Float64) where {T}
    fac = @big 2^(1/3)
    den = @big 1/(2-fac)
    a = Array{T}([ den, -fac*den ])
    TableauSplittingSS(:TripleJumpSplitting, 4, a)
end


@doc raw"""
Suzuki's 4th order "fractal" composition method

The composition reads
```math
\Phi_{\Delta t} = \varphi_{\gamma_5 \Delta t} \circ \varphi_{\gamma_4 \Delta t} \circ \varphi_{\gamma_3 \Delta t} \circ \varphi_{\gamma_2 \Delta t} \circ \varphi_{\gamma_1 \Delta t}
```
with
```math
\gamma_1 = \gamma_2 = \gamma_4 = \gamma_5 = \frac{1}{4 - 4^{1/(p+1)}} , \qquad
\gamma_3 = - \frac{4^{1/(p+1)}}{4 - 4^{1/(p+1)}} .
```

Reference:

    Masuo Suzuki
    Fractal decomposition of exponential operators with applications to many-body theories and Monte Carlo simulations.
    Physics Letters A, Volume 146, Pages 319-323, 1990.
    doi: 10.1016/0375-9601(90)90962-N

"""
struct SuzukiFractal <: SplittingMethod end

function tableau(::SuzukiFractal, ::Type{T}=Float64) where {T}
    fac = @big 4^(1/3)
    den = @big 1/(4-fac)
    a = Array{T}([ den, den, -fac*den ])
    TableauSplittingSS(:SuzukiFractalSplitting, 4, a)
end
