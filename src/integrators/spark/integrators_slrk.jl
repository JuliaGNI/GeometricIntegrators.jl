
"Holds the tableau of an Specialised Partitioned Additive Runge-Kutta method for Variational systems."
struct TableauSLRK{DT <: Number} <: AbstractTableau{DT}
    name::Symbol
    o::Int
    s::Int

    q::CoefficientsRK{DT}
    p::CoefficientsRK{DT}

    q̃::CoefficientsRK{DT}
    p̃::CoefficientsRK{DT}

    ω::Matrix{DT}
    d::Vector{DT}

    function TableauSLRK(name::Symbol, o::Int, s::Int,
                        q::CoefficientsRK{DT}, p::CoefficientsRK{DT},
                        q̃::CoefficientsRK{DT}, p̃::CoefficientsRK{DT},
                        ω::Matrix{DT}, d=DT[]) where {DT}

        @assert s > 0 "Number of stages s must be > 0"

        @assert s==q.s==p.s==q̃.s==p̃.s
        @assert size(ω,1)==s
        @assert size(ω,2)==s+1

        @assert length(d)==0 || length(d)==s

        new{DT}(name, o, s, q, p, q̃, p̃, ω, d)
    end
end


"Parameters for right-hand side function of Specialised Partitioned Additive Runge-Kutta methods for Variational systems."
mutable struct ParametersSLRK{DT,TT,D,S,θT,FT,GT,G̅T,ϕT,ψT,tabType} <: AbstractParametersSPARK{DT,TT}
    f_ϑ::θT
    f_f::FT
    f_g::GT
    f_g̅::G̅T
    f_ϕ::ϕT
    f_ψ::ψT

    Δt::TT

    tab::tabType

    @ParametersSPARK

    function ParametersSLRK{DT,D,S}(f_ϑ::θT, f_f::FT, f_g::GT, f_g̅::G̅T, f_ϕ::ϕT, f_ψ::ψT, Δt::TT, tab::tabType) where {DT,TT,D,S,θT,FT,GT,G̅T,ϕT,ψT,tabType}
        # create solution vectors
        q = zeros(DT,D)
        p = zeros(DT,D)
        λ = zeros(DT,D)

        new{DT,TT,D,S,θT,FT,GT,G̅T,ϕT,ψT,tabType}(f_ϑ, f_f, f_g, f_g̅, f_ϕ, f_ψ, Δt, tab, zero(TT), q, p, λ)
    end
end


@doc raw"""
Specialised Partitioned Additive Runge-Kutta integrator for degenerate
variational systems with projection on secondary constraint.

This integrator solves the following system of equations for the internal stages,
```math
\begin{align}
Q_{n,i} &= q_{n} + h \sum \limits_{j=1}^{s} a^1_{ij} V_{n,j} + h \sum \limits_{j=1}^{\sigma} a^2_{ij} \tilde{\Lambda}_{n,j} , & i &= 1, ..., s , \\
P_{n,i} &= p_{n} + h \sum \limits_{j=1}^{s} \bar{a}^1_{ij} F^1_{n,j} + h \sum \limits_{j=1}^{s} \bar{a}^2_{ij} F^2_{n,j} + h \sum \limits_{j=1}^{\sigma} \bar{a}^3_{ij} \tilde{F}^3_{n,j} , & i &= 1, ..., s , \\
0 &= \Phi_{n,i} , & i &= 1, ..., s ,
\end{align}
```
the projective stages
```math
\begin{align}
\tilde{Q}_{n,i} &= q_{n} + h \sum \limits_{j=1}^{s} \alpha^1_{ij} \tilde{V}_{n,j} + h \sum \limits_{j=1}^{\sigma} \alpha^2_{ij} \tilde{\Lambda}_{n,j} , & i &= 1, ..., \sigma , \\
\tilde{P}_{n,i} &= p_{n} + h \sum \limits_{j=1}^{s} \bar{\alpha}^1_{ij} F^1_{n,j} + h \sum \limits_{j=1}^{s} \bar{\alpha}^2_{ij} \tilde{F}^2_{n,j} + h \sum \limits_{j=1}^{\sigma} \bar{\alpha}^3_{ij} \tilde{F}^3_{n,j} , & i &= 1, ..., \sigma , \\
0 &= \tilde{\Phi}_{n,i} , & i &= 1, ..., \sigma , \\
0 &= \sum \limits_{j=1}^{\sigma} \omega_{ij} \tilde{\Psi}_{n,i} , & i &= 1, ..., \sigma-1 ,
\end{align}
```
and update rule
```math
\begin{align}
q_{n+1} &= q_{n} + h \sum \limits_{i=1}^{s} b^1_{i} V_{n,i} + h \sum \limits_{i=1}^{\sigma} b^2_{i} \tilde{\Lambda}_{n,i} , \\
p_{n+1} &= p_{n} + h \sum \limits_{i=1}^{s} b^1_{i} F^1_{n,i} + h \sum \limits_{i=1}^{s} b^2_{i} F^2_{n,i} + h \sum \limits_{i=1}^{\sigma} b^3_{i} \tilde{F}^3_{n,i} , \\
0 &= \phi (q_{n+1}, p_{n+1}) ,
\end{align}
```
with definitions
```math
\begin{align}
F^1_{n,i} &= \nabla H (Q_{n,i}) , & i &= 1, ..., s , \\
F^2_{n,i} &= \nabla \vartheta (Q_{n,i}) \cdot V_{n,i} , & i &= 1, ..., s , \\
\tilde{F}^3_{n,i} &= \nabla \phi (\tilde{Q}_{n,i}, \tilde{P}_{n,i}) \cdot \tilde{\Lambda}_{n,i} , & i &= 1, ..., \sigma , \\
\Phi_{n,i} &= \phi (Q_{n,i}, P_{n,i}) = P_{n,i} - \vartheta (Q_{n,i}) , & i &= 1, ..., s , \\
\tilde{\Phi}_{n,i} &= \phi (\tilde{Q}_{n,i}, \tilde{P}_{n,i}) = \tilde{P}_{n,i} - \vartheta (\tilde{Q}_{n,i}) , & i &= 1, ..., \sigma , \\
\tilde{\Psi}_{n,i} &= \psi (\tilde{Q}_{n,i}, \tilde{V}_{n,i}, \tilde{P}_{n,i}, \tilde{F}_{n,i}) = \tilde{F}_{n,i} - \tilde{V}_{n,i} \cdot \nabla \vartheta (\tilde{Q}_{n,i}) , & i &= 1, ..., \sigma ,
\end{align}
```
so that
```math
\begin{align}
F^1_{n,i} + F^2_{n,i} &= \frac{\partial L}{\partial q} (Q_{n,i}, V_{n,i}) , & i &= 1, ..., s , \\
\phi (Q_{n,i}, P_{n,i}) &= P_{n,i} - \frac{\partial L}{\partial v} (Q_{n,i}, V_{n,i}) , & i &= 1, ..., s , \\
\psi (\tilde{Q}_{n,i}, \tilde{V}_{n,i}, \tilde{P}_{n,i}, \tilde{F}_{n,i}) &= \tilde{F}_{n,i} - \tilde{V}_{n,i} \cdot \nabla \frac{\partial L}{\partial v} (\tilde{Q}_{n,i}, \tilde{V}_{n,i}) , & i &= 1, ..., \sigma .
\end{align}
```
"""
struct IntegratorSLRK{DT, TT, ET <: VDAE{DT,TT},
                              PT <: ParametersSLRK{DT,TT},
                              ST <: NonlinearSolver{DT},
                              IT <: InitialGuessPODE{DT,TT}, D, S} <: AbstractIntegratorVSPARK{DT, TT}
    equation::ET
    tableau::TableauSLRK{TT}

    params::PT
    solver::ST
    iguess::IT
    cache::IntegratorCacheSPARK{DT,TT,D,S,S}
end

function IntegratorSLRK(equation::VDAE{DT,TT},
                        tableau::TableauSLRK{TT}, Δt::TT) where {DT,TT}
    D = equation.d
    S = tableau.s

    # println(tableau)

    N = 4*D*S

    if isdefined(tableau, :d) && length(tableau.d) > 0
        N += D
    end

    # create params
    params = ParametersSLRK{DT,D,S}(equation.ϑ, equation.f, equation.g, equation.g̅, equation.ϕ, equation.ψ, Δt, tableau)

    # create solver
    solver = create_nonlinear_solver(DT, N, params)

    # create initial guess
    iguess = InitialGuessPODE(get_config(:ig_interpolation), equation, Δt)

    # create cache
    cache = IntegratorCacheSPARK{DT,TT,D,S,S}()

    # create integrator
    IntegratorSLRK{DT, TT, typeof(equation), typeof(params), typeof(solver), typeof(iguess), D, S}(
                              equation, tableau, params, solver, iguess, cache)
end

@inline pstages(int::IntegratorSLRK) = int.tableau.s
@inline Base.ndims(int::IntegratorSLRK{DT,TT,ET,PT,ST,IT,D,S}) where {DT,TT,ET,PT,ST,IT,D,S} = D


function compute_stages!(x::Vector{ST}, cache::IntegratorCacheSPARK{ST,TT,D,S},
                                        params::ParametersSLRK{DT,TT,D,S}) where {ST,DT,TT,D,S}
    local t::TT

    for i in 1:S
        for k in 1:D
            # copy y to Y, Z and Λ
            cache.Yp[i][k] = x[4*(D*(i-1)+k-1)+1]
            cache.Zp[i][k] = x[4*(D*(i-1)+k-1)+2]
            cache.Vp[i][k] = x[4*(D*(i-1)+k-1)+3]
            cache.Λp[i][k] = x[4*(D*(i-1)+k-1)+4]

            # compute Q and P
            cache.Qp[i][k] = params.q[k] + params.Δt * cache.Yp[i][k]
            cache.Pp[i][k] = params.p[k] + params.Δt * cache.Zp[i][k]
        end

        # compute f(X)
        t = params.t + params.Δt * params.tab.p.c[i]
        params.f_f(t, cache.Qp[i], cache.Vp[i], cache.Fp[i])
        params.f_g(t, cache.Qp[i], cache.Λp[i], cache.Gp[i])
        params.f_ϕ(t, cache.Qp[i], cache.Pp[i], cache.Φp[i])
        params.f_ψ(t, cache.Qp[i], cache.Pp[i], cache.Vp[i], cache.Fp[i], cache.Ψp[i])
    end

    if length(params.tab.d) > 0
        for k in 1:D
            cache.μ[k] = x[4*D*S+k]
        end
    end

    # compute q and p
    cache.q̃ .= params.q
    cache.p̃ .= params.p
    for i in 1:S
        cache.q̃ .+= params.Δt .* params.tab.q.b[i] .* cache.Vp[i]
        cache.q̃ .+= params.Δt .* params.tab.q̃.b[i] .* cache.Λp[i]
        cache.p̃ .+= params.Δt .* params.tab.p.b[i] .* cache.Fp[i]
        cache.p̃ .+= params.Δt .* params.tab.p̃.b[i] .* cache.Gp[i]
    end

    # compute ϕ(q,p)
    t = params.t + params.Δt
    params.f_ϕ(t, cache.q̃, cache.p̃, cache.ϕ̃)
end


"Compute stages of specialised partitioned additive Runge-Kutta methods for variational systems."
@generated function Integrators.function_stages!(y::Vector{ST}, b::Vector{ST}, params::ParametersSLRK{DT,TT,D,S}) where {ST,DT,TT,D,S}
    cache = IntegratorCacheSPARK{ST,TT,D,S,S}()

    quote
        compute_stages!(y, $cache, params)

        # compute b = - [(Y-AV-AU), (Z-AF-AG), Φ, ωΨ]
        for i in 1:S
            for k in 1:D
                b[4*(D*(i-1)+k-1)+1] = - $cache.Yp[i][k]
                b[4*(D*(i-1)+k-1)+2] = - $cache.Zp[i][k]
                b[4*(D*(i-1)+k-1)+3] = - $cache.Φp[i][k]
                b[4*(D*(i-1)+k-1)+4] = params.tab.ω[i,S+1] * $cache.ϕ̃[k]
                for j in 1:S
                    b[4*(D*(i-1)+k-1)+1] += params.tab.q.a[i,j] * $cache.Vp[j][k]
                    b[4*(D*(i-1)+k-1)+1] += params.tab.q̃.a[i,j] * $cache.Λp[j][k]
                    b[4*(D*(i-1)+k-1)+2] += params.tab.p.a[i,j] * $cache.Fp[j][k]
                    b[4*(D*(i-1)+k-1)+2] += params.tab.p̃.a[i,j] * $cache.Gp[j][k]
                    b[4*(D*(i-1)+k-1)+4] += params.tab.ω[i,j]   * $cache.Ψp[j][k]
                end
            end
        end

        if length(params.tab.d) > 0
            for i in 1:S
                for k in 1:D
                    b[4*(D*(i-1)+k-1)+2] += $cache.μ[k] * params.tab.d[i] / params.tab.p.b[i]
                    b[4*(D*(i-1)+k-1)+3] += $cache.μ[k] * params.tab.d[i] / params.tab.p.b[i]
                end
            end

            for k in 1:D
                b[4*D*S+k] = 0
                for i in 1:S
                    b[4*D*S+k] -= $cache.Vp[i][k] * params.tab.d[i]
                end
            end
        end
    end
end


function initial_guess!(int::IntegratorSLRK, sol::AtomicSolutionPDAE)
    for i in 1:nstages(int)
        evaluate!(int.iguess, sol.q, sol.p, sol.v, sol.f,
                              sol.q̅, sol.p̅, sol.v̅, sol.f̅,
                              int.cache.q̃, int.cache.p̃, int.cache.ṽ, int.cache.f̃,
                              tableau(int).q̃.c[i], tableau(int).p̃.c[i])

        for k in 1:ndims(int)
            int.solver.x[4*(ndims(int)*(i-1)+k-1)+1] = (int.cache.q̃[k] - sol.q[k])/timestep(int)
            int.solver.x[4*(ndims(int)*(i-1)+k-1)+2] = (int.cache.p̃[k] - sol.p[k])/timestep(int)
            int.solver.x[4*(ndims(int)*(i-1)+k-1)+3] = int.cache.ṽ[k]
            int.solver.x[4*(ndims(int)*(i-1)+k-1)+4] = 0
        end
    end

    if isdefined(tableau(int), :d) && length(tableau(int).d) > 0
        for k in 1:ndims(int)
            int.solver.x[4*ndims(int)*pstages(int)+k] = 0
        end
    end
end


function update_solution!(int::IntegratorSLRK{DT,TT}, sol::AtomicSolutionPDAE{DT,TT}) where {DT,TT}
    # compute final update
    update_solution!(sol.q, int.cache.Vp, int.params.tab.q.b, timestep(int))
    update_solution!(sol.q, int.cache.Λp, int.params.tab.q̃.b, timestep(int))
    update_solution!(sol.p, int.cache.Fp, int.params.tab.p.b, timestep(int))
    update_solution!(sol.p, int.cache.Gp, int.params.tab.p̃.b, timestep(int))
end
