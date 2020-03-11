"Holds the tableau of an variational partitioned additive Runge-Kutta method."
struct TableauVPARK{T} <: AbstractTableau{T}
    name::Symbol
    o::Int
    s::Int
    r::Int

    q::CoefficientsARK{T}
    p::CoefficientsARK{T}

    q̃::CoefficientsPRK{T}
    p̃::CoefficientsPRK{T}

    λ::CoefficientsMRK{T}

    d::Vector{T}

    function TableauVPARK{T}(name, o, s, r, q, p, q̃, p̃, λ, d) where {T}
        @assert isa(name, Symbol)
        @assert isa(s, Integer)
        @assert isa(r, Integer)
        @assert isa(o, Integer)

        @assert s > 0 "Number of stages s must be > 0"
        @assert r > 0 "Number of stages r must be > 0"

        @assert s==q.s==p.s==q̃.s==p̃.s==length(d)
        @assert r==q.r==p.r==q̃.r==p̃.r==λ.r

        new(name, o, s, r, q, p, q̃, p̃, λ, d)
    end

    function TableauVPARK{T}(name, o, s, r, q, p, q̃, p̃, λ) where {T}
        @assert isa(name, Symbol)
        @assert isa(s, Integer)
        @assert isa(r, Integer)
        @assert isa(o, Integer)

        @assert s > 0 "Number of stages s must be > 0"
        @assert r > 0 "Number of stages r must be > 0"

        @assert s==q.s==p.s==q̃.s==p̃.s
        @assert r==q.r==p.r==q̃.r==p̃.r==λ.r

        new(name, o, s, r, q, p, q̃, p̃, λ)
    end
end

function TableauVPARK(name::Symbol, order::Int,
                        a_q::Matrix{T}, a_p::Matrix{T},
                        α_q::Matrix{T}, α_p::Matrix{T},
                        a_q̃::Matrix{T}, a_p̃::Matrix{T},
                        α_q̃::Matrix{T}, α_p̃::Matrix{T},
                        b_q::Vector{T}, b_p::Vector{T},
                        β_q::Vector{T}, β_p::Vector{T},
                        c_q::Vector{T}, c_p::Vector{T},
                        c_λ::Vector{T}, d_λ::Vector{T},
                        d::Vector{T}) where {T <: Real}

    s = length(c_q)
    r = length(c_λ)

    @assert s > 0 "Number of stages s must be > 0"
    @assert r > 0 "Number of stages r must be > 0"

    @assert s==size(a_q,1)==size(a_q,2)==length(b_q)==length(c_q)
    @assert s==size(a_p,1)==size(a_p,2)==length(b_p)==length(c_p)
    @assert s==size(α_q,1)==size(α_p,1)
    @assert r==size(α_q,2)==size(α_p,2)
    @assert s==length(d)
    @assert r==length(c_λ)==length(d_λ)
    @assert r==size(a_q̃,1)==size(a_p̃,1)
    @assert s==size(a_q̃,2)==size(a_p̃,2)
    @assert r==size(α_q̃,1)==size(α_q̃,2)==length(β_q)
    @assert r==size(α_p̃,1)==size(α_p̃,2)==length(β_p)

    q = CoefficientsARK{T}(name, order, s, r, a_q, b_q, c_q, α_q, β_q)
    p = CoefficientsARK{T}(name, order, s, r, a_p, b_p, c_p, α_p, β_p)
    q̃ = CoefficientsPRK{T}(name, order, s, r, a_q̃, c_λ, α_q̃)
    p̃ = CoefficientsPRK{T}(name, order, s, r, a_p̃, c_λ, α_p̃)
    λ = CoefficientsMRK{T}(name, r, d_λ, c_λ)

    TableauVPARK{T}(name, order, s, r, q, p, q̃, p̃, λ, d)
end


function TableauVPARK(name::Symbol, order::Int,
                        a_q::Matrix{T}, a_p::Matrix{T},
                        α_q::Matrix{T}, α_p::Matrix{T},
                        a_q̃::Matrix{T}, a_p̃::Matrix{T},
                        α_q̃::Matrix{T}, α_p̃::Matrix{T},
                        b_q::Vector{T}, b_p::Vector{T},
                        β_q::Vector{T}, β_p::Vector{T},
                        c_q::Vector{T}, c_p::Vector{T},
                        c_λ::Vector{T}, d_λ::Vector{T}) where {T <: Real}

    s = length(c_q)
    r = length(c_λ)

    @assert s > 0 "Number of stages s must be > 0"
    @assert r > 0 "Number of stages r must be > 0"

    @assert s==size(a_q,1)==size(a_q,2)==length(b_q)==length(c_q)
    @assert s==size(a_p,1)==size(a_p,2)==length(b_p)==length(c_p)
    @assert s==size(α_q,1)==size(α_p,1)
    @assert r==size(α_q,2)==size(α_p,2)
    @assert r==length(c_λ)==length(d_λ)
    @assert r==size(a_q̃,1)==size(a_p̃,1)
    @assert s==size(a_q̃,2)==size(a_p̃,2)
    @assert r==size(α_q̃,1)==size(α_q̃,2)==length(β_q)
    @assert r==size(α_p̃,1)==size(α_p̃,2)==length(β_p)

    q = CoefficientsARK{T}(name, order, s, r, a_q, b_q, c_q, α_q, β_q)
    p = CoefficientsARK{T}(name, order, s, r, a_p, b_p, c_p, α_p, β_p)
    q̃ = CoefficientsPRK{T}(name, order, s, r, a_q̃, c_λ, α_q̃)
    p̃ = CoefficientsPRK{T}(name, order, s, r, a_p̃, c_λ, α_p̃)
    λ = CoefficientsMRK{T}(name, r, d_λ, c_λ)

    TableauVPARK{T}(name, order, s, r, q, p, q̃, p̃, λ)
end

# TODO function readTableauVPARKFromFile(dir::AbstractString, name::AbstractString)


"Parameters for right-hand side function of variational partitioned additive Runge-Kutta methods."
mutable struct ParametersVPARK{DT, TT, D, S, R, ET <: NamedTuple} <: AbstractParametersSPARK{DT,TT}
    equs::ET
    tab::TableauVPARK{TT}
    Δt::TT

    @ParametersSPARK

    function ParametersVPARK{DT,D}(equs::ET, tab::TableauVPARK{TT}, Δt::TT) where {DT,TT,D,S,R,ET <: NamedTuple}
        # create solution vectors
        q = zeros(DT,D)
        p = zeros(DT,D)
        λ = zeros(DT,D)

        new{DT,TT,D,tab.s,tab.r,ET}(equs, tab, Δt, zero(TT), q, p, λ)
    end
end


@doc raw"""
Variational partitioned additive Runge-Kutta integrator.

This integrator solves the following system of equations for the internal stages,
```math
\begin{align}
Q_{n,i} &= q_{n} + h \sum \limits_{j=1}^{s} a_{ij} V_{n,j} + h \sum \limits_{j=1}^{r} \alpha_{ij} U_{n,j} , & i &= 1, ..., s , \\
P_{n,i} &= p_{n} + h \sum \limits_{j=1}^{s} a_{ij} F_{n,j} + h \sum \limits_{j=1}^{r} \alpha_{ij} G_{n,j} , & i &= 1, ..., s , \\
\tilde{Q}_{n,i} &= q_{n} + h \sum \limits_{j=1}^{s} \tilde{a}_{ij} V_{n,j} + h \sum \limits_{j=1}^{r} \tilde{\alpha}_{ij} U_{n,j} , & i &= 1, ..., r , \\
\tilde{P}_{n,i} &= p_{n} + h \sum \limits_{j=1}^{s} \tilde{a}_{ij} F_{n,j} + h \sum \limits_{j=1}^{r} \tilde{\alpha}_{ij} G_{n,j} , & i &= 1, ..., r , \\
\tilde{\Phi}_{n,i} &= 0 , & i &= 1, ..., r ,
\end{align}
```
with definitions
```math
\begin{align}
P_{n,i} &= \frac{\partial L}{\partial v} (Q_{n,i}, V_{n,i}) , & i &= 1, ..., s , \\
F_{n,i} &= \frac{\partial L}{\partial q} (Q_{n,i}, V_{n,i}) , & i &= 1, ..., s , \\
U_{n,i} &= \hphantom{-} \frac{\partial \phi}{\partial p} (\tilde{Q}_{n,i}, \tilde{P}_{n,i})^{T} \Lambda_{n,i} , & i &= 1, ..., r , \\
G_{n,i} &=           -  \frac{\partial \phi}{\partial q} (\tilde{Q}_{n,i}, \tilde{P}_{n,i})^{T} \Lambda_{n,i} , & i &= 1, ..., r , \\
\tilde{\Phi}_{n,i} &= \phi(\tilde{Q}_{n,i}, \tilde{P}_{n,i}) , & i &= 1, ..., r ,
\end{align}
```
and update rule
```math
\begin{align}
q_{n+1} &= q_{n} + h \sum \limits_{i=1}^{s} b_{i} V_{n,i} + h \sum \limits_{i=1}^{r} \beta_{i} U_{n,i} , \\
p_{n+1} &= p_{n} + h \sum \limits_{i=1}^{s} b_{i} F_{n,i} + h \sum \limits_{i=1}^{r} \beta_{i} G_{n,i} .
\end{align}
```
"""
struct IntegratorVPARK{DT, TT, D, S, R, PT <: ParametersVPARK{DT,TT},
                                        ST <: NonlinearSolver{DT},
                                        IT <: InitialGuessPODE{DT,TT}} <: AbstractIntegratorVSPARK{DT,TT}
    params::PT
    solver::ST
    iguess::IT
    cache::IntegratorCacheSPARK{DT, TT, D, S, R}

    function IntegratorVPARK(params::ParametersVPARK{DT,TT,D,S,R}, solver::ST, iguess::IT, cache) where {DT,TT,D,S,R,ST,IT}
        new{DT, TT, D, S, R, typeof(params), ST, IT}(params, solver, iguess, cache)
    end

    function IntegratorVPARK{DT,D}(equations::NamedTuple, tableau::TableauVPARK{TT}, Δt::TT) where {DT,TT,D}
        # get number of stages
        S = tableau.s
        R = tableau.r

        N = 3*D*S + 3*D*R

        if isdefined(tableau, :d)
            N += D
        end

        # create params
        params = ParametersVPARK{DT,D}(equations, tableau, Δt)

        # create solver
        solver = create_nonlinear_solver(DT, N, params)

        # create initial guess
        iguess = InitialGuessPODE{DT,D}(get_config(:ig_interpolation), equations[:v], equations[:f], Δt)

        # create cache
        cache = IntegratorCacheSPARK{DT, TT, D, S, R}()

        # create integrator
        IntegratorVPARK(params, solver, iguess, cache)
    end

    function IntegratorVPARK(equation::IDAE{DT,TT}, tableau::TableauVPARK{TT}, Δt::TT; kwargs...) where {DT,TT}
        IntegratorVPARK{DT, equation.d}(get_function_tuple(equation), tableau, Δt; kwargs...)
    end
end


@inline equation(int::IntegratorVPARK, i::Symbol) = int.params.equs[i]
@inline equations(int::IntegratorVPARK) = int.params.equs
@inline tableau(int::IntegratorVPARK) = int.params.tab
@inline nstages(int::IntegratorVPARK{DT,TT,D,S,R}) where {DT,TT,D,S,R} = S
@inline pstages(int::IntegratorVPARK{DT,TT,D,S,R}) where {DT,TT,D,S,R} = R
@inline Base.ndims(int::IntegratorVPARK{DT,TT,D,S,R}) where {DT,TT,D,S,R} = D


function Integrators.initialize!(int::IntegratorVPARK, sol::AtomicSolutionPDAE)
    sol.t̅ = sol.t - timestep(int)

    equation(int, :v)(sol.t, sol.q, sol.p, sol.v)
    equation(int, :f)(sol.t, sol.q, sol.p, sol.f)

    initialize!(int.iguess, sol.t, sol.q, sol.p, sol.v, sol.f,
                            sol.t̅, sol.q̅, sol.p̅, sol.v̅, sol.f̅)
end


function compute_stages!(x::Vector{ST}, cache::IntegratorCacheSPARK{ST,TT,D,S,R},
                                        params::ParametersVPARK{DT,TT,D,S,R}) where {ST,DT,TT,D,S,R}
    local tpᵢ::TT
    local tλᵢ::TT

    for i in 1:S
        for k in 1:D
            # copy x to Y, Z
            cache.Yi[i][k] = x[3*(D*(i-1)+k-1)+1]
            cache.Zi[i][k] = x[3*(D*(i-1)+k-1)+2]
            cache.Vi[i][k] = x[3*(D*(i-1)+k-1)+3]

            # compute Q and P
            cache.Qi[i][k] = params.q[k] + params.Δt * cache.Yi[i][k]
            cache.Pi[i][k] = params.p[k] + params.Δt * cache.Zi[i][k]
        end

        # compute f(X)
        tpᵢ = params.t + params.Δt * params.tab.p.c[i]
        params.equs[:ϑ](tpᵢ, cache.Qi[i], cache.Vi[i], cache.Φi[i])
        params.equs[:f](tpᵢ, cache.Qi[i], cache.Vi[i], cache.Fi[i])

        cache.Φi[i] .-= cache.Pi[i]
    end


    for i in 1:R
        for k in 1:D
            # copy y to Y, Z and Λ
            cache.Yp[i][k] = x[3*D*S+3*(D*(i-1)+k-1)+1]
            cache.Zp[i][k] = x[3*D*S+3*(D*(i-1)+k-1)+2]
            cache.Λp[i][k] = x[3*D*S+3*(D*(i-1)+k-1)+3]

            # compute Q and V
            cache.Qp[i][k] = params.q[k] + params.Δt * cache.Yp[i][k]
            cache.Pp[i][k] = params.p[k] + params.Δt * cache.Zp[i][k]
        end

        # compute f(X)
        tλᵢ = params.t + params.Δt * params.tab.λ.c[i]
        params.equs[:u](tλᵢ, cache.Qp[i], cache.Pp[i], cache.Λp[i], cache.Up[i])
        params.equs[:g](tλᵢ, cache.Qp[i], cache.Pp[i], cache.Λp[i], cache.Gp[i])
        params.equs[:ϕ](tλᵢ, cache.Qp[i], cache.Pp[i], cache.Φp[i])
    end

    if isdefined(params.tab, :d) && length(params.tab.d) > 0
        for k in 1:D
            cache.μ[k] = x[3*D*S+3*D*R+k]
        end
    end
end


"Compute stages of variational partitioned additive Runge-Kutta methods."
@generated function Integrators.function_stages!(y::Vector{ST}, b::Vector{ST}, params::ParametersVPARK{DT,TT,D,S,R}) where {ST,DT,TT,D,S,R}
    cache = IntegratorCacheSPARK{ST,TT,D,S,R}()

    quote
        compute_stages!(y, $cache, params)

        # compute b = - [(Y-AV-AU), (Z-AF-AG), Φ]
        for i in 1:S
            for k in 1:D
                b[3*(D*(i-1)+k-1)+1] = - $cache.Yi[i][k]
                b[3*(D*(i-1)+k-1)+2] = - $cache.Zi[i][k]
                b[3*(D*(i-1)+k-1)+3] = - $cache.Φi[i][k]
                for j in 1:S
                    b[3*(D*(i-1)+k-1)+1] += params.tab.q.a[i,j] * $cache.Vi[j][k]
                    b[3*(D*(i-1)+k-1)+2] += params.tab.p.a[i,j] * $cache.Fi[j][k]
                end
                for j in 1:R
                    b[3*(D*(i-1)+k-1)+1] += params.tab.q.α[i,j] * $cache.Up[j][k]
                    b[3*(D*(i-1)+k-1)+2] += params.tab.p.α[i,j] * $cache.Gp[j][k]
                end
            end
        end

        # compute b = - [(Y-AV-AU), (Z-AF-AG), Φ]
        for i in 1:R
            for k in 1:D
                b[3*D*S+3*(D*(i-1)+k-1)+1] = - $cache.Yp[i][k]
                b[3*D*S+3*(D*(i-1)+k-1)+2] = - $cache.Zp[i][k]
                b[3*D*S+3*(D*(i-1)+k-1)+3] = - $cache.Φp[i][k]
                for j in 1:S
                    b[3*D*S+3*(D*(i-1)+k-1)+1] += params.tab.q̃.a[i,j] * $cache.Vi[j][k]
                    b[3*D*S+3*(D*(i-1)+k-1)+2] += params.tab.p̃.a[i,j] * $cache.Fi[j][k]
                end
                for j in 1:R
                    b[3*D*S+3*(D*(i-1)+k-1)+1] += params.tab.q̃.α[i,j] * $cache.Up[j][k]
                    b[3*D*S+3*(D*(i-1)+k-1)+2] += params.tab.p̃.α[i,j] * $cache.Gp[j][k]
                end
            end
        end

        # compute b = - [Λ₁-λ]
        if params.tab.λ.c[1] == 0
            for k in 1:D
                b[3*D*S+3*(k-1)+3] = - $cache.Λp[1][k] + params.λ[k]
            end
        end

        if isdefined(params.tab, :d) && length(params.tab.d) > 0
            for i in 1:S
                for k in 1:D
                    b[3*(D*(i-1)+k-1)+3] -= $cache.μ[k] * params.tab.d[i]
                end
            end

            for k in 1:D
                b[3*D*S+3*D*R+k] = 0
                for i in 1:S
                    b[3*D*S+3*D*R+k] -= $cache.Vi[i][k] * params.tab.d[i]
                end
            end
        end
    end
end


function initial_guess!(int::IntegratorVPARK, sol::AtomicSolutionPDAE)
    for i in eachstage(int)
        evaluate!(int.iguess, sol.q, sol.p, sol.v, sol.f,
                              sol.q̅, sol.p̅, sol.v̅, sol.f̅,
                              int.cache.q̃, int.cache.p̃, int.cache.ṽ, int.cache.f̃,
                              tableau(int).q.c[i], tableau(int).p.c[i])

        for k in eachdim(int)
            int.solver.x[3*(ndims(int)*(i-1)+k-1)+1] = (int.cache.q̃[k] - sol.q[k])/timestep(int)
            int.solver.x[3*(ndims(int)*(i-1)+k-1)+2] = (int.cache.p̃[k] - sol.p[k])/timestep(int)
            int.solver.x[3*(ndims(int)*(i-1)+k-1)+3] = int.cache.ṽ[k]
        end
    end

    for i in 1:pstages(int)
        evaluate!(int.iguess, sol.q, sol.p, sol.v, sol.f,
                              sol.q̅, sol.p̅, sol.v̅, sol.f̅,
                              int.cache.q̃, int.cache.p̃, int.cache.ṽ, int.cache.f̃,
                              tableau(int).q̃.c[i], tableau(int).p̃.c[i])

        for k in 1:ndims(int)
            int.solver.x[3*ndims(int)*nstages(int)+3*(ndims(int)*(i-1)+k-1)+1] = (int.cache.q̃[k] - sol.q[k])/timestep(int)
            int.solver.x[3*ndims(int)*nstages(int)+3*(ndims(int)*(i-1)+k-1)+2] = (int.cache.p̃[k] - sol.p[k])/timestep(int)
            int.solver.x[3*ndims(int)*nstages(int)+3*(ndims(int)*(i-1)+k-1)+3] = 0
        end
    end

    if int.params.tab.λ.c[1] == 0
        for k in eachdim(int)
            int.solver.x[3*ndims(int)*nstages(int)+3*(k-1)+3] = sol.λ[k]
        end
    end

    if isdefined(tableau(int), :d)
        for k in eachdim(int)
            int.solver.x[3*ndims(int)*nstages(int)+3*ndims(int)*pstages(int)+k] = 0
        end
    end
end


function update_solution!(int::IntegratorVPARK{DT,TT}, sol::AtomicSolutionPDAE{DT,TT}) where {DT,TT}
    # compute final update
    update_solution!(sol.q, sol.q̃, int.cache.Vi, int.params.tab.q.b, timestep(int))
    update_solution!(sol.p, sol.p̃, int.cache.Fi, int.params.tab.p.b, timestep(int))

    # compute projection
    update_solution!(sol.q, sol.q̃, int.cache.Up, int.params.tab.q.β, timestep(int))
    update_solution!(sol.p, sol.p̃, int.cache.Gp, int.params.tab.p.β, timestep(int))
    update_multiplier!(sol.λ, int.cache.Λp, int.params.tab.λ.b)
end
