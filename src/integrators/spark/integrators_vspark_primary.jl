
"Parameters for right-hand side function of Specialised Partitioned Additive Runge-Kutta methods for Variational systems."
mutable struct ParametersVSPARKprimary{DT,TT,D,S,R,FT,PT,UT,GT,ϕT} <: Parameters{DT,TT}
    f_f::FT
    f_p::PT
    f_u::UT
    f_g::GT
    f_ϕ::ϕT

    Δt::TT

    t_q::CoefficientsARK{TT}
    t_p::CoefficientsARK{TT}
    t_q̃::CoefficientsPRK{TT}
    t_p̃::CoefficientsPRK{TT}
    t_λ::CoefficientsMRK{TT}
    t_ω::Matrix{TT}
    d_v::Vector{TT}

    t::TT
    q::Vector{DT}
    p::Vector{DT}
    λ::Vector{DT}

    function ParametersVSPARKprimary{DT,TT,D,S,R,FT,PT,UT,GT,ϕT}(f_f, f_p, f_u, f_g, f_ϕ, Δt, t_q, t_p, t_q̃, t_p̃, t_λ, t_ω, d_v) where {DT,TT,D,S,R,FT,PT,UT,GT,ϕT}
        # create solution vectors
        q = zeros(DT,D)
        p = zeros(DT,D)
        λ = zeros(DT,D)

        new(f_f, f_p, f_u, f_g, f_ϕ, Δt,
            t_q, t_p, t_q̃, t_p̃, t_λ, t_ω, d_v,
            zero(TT), q, p, λ)
    end
end


function compute_stages!(x::Vector{ST}, cache::IntegratorCacheVSPARK{ST,TT,D,S,R},
                                        params::ParametersVSPARKprimary{DT,TT,D,S,R}) where {ST,DT,TT,D,S,R}
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
        tpᵢ = params.t + params.Δt * params.t_p.c[i]
        params.f_f(tpᵢ, cache.Qi[i], cache.Vi[i], cache.Fi[i])
        params.f_p(tpᵢ, cache.Qi[i], cache.Vi[i], cache.Φi[i])

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
        tλᵢ = params.t + params.Δt * params.t_λ.c[i]
        params.f_u(tλᵢ, cache.Qp[i], cache.Pp[i], cache.Λp[i], cache.Up[i])
        params.f_g(tλᵢ, cache.Qp[i], cache.Pp[i], cache.Λp[i], cache.Gp[i])
        params.f_ϕ(tλᵢ, cache.Qp[i], cache.Pp[i], cache.Φp[i])
    end

    if length(params.d_v) > 0
        for k in 1:D
            cache.μ[k] = x[3*D*S+3*D*R+k]
        end
    end
end


"Compute stages of specialised partitioned additive Runge-Kutta methods for variational systems."
@generated function function_stages!(y::Vector{ST}, b::Vector{ST}, params::ParametersVSPARKprimary{DT,TT,D,S,R}) where {ST,DT,TT,D,S,R}
    cache = IntegratorCacheVSPARK{ST,TT,D,S,R}()

    quote
        compute_stages!(y, $cache, params)

        # compute b = - [(Y-AV-AU), (Z-AF-AG), Φ]
        for i in 1:S
            for k in 1:D
                b[3*(D*(i-1)+k-1)+1] = - $cache.Yi[i][k]
                b[3*(D*(i-1)+k-1)+2] = - $cache.Zi[i][k]
                b[3*(D*(i-1)+k-1)+3] = - $cache.Φi[i][k]
                for j in 1:S
                    b[3*(D*(i-1)+k-1)+1] += params.t_q.a[i,j] * $cache.Vi[j][k]
                    b[3*(D*(i-1)+k-1)+2] += params.t_p.a[i,j] * $cache.Fi[j][k]
                end
                for j in 1:R
                    b[3*(D*(i-1)+k-1)+1] += params.t_q.α[i,j] * $cache.Up[j][k]
                    b[3*(D*(i-1)+k-1)+2] += params.t_p.α[i,j] * $cache.Gp[j][k]
                end
            end
        end

        # compute b = - [(Y-AV-AU), (Z-AF-AG), ωΦ]
        for i in 1:R
            for k in 1:D
                b[3*D*S+3*(D*(i-1)+k-1)+1] = - $cache.Yp[i][k]
                b[3*D*S+3*(D*(i-1)+k-1)+2] = - $cache.Zp[i][k]
                for j in 1:S
                    b[3*D*S+3*(D*(i-1)+k-1)+1] += params.t_q̃.a[i,j] * $cache.Vi[j][k]
                    b[3*D*S+3*(D*(i-1)+k-1)+2] += params.t_p̃.a[i,j] * $cache.Fi[j][k]
                end
                for j in 1:R
                    b[3*D*S+3*(D*(i-1)+k-1)+1] += params.t_q̃.α[i,j] * $cache.Up[j][k]
                    b[3*D*S+3*(D*(i-1)+k-1)+2] += params.t_p̃.α[i,j] * $cache.Gp[j][k]
                end
            end
        end
        for i in 1:R-1
            for k in 1:D
                b[3*D*S+3*(D*(i-1)+k-1)+3] = 0
                for j in 1:R
                    b[3*D*S+3*(D*(i-1)+k-1)+3] -= params.t_ω[i,j] * $cache.Φp[j][k]
                end
            end
        end

        # compute b = d_λ ⋅ Λ
        for k in 1:D
            b[3*D*S+3*(D*(R-1)+k-1)+3] = 0
            for j in 1:R
                b[3*D*S+3*(D*(R-1)+k-1)+3] -= params.t_λ.b[j] * $cache.Λp[j][k]
            end
        end

        if length(params.d_v) > 0
            for i in 1:S
                for k in 1:D
                    b[3*(D*(i-1)+k-1)+3] -= $cache.μ[k] * params.d_v[i]
                end
            end

            for k in 1:D
                b[3*D*S+3*D*R+k] = 0
                for i in 1:S
                    b[3*D*S+3*D*R+k] -= $cache.Vi[i][k] * params.d_v[i]
                end
            end
        end
    end
end


@doc raw"""
Specialised Partitioned Additive Runge-Kutta integrator for Variational systems.

This integrator solves the following system of equations for the internal stages,
```math
\begin{align}
Q_{n,i} &= q_{n} + h \sum \limits_{j=1}^{s} a_{ij} V_{n,j} + h \sum \limits_{j=1}^{r} \alpha_{ij} U_{n,j} , & i &= 1, ..., s , \\
P_{n,i} &= p_{n} + h \sum \limits_{j=1}^{s} a_{ij} F_{n,j} + h \sum \limits_{j=1}^{r} \alpha_{ij} G_{n,j} , & i &= 1, ..., s , \\
\tilde{Q}_{n,i} &= q_{n} + h \sum \limits_{j=1}^{s} \tilde{a}_{ij} V_{n,j} + h \sum \limits_{j=1}^{r} \tilde{\alpha}_{ij} U_{n,j} , & i &= 1, ..., r , \\
\tilde{P}_{n,i} &= p_{n} + h \sum \limits_{j=1}^{s} \tilde{a}_{ij} F_{n,j} + h \sum \limits_{j=1}^{r} \tilde{\alpha}_{ij} G_{n,j} , & i &= 1, ..., r , \\
0 &= \sum \limits_{j=1}^{r} \omega_{ij} \tilde{\Phi}_{n,j} , & i &= 1, ..., r-1 , \\
0 &= \sum \limits_{i=1}^{r} \tilde{d}_i \, \Lambda_{n,i} ,
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
struct IntegratorVSPARKprimary{DT, TT, ET <: IDAE{DT,TT},
                                       PT <: ParametersVSPARKprimary{DT,TT},
                                       ST <: NonlinearSolver{DT},
                                       IT <: InitialGuessPODE{DT,TT}} <: AbstractIntegratorSPARK{DT, TT}
    equation::ET
    tableau::TableauVSPARK{TT}

    params::PT
    solver::ST
    iguess::IT
end

function IntegratorVSPARKprimary(equation::IDAE{DT,TT,FT,PT,UT,GT,ϕT,VT},
                                 tableau::TableauVSPARK{TT}, Δt::TT) where {DT,TT,FT,PT,UT,GT,ϕT,VT}
    D = equation.d
    S = tableau.s
    R = tableau.r

    @assert tableau.ρ == tableau.r-1

    N = 3*D*S + 3*D*R

    if isdefined(tableau, :d)
        N += D
        d_v = tableau.d
    else
        d_v = DT[]
    end

    # create params
    params = ParametersVSPARKprimary{DT,TT,D,S,R,FT,PT,UT,GT,ϕT}(
                                equation.f, equation.p, equation.u, equation.g, equation.ϕ, Δt,
                                tableau.q, tableau.p, tableau.q̃, tableau.p̃, tableau.λ, tableau.ω, d_v)

    # create solver
    solver = create_nonlinear_solver(DT, N, params)

    # create initial guess
    iguess = InitialGuessPODE(get_config(:ig_interpolation), equation, Δt)

    # create integrator
    IntegratorVSPARKprimary{DT, TT, typeof(equation), typeof(params), typeof(solver), typeof(iguess)}(
                                        equation, tableau, params, solver, iguess)
end

equation(int::IntegratorVSPARKprimary) = int.equation
timestep(int::IntegratorVSPARKprimary) = int.params.Δt
tableau(int::IntegratorVSPARKprimary) = int.tableau
nstages(int::IntegratorVSPARKprimary) = int.tableau.s
pstages(int::IntegratorVSPARKprimary) = int.tableau.r


function create_integrator_cache(int::IntegratorVSPARKprimary{DT,TT}) where {DT,TT}
    IntegratorCacheVSPARK{DT, TT, ndims(int), nstages(int), pstages(int)}()
end


function update_params!(params::ParametersVSPARKprimary, cache::IntegratorCacheVSPARK)
    # set time for nonlinear solver and copy previous solution
    params.t  = cache.t
    params.q .= cache.q
    params.p .= cache.p
    params.λ .= cache.λ
end


function initialize!(int::IntegratorVSPARKprimary, cache::IntegratorCacheVSPARK)
    cache.t̅ = cache.t - timestep(int)

    equation(int).v(cache.t, cache.q, cache.p, cache.v)
    equation(int).f(cache.t, cache.q, cache.p, cache.f)

    initialize!(int.iguess, cache.t, cache.q, cache.p, cache.v, cache.f,
                            cache.t̅, cache.q̅, cache.p̅, cache.v̅, cache.f̅)
end


function initial_guess!(int::IntegratorVSPARKprimary, cache::IntegratorCacheVSPARK)
    for i in 1:nstages(int)
        evaluate!(int.iguess, cache.q, cache.p, cache.v, cache.f,
                              cache.q̅, cache.p̅, cache.v̅, cache.f̅,
                              cache.q̃, cache.p̃, cache.ṽ, cache.f̃,
                              tableau(int).q.c[i], tableau(int).p.c[i])

        for k in 1:ndims(int)
            int.solver.x[3*(ndims(int)*(i-1)+k-1)+1] = (cache.q̃[k] - cache.q[k])/timestep(int)
            int.solver.x[3*(ndims(int)*(i-1)+k-1)+2] = (cache.p̃[k] - cache.p[k])/timestep(int)
            int.solver.x[3*(ndims(int)*(i-1)+k-1)+3] = cache.ṽ[k]
        end
    end

    for i in 1:pstages(int)
        evaluate!(int.iguess, cache.q, cache.p, cache.v, cache.f,
                              cache.q̅, cache.p̅, cache.v̅, cache.f̅,
                              cache.q̃, cache.p̃, cache.ṽ, cache.f̃,
                              tableau(int).q̃.c[i], tableau(int).p̃.c[i])

        for k in 1:ndims(int)
            int.solver.x[3*ndims(int)*nstages(int)+3*(ndims(int)*(i-1)+k-1)+1] = (cache.q̃[k] - cache.q[k])/timestep(int)
            int.solver.x[3*ndims(int)*nstages(int)+3*(ndims(int)*(i-1)+k-1)+2] = (cache.p̃[k] - cache.p[k])/timestep(int)
            int.solver.x[3*ndims(int)*nstages(int)+3*(ndims(int)*(i-1)+k-1)+3] = 0
        end
    end

    if int.params.t_λ.c[1] == 0
        for k in 1:ndims(int)
            int.solver.x[3*ndims(int)*nstages(int)+3*(k-1)+3] = cache.λ[k]
        end
    end

    if isdefined(tableau(int), :d)
        for k in 1:ndims(int)
            int.solver.x[3*ndims(int)*nstages(int)+3*ndims(int)*pstages(int)+k] = 0
        end
    end
end


"Integrate an implicit DAE with a partitioned additive Runge-Kutta integrator for variational systems."
function integrate_step!(int::IntegratorVSPARKprimary{DT,TT}, cache::IntegratorCacheVSPARK{DT,TT}) where {DT,TT}
    # update nonlinear solver parameters from cache
    update_params!(int.params, cache)

    # compute initial guess
    initial_guess!(int, cache)

    # reset cache
    reset!(cache, timestep(int))

    # call nonlinear solver
    solve!(int.solver)

    # print solver status
    print_solver_status(int.solver.status, int.solver.params, cache.n)

    # check if solution contains NaNs or error bounds are violated
    check_solver_status(int.solver.status, int.solver.params, cache.n)

    # compute vector fields at internal stages
    compute_stages!(int.solver.x, cache, int.params)

    # compute final update
    update_solution!(cache.q, cache.Vi, int.params.t_q.b, timestep(int))
    update_solution!(cache.p, cache.Fi, int.params.t_p.b, timestep(int))

    # compute projection
    update_solution!(cache.q, cache.Up, int.params.t_q.β, timestep(int))
    update_solution!(cache.p, cache.Gp, int.params.t_p.β, timestep(int))
    # TODO # update_multiplier!(cache.λ, cache.Λp, int.params.t_λ.b)

    # copy solution to initial guess
    update!(int.iguess, cache.t, cache.q, cache.p, cache.v, cache.f)

    # take care of periodic solutions
    cut_periodic_solution!(cache, equation(int).periodicity)
end
