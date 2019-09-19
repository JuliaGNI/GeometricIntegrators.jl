
TableauHPARK = TableauVPARK


"Parameters for right-hand side function of Hamiltonian partitioned additive Runge-Kutta methods."
mutable struct ParametersHPARK{DT,TT,D,S,R,VT,FT,UT,GT,ϕT} <: Parameters{DT,TT}
    f_v::VT
    f_f::FT
    f_u::UT
    f_g::GT
    f_ϕ::ϕT

    Δt::TT

    t_q::CoefficientsARK{TT}
    t_p::CoefficientsARK{TT}
    t_q̃::CoefficientsPRK{TT}
    t_p̃::CoefficientsPRK{TT}
    t_λ::CoefficientsMRK{TT}
    d_v::Vector{TT}

    t::TT
    q::Vector{DT}
    p::Vector{DT}
    λ::Vector{DT}

    function ParametersHPARK{DT,TT,D,S,R,VT,FT,UT,GT,ϕT}(f_v, f_f, f_u, f_g, f_ϕ, Δt, t_q, t_p, t_q̃, t_p̃, t_λ, d_v) where {DT,TT,D,S,R,VT,FT,UT,GT,ϕT}
        # create solution vectors
        q = zeros(DT,D)
        p = zeros(DT,D)
        λ = zeros(DT,D)

        new(f_v, f_f, f_u, f_g, f_ϕ, Δt,
            t_q, t_p, t_q̃, t_p̃, t_λ, d_v,
            zero(TT), q, p, λ)
    end
end


function compute_stages!(x::Vector{ST}, cache::IntegratorCacheVPARK{ST,TT,D,S,R},
                                        params::ParametersHPARK{DT,TT,D,S,R}) where {ST,DT,TT,D,S,R}
    local tqᵢ::TT
    local tpᵢ::TT
    local tλᵢ::TT

    for i in 1:S
        for k in 1:D
            # copy x to Y, Z
            cache.Yi[i][k] = x[2*(D*(i-1)+k-1)+1]
            cache.Zi[i][k] = x[2*(D*(i-1)+k-1)+2]

            # compute Q and P
            cache.Qi[i][k] = params.q[k] + params.Δt * cache.Yi[i][k]
            cache.Pi[i][k] = params.p[k] + params.Δt * cache.Zi[i][k]
        end

        # compute f(X)
        tqᵢ = params.t + params.Δt * params.t_q.c[i]
        tpᵢ = params.t + params.Δt * params.t_p.c[i]
        params.f_v(tqᵢ, cache.Qi[i], cache.Pi[i], cache.Vi[i])
        params.f_f(tpᵢ, cache.Qi[i], cache.Pi[i], cache.Fi[i])
    end

    for i in 1:R
        for k in 1:D
            # copy y to Y, Z and Λ
            cache.Yp[i][k] = x[2*D*S+3*(D*(i-1)+k-1)+1]
            cache.Zp[i][k] = x[2*D*S+3*(D*(i-1)+k-1)+2]
            cache.Λp[i][k] = x[2*D*S+3*(D*(i-1)+k-1)+3]

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
            cache.μ[k] = x[2*D*S+3*D*R+k]
        end
    end
end


"Compute stages of variational partitioned additive Runge-Kutta methods."
@generated function function_stages!(y::Vector{ST}, b::Vector{ST}, params::ParametersHPARK{DT,TT,D,S,R}) where {ST,DT,TT,D,S,R}
    cache = IntegratorCacheVPARK{ST,TT,D,S,R}()

    quote
        compute_stages!(y, $cache, params)

        # compute b = - [(Y-AV-AU), (Z-AF-AG), Φ]
        for i in 1:S
            for k in 1:D
                b[2*(D*(i-1)+k-1)+1] = - $cache.Yi[i][k]
                b[2*(D*(i-1)+k-1)+2] = - $cache.Zi[i][k]
                for j in 1:S
                    b[2*(D*(i-1)+k-1)+1] += params.t_q.a[i,j] * $cache.Vi[j][k]
                    b[2*(D*(i-1)+k-1)+2] += params.t_p.a[i,j] * $cache.Fi[j][k]
                end
                for j in 1:R
                    b[2*(D*(i-1)+k-1)+1] += params.t_q.α[i,j] * $cache.Up[j][k]
                    b[2*(D*(i-1)+k-1)+2] += params.t_p.α[i,j] * $cache.Gp[j][k]
                end
            end
        end

        # compute b = - [(Y-AV-AU), (Z-AF-AG), Φ]
        for i in 1:R
            for k in 1:D
                b[2*D*S+3*(D*(i-1)+k-1)+1] = - $cache.Yp[i][k]
                b[2*D*S+3*(D*(i-1)+k-1)+2] = - $cache.Zp[i][k]
                b[2*D*S+3*(D*(i-1)+k-1)+3] = - $cache.Φp[i][k]
                for j in 1:S
                    b[2*D*S+3*(D*(i-1)+k-1)+1] += params.t_q̃.a[i,j] * $cache.Vi[j][k]
                    b[2*D*S+3*(D*(i-1)+k-1)+2] += params.t_p̃.a[i,j] * $cache.Fi[j][k]
                end
                for j in 1:R
                    b[2*D*S+3*(D*(i-1)+k-1)+1] += params.t_q̃.α[i,j] * $cache.Up[j][k]
                    b[2*D*S+3*(D*(i-1)+k-1)+2] += params.t_p̃.α[i,j] * $cache.Gp[j][k]
                end
            end
        end

        # compute b = - [Λ₁-λ]
        if params.t_λ.c[1] == 0
            for k in 1:D
                b[2*D*S+3*(k-1)+3] = - $cache.Λp[1][k] + params.λ[k]
            end
        end

        if length(params.d_v) > 0
            for i in 1:S
                for k in 1:D
                    b[2*(D*(i-1)+k-1)+3] -= $cache.μ[k] * params.d_v[i]
                end
            end

            for k in 1:D
                b[2*D*S+3*D*R+k] = 0
                for i in 1:S
                    b[2*D*S+3*D*R+k] -= $cache.Vi[i][k] * params.d_v[i]
                end
            end
        end
    end
end


@doc raw"""
Partitioned Additive Runge-Kutta integrator for Hamiltonian systems subject
to Dirac constraints.

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
V_{n,i} &= \hphantom{-} \frac{\partial H}{\partial p} (Q_{n,i}, P_{n,i}) , & i &= 1, ..., s , \\
F_{n,i} &=           -  \frac{\partial H}{\partial q} (Q_{n,i}, P_{n,i}) , & i &= 1, ..., s , \\
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
struct IntegratorHPARK{DT, TT, ET <: PDAE{DT,TT},
                               PT <: ParametersHPARK{DT,TT},
                               ST <: NonlinearSolver{DT},
                               IT <: InitialGuessPODE{DT,TT}} <: AbstractIntegratorSPARK{DT, TT}
    equation::ET
    tableau::TableauHPARK{TT}

    params::PT
    solver::ST
    iguess::IT
end

function IntegratorHPARK(equation::PDAE{DT,TT,FT,PT,UT,GT,ϕT,VT},
                         tableau::TableauHPARK{TT}, Δt::TT) where {DT,TT,FT,PT,UT,GT,ϕT,VT}
    D = equation.d
    S = tableau.s
    R = tableau.r

    N = 2*D*S + 3*D*R

    if isdefined(tableau, :d)
        N += D
        d_v = tableau.d
    else
        d_v = DT[]
    end

    # create params
    params = ParametersHPARK{DT,TT,D,S,R,FT,PT,UT,GT,ϕT}(
                                equation.v, equation.f, equation.u, equation.g, equation.ϕ, Δt,
                                tableau.q, tableau.p, tableau.q̃, tableau.p̃, tableau.λ, d_v)

    # create solver
    solver = create_nonlinear_solver(DT, N, params)

    # create initial guess
    iguess = InitialGuessPODE(get_config(:ig_interpolation), equation, Δt)

    # create integrator
    IntegratorHPARK{DT, TT, typeof(equation), typeof(params), typeof(solver), typeof(iguess)}(
                                        equation, tableau, params, solver, iguess)
end

equation(int::IntegratorHPARK) = int.equation
timestep(int::IntegratorHPARK) = int.params.Δt
tableau(int::IntegratorHPARK) = int.tableau
nstages(int::IntegratorHPARK) = int.tableau.s
pstages(int::IntegratorHPARK) = int.tableau.r


function create_integrator_cache(int::IntegratorHPARK{DT,TT}) where {DT,TT}
    IntegratorCacheVPARK{DT, TT, ndims(int), nstages(int), pstages(int)}()
end


function update_params!(params::ParametersHPARK, cache::IntegratorCacheVPARK)
    # set time for nonlinear solver and copy previous solution
    params.t  = cache.t
    params.q .= cache.q
    params.p .= cache.p
    params.λ .= cache.λ
end


function initialize!(int::IntegratorHPARK, cache::IntegratorCacheVPARK)
    cache.t̅ = cache.t - timestep(int)

    equation(int).v(cache.t, cache.q, cache.p, cache.v)
    equation(int).f(cache.t, cache.q, cache.p, cache.f)

    initialize!(int.iguess, cache.t, cache.q, cache.p, cache.v, cache.f,
                            cache.t̅, cache.q̅, cache.p̅, cache.v̅, cache.f̅)
end


function initial_guess!(int::IntegratorHPARK, cache::IntegratorCacheVPARK)
    for i in 1:nstages(int)
        evaluate!(int.iguess, cache.q, cache.p, cache.v, cache.f,
                              cache.q̅, cache.p̅, cache.v̅, cache.f̅,
                              cache.q̃, cache.p̃, cache.ṽ, cache.f̃,
                              tableau(int).q.c[i], tableau(int).p.c[i])

        for k in 1:ndims(int)
            int.solver.x[2*(ndims(int)*(i-1)+k-1)+1] = (cache.q̃[k] - cache.q[k])/timestep(int)
            int.solver.x[2*(ndims(int)*(i-1)+k-1)+2] = (cache.p̃[k] - cache.p[k])/timestep(int)
        end
    end

    for i in 1:pstages(int)
        evaluate!(int.iguess, cache.q, cache.p, cache.v, cache.f,
                              cache.q̅, cache.p̅, cache.v̅, cache.f̅,
                              cache.q̃, cache.p̃, cache.ṽ, cache.f̃,
                              tableau(int).q̃.c[i], tableau(int).p̃.c[i])

        for k in 1:ndims(int)
            int.solver.x[2*ndims(int)*nstages(int)+3*(ndims(int)*(i-1)+k-1)+1] = (cache.q̃[k] - cache.q[k])/timestep(int)
            int.solver.x[2*ndims(int)*nstages(int)+3*(ndims(int)*(i-1)+k-1)+2] = (cache.p̃[k] - cache.p[k])/timestep(int)
            int.solver.x[2*ndims(int)*nstages(int)+3*(ndims(int)*(i-1)+k-1)+3] = 0
        end
    end

    if int.params.t_λ.c[1] == 0
        for k in 1:ndims(int)
            int.solver.x[2*ndims(int)*nstages(int)+3*(k-1)+3] = cache.λ[k]
        end
    end

    if isdefined(tableau(int), :d)
        for k in 1:ndims(int)
            int.solver.x[2*ndims(int)*nstages(int)+3*ndims(int)*pstages(int)+k] = 0
        end
    end
end


"Integrate DAE with variational partitioned additive Runge-Kutta integrator."
function integrate_step!(int::IntegratorHPARK{DT,TT}, cache::IntegratorCacheVPARK{DT,TT}) where {DT,TT}
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
    update_multiplier!(cache.λ, cache.Λp, int.params.t_λ.b)

    # copy solution to initial guess
    update!(int.iguess, cache.t, cache.q, cache.p, cache.v, cache.f)

    # take care of periodic solutions
    cut_periodic_solution!(cache, equation(int).periodicity)
end
