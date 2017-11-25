"""
`ParametersDGVI`: Parameters for right-hand side function of discontinuous Galerkin variational Integrator.

### Parameters

* `Θ`: function of the noncanonical one-form (∂L/∂v)
* `f`: function of the force (∂L/∂q)
* `Δt`: time step
* `b`: quadrature weights
* `c`: quadrature nodes
* `m`: mass matrix
* `a`: derivative matrix
* `r₀`: reconstruction coefficients, left-hand side
* `r₁`: reconstruction coefficients, right-hand side
* `β`: weights of the quadrature rule for the flux
* `γ`: nodes of the quadrature rule for the flux
* `μ₀`: mass vector for the lhs flux
* `μ₁`: mass vector for the rhs flux
* `α₀`: derivative vector for the lhs flux
* `α₁`: derivative vector for the rhs flux
* `t`: current time
* `q`: current solution of q
* `p`: current solution of p
* `D`: dimension of the system
* `S`: number of basis nodes
* `R`: number of quadrature nodes
* `P`: number of quadrature nodes for the flux
"""
mutable struct ParametersDGVI{DT,TT,ΘT,FT,GT,D,S,QR,FR} <: Parameters{DT,TT}
    Θ::ΘT
    f::FT
    g::GT

    Δt::TT

    b::Vector{TT}
    c::Vector{TT}
    m::Matrix{TT}
    a::Matrix{TT}
    r₀::Vector{TT}
    r₁::Vector{TT}

    β::Vector{TT}
    γ::Vector{TT}
    μ₀::Vector{TT}
    μ₁::Vector{TT}
    α₀::Vector{TT}
    α₁::Vector{TT}

    t::TT

    q::Vector{DT}
    p::Vector{DT}
end

function ParametersDGVI(Θ::ΘT, f::FT, g::GT, Δt::TT,
                b::Vector{TT}, c::Vector{TT}, m::Matrix{TT}, a::Matrix{TT},
                r₀::Vector{TT}, r₁::Vector{TT}, β::Vector{TT}, γ::Vector{TT},
                μ₀::Vector{TT}, μ₁::Vector{TT}, α₀::Vector{TT}, α₁::Vector{TT},
                q::Vector{DT}, p::Vector{DT}) where {DT,TT,ΘT,FT,GT}

    @assert length(q)  == length(p)
    @assert length(b)  == length(c)
    @assert length(β)  == length(γ)
    @assert length(r₀) == length(r₁)

    D  = length(q)
    S  = length(r₀)
    QR = length(c)
    FR = length(γ)

    println()
    println("  Discontinuous Galerkin Variational Integrator")
    println("  =============================================")
    println()
    println("    c = ", c)
    println("    b = ", b)
    println("    m = ", m)
    println("    a = ", a)
    println("    r₀= ", r₀)
    println("    r₁= ", r₁)
    println("    β = ", β)
    println("    γ = ", γ)
    println("    μ₀= ", μ₀)
    println("    μ₁= ", μ₁)
    println("    α₀= ", α₀)
    println("    α₁= ", α₁)
    println()

    ParametersDGVI{DT,TT,ΘT,FT,GT,D,S,QR,FR}(
                Θ, f, g, Δt, b, c, m, a, r₀, r₁, β, γ, μ₀, μ₁, α₀, α₁, 0, q, p)
end


function ParametersDGVI(Θ::ΘT, f::FT, g::GT, Δt::TT,
                basis::Basis{TT}, quadrature::Quadrature{TT}, flux::Flux{TT},
                q::Vector{DT}, p::Vector{DT}) where {DT,TT,ΘT,FT,GT}

    # compute coefficients
    m = zeros(TT, nnodes(quadrature), nbasis(basis))
    a = zeros(TT, nnodes(quadrature), nbasis(basis))
    r₀= zeros(TT, nbasis(basis))
    r₁= zeros(TT, nbasis(basis))

    for i in 1:nbasis(basis)
        for j in 1:nnodes(quadrature)
            m[j,i] = evaluate(basis, i, nodes(quadrature)[j])
            a[j,i] = derivative(basis, i, nodes(quadrature)[j])
        end
        r₀[i] = evaluate(basis, i, zero(TT))
        r₁[i] = evaluate(basis, i, one(TT))
        # r₀[i] = evaluate(basis, i, one(TT)/2)
        # r₁[i] = evaluate(basis, i, one(TT)/2)
    end

    μ₀ = zeros(TT, nnodes(flux.quadrature))
    μ₁ = zeros(TT, nnodes(flux.quadrature))
    α₀ = zeros(TT, nnodes(flux.quadrature))
    α₁ = zeros(TT, nnodes(flux.quadrature))

    for i in 1:nnodes(flux.quadrature)
        μ₀[i] = evaluate_l(flux.path, nodes(flux.quadrature)[i])
        μ₁[i] = evaluate_r(flux.path, nodes(flux.quadrature)[i])
        α₀[i] = derivative_l(flux.path, nodes(flux.quadrature)[i])
        α₁[i] = derivative_r(flux.path, nodes(flux.quadrature)[i])
    end

    ParametersDGVI(
                    Θ, f, g, Δt, weights(quadrature), nodes(quadrature), m, a, r₀, r₁,
                    weights(flux.quadrature), nodes(flux.quadrature), μ₀, μ₁, α₀, α₁, q, p)
end


struct NonlinearFunctionCacheDGVI{ST}
    X::Matrix{ST}
    Q::Matrix{ST}
    V::Matrix{ST}
    P::Matrix{ST}
    F::Matrix{ST}
    q̅::Vector{ST}
    p̅::Vector{ST}
    λ::Matrix{ST}
    λ̅::Matrix{ST}
    ϕ::Matrix{ST}
    ϕ̅::Matrix{ST}

    function NonlinearFunctionCacheDGVI{ST}(D,S,QR,FR) where {ST}
        # create internal stage vectors
        X = zeros(ST,D,S)
        Q = zeros(ST,D,QR)
        V = zeros(ST,D,QR)
        P = zeros(ST,D,QR)
        F = zeros(ST,D,QR)
        q̅ = zeros(ST,D)
        p̅ = zeros(ST,D)
        λ = zeros(ST,D,FR)
        λ̅ = zeros(ST,D,FR)
        ϕ = zeros(ST,D,FR)
        ϕ̅ = zeros(ST,D,FR)

        new(X, Q, V, P, F, q̅, p̅, λ, λ̅, ϕ, ϕ̅)
    end
end


"Compute stages of variational partitioned Runge-Kutta methods."
@generated function function_stages!(x::Vector{ST}, b::Vector{ST}, params::ParametersDGVI{DT,TT,ΘT,FT,GT,D,S,QR,FR}) where {ST,DT,TT,ΘT,FT,GT,D,S,QR,FR}
    cache = NonlinearFunctionCacheDGVI{ST}(D, S, QR, FR)

    quote
        @assert length(x) == length(b)

        compute_stages!(x, $cache.X, $cache.Q, $cache.V, $cache.P, $cache.F,
                           $cache.q̅, $cache.p̅, $cache.λ, $cache.λ̅,
                           $cache.ϕ, $cache.ϕ̅, params)

        compute_rhs!(b, $cache.X, $cache.P, $cache.F,
                        $cache.λ, $cache.λ̅, $cache.ϕ, $cache.ϕ̅, params)
    end
end


function compute_stages!(x, X, Q, V, P, F, q̅, p̅, λ, λ̅, ϕ, ϕ̅, params::ParametersDGVI{DT,TT,ΘT,FT,GT,D,S,QR,FR}) where {DT,TT,ΘT,FT,GT,D,S,QR,FR}

    # copy x to X
    for i in 1:S
        for k in 1:D
            X[k,i] = x[D*(i-1)+k]
        end
    end

    # copy x to q̅=qₙ+₁
    for k in 1:D
        q̅[k] = x[D*S+k]
    end

    # compute Q
    compute_stages_q!(X, Q, p̅, params)

    # compute V
    compute_stages_v!(X, V, params)

    # compute P and F
    compute_stages_p!(Q, V, P, F, params)

    # compute jump
    compute_stages_λ!(X, q̅, λ, λ̅, ϕ, ϕ̅, params)

    # debug output
    # println()
    # for k in 1:D
    #     println(params.p[k], ",  ", params.q[k], ",  ", ϕ[k], ",  ", p̅[k], ",  ", q̅[k], ",  ", ϕ̅[k])
    # end
    # println()
end


function compute_stages_q!(X::Matrix{ST}, Q::Matrix{ST}, p̅::Vector{ST},
            params::ParametersDGVI{DT,TT,AT,FT,GT,D,S,R}) where {ST,DT,TT,AT,FT,GT,D,S,R}

    @assert D == size(Q,1) == size(X,1)
    @assert R == size(Q,2)
    # @assert S == size(X,2)

    local y::ST

    # compute Q
    for i in 1:size(Q,2)
        for k in 1:size(Q,1)
            y = 0
            for j in 1:size(X,2)
                y += params.m[i,j] * X[k,j]
            end
            Q[k,i] = y
        end
    end

    # compute p̅=q_n+1^-
    for k in 1:size(X,1)
        y = 0
        for i in 1:size(X,2)
            y += params.r₁[i] * X[k,i]
        end
        p̅[k] = y
    end
end


function compute_stages_v!(X::Matrix{ST}, V::Matrix{ST}, params::ParametersDGVI{DT,TT,AT,FT,GT,D,S,R}) where {ST,DT,TT,AT,FT,GT,D,S,R}
    @assert D == size(V,1) == size(X,1)
    @assert R == size(V,2)
    @assert S == size(X,2)

    local y::ST

    # compute V
    for i in 1:R
        for k in 1:D
            y = 0
            for j in 1:S
                y += params.a[i,j] * X[k,j]
            end
            V[k,i] = y / params.Δt
        end
    end
end


@generated function compute_stages_p!(Q::Matrix{ST}, V::Matrix{ST}, P::Matrix{ST}, F::Matrix{ST},
            params::ParametersDGVI{DT,TT,AT,FT,GT,D,S,R}) where {ST,DT,TT,AT,FT,GT,D,S,R}

    # create temporary vectors
    tQ = zeros(ST,D)
    tV = zeros(ST,D)
    tP = zeros(ST,D)
    tF = zeros(ST,D)

    quote
        @assert D == size(Q,1) == size(V,1) == size(P,1) == size(F,1)
        @assert R == size(Q,2) == size(V,2) == size(P,2) == size(F,2)

        local tᵢ::TT

        # compute P=α(Q) and F=f(Q)
        for i in 1:R
            tᵢ = params.t + params.Δt * params.c[i]
            simd_copy_xy_first!($tQ, Q, i)
            simd_copy_xy_first!($tV, V, i)
            params.Θ(tᵢ, $tQ, $tV, $tP)
            params.f(tᵢ, $tQ, $tV, $tF)
            simd_copy_yx_first!($tP, P, i)
            simd_copy_yx_first!($tF, F, i)
        end
    end
end


function compute_stages_λ!(X::Matrix{ST}, q̅::Vector{ST}, λ::Matrix{ST}, λ̅::Matrix{ST}, ϕ::Matrix{ST}, ϕ̅::Matrix{ST},
                params::ParametersDGVI{DT,TT,AT,FT,GT,D,S,QR,FR}) where {ST,DT,TT,AT,FT,GT,D,S,QR,FR}

    @assert D == size(X,1)
    @assert S == size(X,2)

    local q = params.q

    local y₀::ST # qₙ^+
    local y₁::ST # qₙ+₁^-

    # compute λ and λ̅
    for k in 1:D
        y₀ = 0
        y₁ = 0
        for j in 1:S
            y₀ += params.r₀[j] * X[k,j]
            y₁ += params.r₁[j] * X[k,j]
        end

        for i in 1:FR
            λ[k,i] = params.α₀[i] * q[k] + params.α₁[i] * y₀
            λ̅[k,i] = params.α₀[i] * y₁ + params.α₁[i] * q̅[k]
        end

        for i in 1:FR
            ϕ[k,i] = 2 * params.μ₁[i] * q[k] + (params.μ₀[i] - params.μ₁[i]) * y₀
            ϕ̅[k,i] = 2 * params.μ₀[i] * q̅[k] + (params.μ₁[i] - params.μ₀[i]) * y₁
            # ϕ[k,i] = y₀
            # ϕ̅[k,i] = y₁
            # ϕ[k,i] = q[k]
            # ϕ̅[k,i] = q̅[k]
        end
    end
end


@generated function compute_rhs!(b::Vector{ST}, X::Matrix{ST}, P::Matrix{ST}, F::Matrix{ST},
                λ::Matrix{ST}, λ̅::Matrix{ST}, ϕ::Matrix{ST}, ϕ̅::Matrix{ST},
                params::ParametersDGVI{DT,TT,AT,FT,GT,D,S,QR,FR}) where {ST,DT,TT,AT,FT,GT,D,S,QR,FR}

    local θ::Matrix{ST} = zeros(ST,D,FR)
    local Θ̅::Matrix{ST} = zeros(ST,D,FR)
    local g::Matrix{ST} = zeros(ST,D,FR)
    local g̅::Matrix{ST} = zeros(ST,D,FR)

    local tϕ = zeros(ST,D)
    local tλ = zeros(ST,D)
    local tθ = zeros(ST,D)
    local tg = zeros(ST,D)

    quote
        local y::ST
        local z::ST

        local t₀ = params.t
        local t₁ = params.t + params.Δt

        for j in 1:FR
            simd_copy_xy_first!($tϕ, ϕ, j)
            simd_copy_xy_first!($tλ, λ, j)
            params.Θ(t₀, $tϕ, $tλ, $tθ)
            params.g(t₀, $tϕ, $tλ, $tg)
            simd_copy_yx_first!($tθ, $θ, j)
            simd_copy_yx_first!($tg, $g, j)

            simd_copy_xy_first!($tϕ, ϕ̅, j)
            simd_copy_xy_first!($tλ, λ̅, j)
            params.Θ(t₁, $tϕ, $tλ, $tθ)
            params.g(t₁, $tϕ, $tλ, $tg)
            simd_copy_yx_first!($tθ, $Θ̅, j)
            simd_copy_yx_first!($tg, $g̅, j)
        end

        # compute b = - [(P-AF)]
        for i in 1:S
            for k in 1:D
                z = 0
                for j in 1:QR
                    z += params.b[j] * params.m[j,i] * F[k,j] * params.Δt
                    z += params.b[j] * params.a[j,i] * P[k,j]
                end
                for j in 1:FR
                    z += 1 * params.β[j] * params.r₀[i] * params.α₁[j] * $θ[k,j]
                    z += 1 * params.β[j] * params.r₁[i] * params.α₀[j] * $Θ̅[k,j]
                    z += 2 * params.β[j] * params.r₀[i] * params.μ₀[j] * $g[k,j]
                    z += 2 * params.β[j] * params.r₁[i] * params.μ₁[j] * $g̅[k,j]

                    # z += params.β[j] * params.r₀[i] * params.α₁[j] * $θ[k,j]
                    # z += params.β[j] * params.r₁[i] * params.α₀[j] * $Θ̅[k,j]
                    # z += params.β[j] * params.r₀[i] * params.μ₀[j] * $g[k,j]
                    # z += params.β[j] * params.r₁[i] * params.μ₁[j] * $g̅[k,j]

                    # z += params.β[j] * params.r₀[i] * $θ[k,j]
                    # z -= params.β[j] * params.r₁[i] * $Θ̅[k,j]
                    # z += params.β[j] * params.r₀[i] * $g[k,j]
                    # z -= params.β[j] * params.r₁[i] * $g̅[k,j]
                end

                b[D*(i-1)+k] = z
            end
        end

        # compute b = - [(2qₙ - qₙ^- - qₙ^+)]
        for k in 1:size(X,1)
            y = 0
            for j in 1:size(X,2)
                y += params.r₀[j] * X[k,j]
            end
            b[D*S+k] = 2 * params.q[k] - params.p[k] - y
        end
    end
end


"Discontinuous Galerkin Variational Integrator."
struct IntegratorDGVI{DT,TT,ΘT,FT,GT,VT,FPT,ST,IT,BT<:Basis,FLT<:Flux,D,S,R} <: Integrator{DT,TT}
    equation::IODE{DT,TT,ΘT,FT,GT,VT}

    basis::BT
    quadrature::Quadrature{TT,R}
    flux::FLT

    Δt::TT

    params::FPT
    solver::ST
    iguess::InitialGuessPODE{DT,TT,VT,FT,IT}

    q::Vector{DT}
    p::Vector{DT}

    cache::NonlinearFunctionCacheDGVI{DT}
end

function IntegratorDGVI(equation::IODE{DT,TT,ΘT,FT,GT,VT}, basis::Basis{TT,P},
                quadrature::Quadrature{TT,R}, flux::Flux{TT,PT,QN}, Δt::TT;
                interpolation=HermiteInterpolation{DT}) where {DT,TT,ΘT,FT,GT,VT,P,R,PT,QN}

    D = equation.d
    S = nbasis(basis)

    N = D*(S+1)

    # create solution vector for nonlinear solver
    x = zeros(DT,N)

    # create solution vectors
    q = zeros(DT,D)
    p = zeros(DT,D)

    # create cache for internal stage vectors and update vectors
    cache = NonlinearFunctionCacheDGVI{DT}(D,S,R,QN)

    # create params
    params = ParametersDGVI(equation.α, equation.f, equation.g,
                Δt, basis, quadrature, flux, q, p)

    # create rhs function for nonlinear solver
    function_stages = (x,b) -> function_stages!(x, b, params)

    # create nonlinear solver
    solver = get_config(:nls_solver)(x, function_stages)

    # create initial guess
    iguess = InitialGuessPODE(interpolation, equation, Δt)

    # create integrator
    IntegratorDGVI{DT, TT, ΘT, FT, GT, VT, typeof(params), typeof(solver),
                typeof(iguess.int), typeof(basis), typeof(flux), D, S, R}(
                equation, basis, quadrature, flux, Δt, params, solver, iguess,
                q, p, cache)
end



function initialize!(int::IntegratorDGVI, sol::Union{SolutionPODE, SolutionPDAE}, m::Int)
    @assert m ≥ 1
    @assert m ≤ sol.ni

    # copy initial conditions from solution
    get_initial_conditions!(sol, int.q, int.p, m)

    # initialise initial guess
    initialize!(int.iguess, m, sol.t[0], int.q, int.p)
end


function update_solution!(int::IntegratorDGVI{DT,TT}, cache::NonlinearFunctionCacheDGVI{DT}) where {DT,TT}
    for k in eachindex(int.q, cache.q̅)
        int.q[k] = cache.q̅[k]
    end
    for k in eachindex(int.p, cache.p̅)
        int.p[k] = cache.p̅[k]
    end
end


@generated function initial_guess!(int::IntegratorDGVI{DT,TT,ΘT,FT,GT,VT,FPT,ST,IT,BT,FLT,D,S,R}, m::Int) where {DT,TT,ΘT,FT,GT,VT,FPT,ST,IT,BT,FLT,D,S,R}
    v = zeros(DT,D)
    y = zeros(DT,D)
    z = zeros(DT,D)

    quote
        # compute initial guess
        if nnodes(int.basis) > 0
            for i in 1:S
                evaluate!(int.iguess, m, $y, $z, $v, nodes(int.basis)[i], nodes(int.basis)[i])
                for k in 1:D
                    int.solver.x[D*(i-1)+k] = $y[k]
                end
            end
        else
            for i in 1:S
                for k in 1:D
                    int.solver.x[D*(i-1)+k] = 0
                end
            end
        end

        evaluate!(int.iguess, m, $y, $z, $v, one(TT), one(TT))
        for k in 1:D
            int.solver.x[D*S+k] = $z[k]
        end
    end
end


"Integrate ODE with variational partitioned Runge-Kutta integrator."
function integrate_step!(int::IntegratorDGVI{DT,TT}, sol::Union{SolutionPODE{DT,TT}, SolutionPDAE{DT,TT}}, m::Int, n::Int) where {DT,TT}
    # set time for nonlinear solver
    int.params.t = sol.t[0] + (n-1)*int.Δt

    # compute initial guess
    initial_guess!(int, m)

    # call nonlinear solver
    solve!(int.solver)

    # print solver status
    print_solver_status(int.solver.status, int.solver.params, n)

    # check if solution contains NaNs or error bounds are violated
    check_solver_status(int.solver.status, int.solver.params, n)

    # compute final update
    compute_stages!(int.solver.x, int.cache.X, int.cache.Q, int.cache.V, int.cache.P, int.cache.F,
                    int.cache.q̅, int.cache.p̅, int.cache.λ, int.cache.λ̅, int.cache.ϕ, int.cache.ϕ̅, int.params)

    update_solution!(int, int.cache)

    # debug output
    # println(int.q)
    # println(int.p)
    # println()

    # copy solution to initial guess for next time step
    update!(int.iguess, m, sol.t[0] + n*int.Δt, int.q, int.p)

    # take care of periodic solutions
    cut_periodic_solution!(int.q, int.equation.periodicity)
    cut_periodic_solution!(int.p, int.equation.periodicity)

    # copy to solution
    copy_solution!(sol, int.q, int.p, n, m)
end
