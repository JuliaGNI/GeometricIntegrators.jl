
struct TableauVPRKpLegendre{T} <: AbstractTableauPRK{T}
    @HeaderTableau

    q::CoefficientsRK{T}
    p::CoefficientsRK{T}

    R∞::Int

    ω::Matrix{T}
    d::Vector{T}

    function TableauVPRKpLegendre{T}(name, o, q, p, R∞, ω, d) where {T}
        @assert q.s == p.s == length(d)
        new(name, o, q.s, q, p, R∞, ω, d)
    end

    function TableauVPRKpLegendre{T}(name, o, q, p, R∞, ω) where {T}
        @assert q.s == p.s
        new(name, o, q.s, q, p, R∞, ω)
    end
end

function TableauVPRKpLegendre(name::Symbol, order::Int, q::CoefficientsRK{T}, p::CoefficientsRK{T}, R∞::Int, ω::Matrix{T}, d::Vector{T}) where {T}
    TableauVPRKpLegendre{T}(name, order, q, p, R∞, ω, d)
end

function TableauVPRKpLegendre(name::Symbol, order::Int, q::CoefficientsRK{T}, p::CoefficientsRK{T}, R∞::Int, ω::Matrix{T}) where {T}
    TableauVPRKpLegendre{T}(name, order, q, p, R∞, ω)
end

function TableauVPRKpLegendre(name::Symbol, order::Int, q::CoefficientsRK{T}, R∞::Int, ω::Matrix{T}, d::Vector{T}) where {T}
    TableauVPRKpLegendre{T}(name, order, q, get_symplectic_conjugate_coefficients(q), R∞, ω, d)
end

function TableauVPRKpLegendre(name::Symbol, order::Int, q::CoefficientsRK{T}, R∞::Int, ω::Matrix{T}) where {T}
    TableauVPRKpLegendre{T}(name, order, q, get_symplectic_conjugate_coefficients(q), R∞, ω)
end


"Parameters for right-hand side function of variational special partitioned additive Runge-Kutta methods."
mutable struct ParametersVPRKpLegendre{DT,TT,ΘT,FT,D,S} <: Parameters{DT,TT}
    f_θ::ΘT
    f_f::FT

    Δt::TT

    t_q::CoefficientsRK{TT}
    t_p::CoefficientsRK{TT}
    d_v::Vector{TT}

    t::TT

    q::Vector{DT}
    p::Vector{DT}


    function ParametersVPRKpLegendre{DT,TT,ΘT,FT,D,S}(f_θ, f_f, Δt, t_q, t_p, d_v, q, p) where {DT,TT,ΘT,FT,D,S}
        new(f_θ, f_f, Δt, t_q, t_p, d_v, 0, q, p)
    end
end


struct NonlinearFunctionCacheVPRKpLegendre{ST}

    q̅::Vector{ST}
    p̅::Vector{ST}
    ϕ::Vector{ST}
    μ::Vector{ST}

    v::Vector{ST}
    y::Vector{ST}
    z::Vector{ST}

    Q::Matrix{ST}
    P::Matrix{ST}
    V::Matrix{ST}
    F::Matrix{ST}
    Φ::Matrix{ST}

    function NonlinearFunctionCacheVPRKpLegendre{ST}(d,s) where {ST}
        # create solution vectors
        q̅ = zeros(ST,d)
        p̅ = zeros(ST,d)
        ϕ = zeros(ST,d)
        μ = zeros(ST,d)

        v = zeros(ST,d)
        y = zeros(ST,d)
        z = zeros(ST,d)

        # create internal stage vectors
        Q = zeros(ST,d,s)
        P = zeros(ST,d,s)
        V = zeros(ST,d,s)
        F = zeros(ST,d,s)
        Φ = zeros(ST,d,s)

        new(q̅, p̅, ϕ, μ, v, y, z, Q, P, V, F, Φ)
    end
end


"Compute stages of variational special partitioned additive Runge-Kutta methods."
@generated function function_stages!(y::Vector{ST}, b::Vector{ST}, params::ParametersVPRKpLegendre{DT,TT,ΘT,FT,D,S}) where {ST,DT,TT,ΘT,FT,D,S}
    cache = NonlinearFunctionCacheVPRKpLegendre{ST}(D,S)

    quote
        @assert length(y) == length(b)

        compute_stages!(y, $cache.Q, $cache.V, $cache.P, $cache.F, $cache.Φ, $cache.q̅, $cache.p̅, $cache.ϕ, $cache.μ, params)

        compute_rhs!(b, $cache.Q, $cache.V, $cache.P, $cache.F, $cache.Φ, $cache.q̅, $cache.p̅, $cache.ϕ, $cache.μ, params)

        # debug output
        # println()
        # for k in 1:D
        #     println(params.q[k], ",  ", params.p[k], ",  ", $cache.Q[k,1], ",  ", $cache.P[k,1], ",  ", $cache.q̅[k], ",  ", $cache.p̅[k])
        # end
        # println()
    end
end

@generated function compute_stages!(y::Vector{ST}, Q, V, P, F, Φ, q̅, p̅, ϕ, μ, params::ParametersVPRKpLegendre{DT,TT,ΘT,FT,D,S}) where {ST,DT,TT,ΘT,FT,D,S}
    # create temporary vectors
    Qt = zeros(ST,D)
    Pt = zeros(ST,D)
    Vt = zeros(ST,D)
    Ft = zeros(ST,D)

    quote
        local offset::Int
        local tqᵢ::TT

        for i in 1:S
            # copy y to Q and V
            for k in 1:D
                offset = 3*(D*(i-1)+k-1)
                Q[k,i] = y[offset+1]
                V[k,i] = y[offset+2]
                P[k,i] = y[offset+3]
            end

            # compute P=p(t,Q) and F=f(t,Q,V)
            tqᵢ = params.t + params.Δt * params.t_q.c[i]
            simd_copy_xy_first!($Qt, Q, i)
            simd_copy_xy_first!($Vt, V, i)
            params.f_θ(tqᵢ, $Qt, $Vt, $Pt)
            params.f_f(tqᵢ, $Qt, $Vt, $Ft)
            simd_copy_yx_first!($Pt, Φ, i)
            simd_copy_yx_first!($Ft, F, i)
        end

        # copy y to q̅ and p̅
        for k in 1:D
            q̅[k] = y[(3*S+0)*D+k]
            p̅[k] = y[(3*S+1)*D+k]
        end

        # compute p̅=p(t,q̅)
        $Vt .= 0
        params.f_θ(params.t + params.Δt, q̅, $Vt, ϕ)

        # for Lobatto-type methods, copy y to μ
        # if length(params.d_v) > 0
        #     offset = (3*S+2)*D
        #     for k in 1:D
        #         μ[k] = y[offset+k]
        #     end
        # end
    end
end


function compute_rhs!(b::Vector{ST}, Q::Matrix{ST}, V::Matrix{ST}, P::Matrix{ST}, F::Matrix{ST}, Φ::Matrix{ST}, q̅::Vector{ST}, p̅::Vector{ST}, ϕ::Vector{ST}, μ::Vector{ST}, params::ParametersVPRKpLegendre{DT,TT,ΘT,FT,D,S}) where {ST,DT,TT,ΘT,FT,D,S}
    local ω::ST
    local z::ST

    # compute b = - (Q-q-AV)
    for i in 1:S
        offset = 0*S*D + D*(i-1)
        for k in 1:D
            z = 0
            for j in 1:S
                z += params.t_q.a[i,j] * V[k,j]
            end
            b[offset+k] = (params.q[k] - Q[k,i]) + params.Δt * z
        end
    end

    # compute b = - (P-p-AF)
    for i in 1:S
        offset = 1*S*D + D*(i-1)
        for k in 1:D
            z = 0
            for j in 1:S
                z += params.t_p.a[i,j] * F[k,j]
            end
            b[offset+k] = (params.p[k] - P[k,i]) + params.Δt * z
        end
    end

    # compute b = - (P-Φ)
    for i in 1:S-1
        offset = 2*S*D + D*(i-1)
        for k in 1:D
            z = 0
            for j in 1:S
                ω = params.t_p.b[j] * params.t_p.c[j]^(i-1)
                # ω = params.t_q.a[i+1,j]
                z += ω * (Φ[k,j] - P[k,j])
            end
            b[offset+k] = z
        end
    end

    # compute b = - (q̅-q-AV)
    offset = (3*S-1)*D
    for k in 1:D
        z = 0
        for j in 1:S
            z += params.t_q.b[j] * V[k,j]
        end
        b[offset+k] = (params.q[k] - q̅[k]) + params.Δt * z
    end

    # compute b = - (p̅-p-AF)
    offset = (3*S+0)*D
    for k in 1:D
        z = 0
        for j in 1:S
            z += params.t_p.b[j] * F[k,j]
        end
        b[offset+k] = (params.p[k] - p̅[k]) + params.Δt * z
    end

    # compute b = - (p̅-ϕ)
    offset = (3*S+1)*D
    for k in 1:D
        b[offset+k] = ϕ[k] - p̅[k]
    end

    # if length(params.d_v) > 0
    #     offset = (3*S+2)*D
    #     for k in 1:D
    #         b[offset+k] = 0
    #         for i in 1:S
    #             b[offset+k] -= V[k,i] * params.d_v[i]
    #         end
    #     end
    # end

    if length(params.d_v) > 0
        sl     = div(S+1, 2)
        offset = D*(S+sl-1)

        # compute μ
        z = params.t_p.b[sl] / params.d_v[sl]
        for k in 1:D
            μ[k] = z * b[offset+k]
        end

        # replace equation for Pₗ with constraint on V
        for k in 1:D
            z = 0
            for i in 1:S
                z += V[k,i] * params.d_v[i]
            end
            b[offset+k] = z
        end

        # modify P₁, ..., Pₛ except for Pₗ
        for i in 1:S
            if i ≠ sl
                offset = D*(S+i-1)
                z = params.d_v[i] / params.t_p.b[i]
                for k in 1:D
                    b[offset+k] -= z * μ[k]
                end
            end
        end
    end
end



"Variational special partitioned additive Runge-Kutta integrator."
mutable struct IntegratorVPRKpLegendre{DT, TT, ΘT, FT, GT, VT, SPT, ST, IT} <: Integrator{DT, TT}
    equation::IODE{DT,TT,ΘT,FT,GT,VT}
    tableau::TableauVPRK{TT}
    Δt::TT

    params::SPT
    solver::ST
    iguess::InitialGuessPODE{DT, TT, VT, FT, IT}

    q::Vector{DT}
    p::Vector{DT}

    cache::NonlinearFunctionCacheVPRKpLegendre{DT}
end

function IntegratorVPRKpLegendre(equation::IODE{DT,TT,ΘT,FT,GT,VT}, tableau::TableauVPRK{TT}, Δt::TT) where {DT,TT,ΘT,FT,GT,VT}
    D = equation.d
    S = tableau.s

    N = (3*S+2)*D

    if @isdefined(tableau, d)
        d_v = tableau.d
    else
        d_v = DT[]
    end

    # create solution vectors
    q = zeros(DT,D)
    p = zeros(DT,D)

    # create solution vector for internal stages / nonlinear solver
    z = zeros(DT,N)

    # create cache for internal stage vectors and update vectors
    cache = NonlinearFunctionCacheVPRKpLegendre{DT}(D,S)

    # create params
    params = ParametersVPRKpLegendre{DT,TT,ΘT,FT,D,S}(
                                                equation.α, equation.f,
                                                Δt, tableau.q, tableau.p, d_v, q, p)

    # create rhs function for nonlinear solver
    function_stages = (x,b) -> function_stages!(x, b, params)

    # create solver
    solver = get_config(:nls_solver)(z, function_stages)

    # create initial guess
    iguess = InitialGuessPODE(get_config(:ig_interpolation), equation, Δt)

    # create integrator
    IntegratorVPRKpLegendre{DT, TT, ΘT, FT, GT, VT, typeof(params), typeof(solver), typeof(iguess.int)}(
                equation, tableau, Δt, params, solver, iguess, q, p, cache)
end


function initialize!(int::IntegratorVPRKpLegendre, sol::SolutionPDAE, m::Int)
    @assert m ≥ 1
    @assert m ≤ sol.ni

    # copy initial conditions from solution
    get_initial_conditions!(sol, int.q, int.p, m)

    # initialise initial guess
    initialize!(int.iguess, m, sol.t[0], int.q, int.p)
end


function update_solution!(int::IntegratorVPRKpLegendre{DT,TT}, cache::NonlinearFunctionCacheVPRKpLegendre{DT}) where {DT,TT}
    update_solution!(int.q, cache.V, int.tableau.q.b, int.tableau.q.b̂, int.Δt)
    update_solution!(int.p, cache.F, int.tableau.p.b, int.tableau.p.b̂, int.Δt)
end


"Integrate DAE with variational special partitioned additive Runge-Kutta integrator."
function integrate_step!(int::IntegratorVPRKpLegendre{DT,TT,ΘT,FT,VT}, sol::SolutionPDAE{DT,TT,N}, m::Int, n::Int) where {DT,TT,ΘT,FT,VT,N}
    local offset::Int

    # set time for nonlinear solver
    int.params.t = sol.t[0] + (n-1)*int.Δt

    # compute initial guess
    for i in 1:int.tableau.s
        evaluate!(int.iguess, m, int.cache.y, int.cache.z, int.cache.v, int.tableau.q.c[i], int.tableau.p.c[i])
        for k in 1:int.equation.d
            offset = 3*(int.equation.d*(i-1)+k-1)
            int.solver.x[offset+1] = int.cache.y[k]
            int.solver.x[offset+2] = int.cache.v[k]
            int.solver.x[offset+3] = int.cache.z[k]
        end
    end

    evaluate!(int.iguess, m, int.cache.y, int.cache.z, int.cache.v, one(TT), one(TT))
    for k in 1:int.equation.d
        int.solver.x[(3*int.tableau.s+0)*int.equation.d+k] = int.cache.y[k]
        int.solver.x[(3*int.tableau.s+1)*int.equation.d+k] = int.cache.z[k]
    end

    # if @isdefined(int.tableau, d)
    #     offset = (3*int.tableau.s+2)*int.equation.d
    #     for k in 1:int.equation.d
    #         int.solver.x[offset+k] = 0
    #     end
    # end

    # debug output
    # println()
    # println(int.solver.x)

    # call nonlinear solver
    solve!(int.solver)

    # print solver status
    print_solver_status(int.solver.status, int.solver.params, n)

    # check if solution contains NaNs or error bounds are violated
    check_solver_status(int.solver.status, int.solver.params, n)

    # compute final update
    compute_stages!(int.solver.x, int.cache.Q, int.cache.V, int.cache.P, int.cache.F, int.cache.Φ, int.cache.q̅, int.cache.p̅, int.cache.ϕ, int.cache.μ, int.params)
    update_solution!(int, int.cache)

    # debug output
    # println(int.solver.x)
    # println()

    # copy solution to initial guess for next time step
    update!(int.iguess, m, sol.t[0] + n*int.Δt, int.q, int.p)

    # take care of periodic solutions
    # cut_periodic_solution!(int.q, int.equation.periodicity)

    # copy to solution
    copy_solution!(sol, int.q, int.p, n, m)
end
