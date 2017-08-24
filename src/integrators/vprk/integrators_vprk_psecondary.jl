
"Parameters for right-hand side function of variational partitioned Runge-Kutta methods."
mutable struct NonlinearFunctionParametersVPRKpSecondary{DT,TT,ΑT,FT,GT,VT,ΩT,HT,D,S} <: AbstractNonlinearFunctionParametersVPRK{DT,TT,ΑT,FT,D,S}
    α::ΑT
    f::FT
    g::GT
    v::VT
    ω::ΩT
    dH::HT

    Δt::TT

    o::Int
    t_q::CoefficientsRK{TT}
    t_p::CoefficientsRK{TT}
    d_v::Vector{TT}
    ω_λ::Matrix{TT}
    R::Vector{TT}

    t::TT

    q::Vector{DT}
    p::Vector{DT}

    function NonlinearFunctionParametersVPRKpSecondary{DT,TT,ΑT,FT,GT,VT,ΩT,HT,D,S}(α, f, g, v, ω, dH, Δt, o, t_q, t_p, d_v, ω_λ, R∞, q, p) where {DT,TT,ΑT,FT,GT,VT,ΩT,HT,D,S}
        R = convert(Vector{TT}, [1, R∞])
        new(α, f, g, v, ω, dH, Δt, o, t_q, t_p, d_v, ω_λ, R, 0, q, p)
    end
end


function compute_stages_vprk!(x, q̅, p̅, λ, Q, V, Λ, U, P, F, R, G, Φ, params)
    # copy x to V
    compute_stages_v_vprk!(x, V, params)

    # copy x to Λ
    compute_stages_λ_vprk!(x, Λ, params)

    # compute U, G and p̅
    compute_projection_vprk!(x, q̅, p̅, λ, V, U, G, params)

    # compute Q
    compute_stages_q_vprk!(Q, V, U, params)

    # compute R and Φ
    compute_internal_projection_vprk!(x, Q, V, Λ, R, Φ, params)

    # compute P and F
    compute_stages_p_vprk!(Q, V, P, F, params)
end


@generated function compute_projection_vprk!(x::Vector{ST}, q̅::Vector{ST}, p̅::Vector{ST}, λ::Vector{ST}, V::Matrix{ST}, U::Matrix{ST}, G::Matrix{ST}, params::NonlinearFunctionParametersVPRKpSecondary{DT,TT,ΑT,FT,GT,VT,ΩT,HT,D,S}) where {ST,DT,TT,ΑT,FT,GT,VT,ΩT,HT,D,S}
    # create temporary vectors
    tG = zeros(ST,D)

    compute_projection_vprk = quote
        local t₀::TT = params.t
        local t₁::TT = params.t + params.Δt
        local y::ST

        # compute q̅
        for k in 1:D
            y = params.R[1] * U[k,1] + params.R[2] * U[k,2]
            for j in 1:S
                y += params.t_q.b[j] * V[k,j]
            end
            q̅[k] = params.q[k] + params.Δt * y
        end

        # compute U=λ at tₙ
        params.v(t₀, params.q, params.p, λ)
        simd_copy_yx_first!(λ, U, 1)

        # compute G=g(q,λ) at tₙ
        params.g(t₀, params.q, λ, $tG)
        simd_copy_yx_first!($tG, G, 1)

        # compute U=λ at tₙ+₁
        params.v(t₁, q̅, params.p, λ)
        simd_copy_yx_first!(λ, U, 2)

        # compute G=g(q,λ) at tₙ+₁
        params.g(t₁, q̅, λ, $tG)
        simd_copy_yx_first!($tG, G, 2)

        # compute p̅=α(q̅)
        params.α(t₁, q̅, λ, p̅)
    end

    return compute_projection_vprk
end


@generated function compute_internal_projection_vprk!(x::Vector{ST}, Q::Matrix{ST}, V::Matrix{ST}, Λ::Matrix{ST}, R::Matrix{ST}, Φ::Matrix{ST}, params::NonlinearFunctionParametersVPRKpSecondary{DT,TT,ΑT,FT,GT,VT,ΩT,HT,D,S}) where {ST,DT,TT,ΑT,FT,GT,VT,ΩT,HT,D,S}
    # create temporary vectors
    tQ = zeros(ST,D)
    tΛ = zeros(ST,D)
    tV = zeros(ST,D)
    tR = zeros(ST,D)
    tΦ = zeros(ST,D)
    tH = zeros(ST,D)
    tΩ = zeros(ST,D,D)

    compute_internal_projection_vprk = quote
        for i in 1:S
            simd_copy_xy_first!($tQ, Q, i)
            simd_copy_xy_first!($tV, V, i)
            simd_copy_xy_first!($tΛ, Λ, i)

            params.ω(params.t, $tQ, $tΩ)
            params.dH(params.t, $tQ, $tH)

            simd_mult!($tR, $tΩ, $tΛ)
            simd_mult!($tΦ, $tΩ, $tV)
            simd_axpy!(-one(ST), $tH, $tΦ)

            simd_copy_yx_first!($tR, R, i)
            simd_copy_yx_first!($tΦ, Φ, i)
        end
    end

    return compute_internal_projection_vprk
end


function compute_rhs_vprk_symplectic_projection!(b::Vector{ST}, p̅::Vector{ST}, F::Matrix{ST}, R::Matrix{ST}, G::Matrix{ST}, Φ::Matrix{ST}, offset::Int, params::AbstractNonlinearFunctionParametersVPRK{DT,TT,AT,FT,D,S}) where {ST,DT,TT,AT,FT,D,S}
    local z::ST

    # for i in 1:S
    #     for k in 1:D
    #         b[offset+D*(i-1)+k] = 0
    #         for j in 1:S
    #             b[offset+D*(i-1)+k] -= Φ[k,j]
    #         end
    #     end
    # end
    #
    # for i in 1:(S-1)
    #     for k in 1:D
    #         b[offset+D*(i-1)+k] = 0
    #         for j in 1:S
    #             b[offset+D*(i-1)+k] -= params.ω_λ[i,j] * Φ[k,j]
    #         end
    #     end
    # end

    for k in 1:D
        z = params.R[1] * G[k,1] + params.R[2] * G[k,2]
        for j in 1:S
            z += params.t_p.b[j] * (F[k,j] + R[k,j])
        end
        b[offset+D*(S-1)+k] = - (p̅[k] - params.p[k]) + params.Δt * z
    end
end


"Compute stages of variational partitioned Runge-Kutta methods."
@generated function function_stages!(x::Vector{ST}, b::Vector{ST}, params::NonlinearFunctionParametersVPRKpSecondary{DT,TT,ΑT,FT,GT,VT,ΩT,HT,D,S}) where {ST,DT,TT,ΑT,FT,GT,VT,ΩT,HT,D,S}
    scache = NonlinearFunctionCacheVPRK{ST}(D,S)
    pcache = NonlinearFunctionCacheVPRKprojection{ST}(D,S)

    function_stages = quote
        compute_stages_vprk!(x, $pcache.q̅, $pcache.p̅, $pcache.λ, $scache.Q, $scache.V, $pcache.Λ, $pcache.U, $scache.P, $scache.F, $pcache.R, $pcache.G, $pcache.Φ, params)

        # compute b = - [P-AF-U]
        compute_rhs_vprk!(b, $scache.P, $scache.F, $pcache.R, $pcache.G, params)

        # compute b = - [p-bF-G]
        compute_rhs_vprk_symplectic_projection!(b, $pcache.p̅, $scache.F, $pcache.R, $pcache.G, $pcache.Φ, D*S, params)

        compute_rhs_vprk_correction!(b, $scache.V, params)
    end

    return function_stages
end


"Variational partitioned Runge-Kutta integrator."
struct IntegratorVPRKpSecondary{DT,TT,ΑT,FT,GT,VT,ΩT,HT,FPT,ST,IT} <: AbstractIntegratorVPRK{DT,TT}
    equation::VODE{DT,TT,ΑT,FT,GT,VT,ΩT,HT}
    tableau::TableauVPRK{TT}
    Δt::TT

    params::FPT
    solver::ST

    scache::NonlinearFunctionCacheVPRK{DT}
    pcache::NonlinearFunctionCacheVPRKprojection{DT}

    iguess::InitialGuessPODE{DT,TT,VT,FT,IT}

    q::Array{DT,1}
    p::Array{DT,1}

    qₑᵣᵣ::Vector{DT}
    pₑᵣᵣ::Vector{DT}
end

function IntegratorVPRKpSecondary(equation::VODE{DT,TT,ΑT,FT,GT,VT,ΩT,HT}, tableau::TableauVPRK{TT}, Δt::TT;
                                        nonlinear_solver=DEFAULT_NonlinearSolver,
                                        nmax=DEFAULT_nmax, atol=DEFAULT_atol, rtol=DEFAULT_rtol, stol=DEFAULT_stol,
                                        interpolation=HermiteInterpolation{DT}) where {DT,TT,ΑT,FT,GT,VT,ΩT,HT}
    D = equation.d
    S = tableau.s

    N = D*S*2

    if isdefined(tableau, :d)
        tableau_d = tableau.d
    else
        tableau_d = DT[]
    end

    # compute reduction matrix
    ω_λ = zeros(TT, S-1, S)

    for i in 1:(S-1)
        for j in 1:S
            ω_λ[i,j] = tableau.b[j] * tableau.c[j]^(i-1)
        end
    end

    # create solution vectors
    q = zeros(DT,D)
    p = zeros(DT,D)

    # create compensated summation error vectors
    qₑᵣᵣ = zeros(DT,D)
    pₑᵣᵣ = zeros(DT,D)

    # create cache for internal stage vectors and update vectors
    scache = NonlinearFunctionCacheVPRK{DT}(D,S)
    pcache = NonlinearFunctionCacheVPRKprojection{DT}(D,S)

    # create params
    params = NonlinearFunctionParametersVPRKpSecondary{DT,TT,ΑT,FT,GT,VT,ΩT,HT,D,S}(
                                                equation.α, equation.f, equation.g, equation.v, equation.ω, equation.dH,
                                                Δt, tableau.o, tableau.q, tableau.p, tableau_d, ω_λ, tableau.R∞,
                                                q, p)

    # create rhs function for nonlinear solver
    function_stages_solver = (x,b) -> function_stages!(x, b, params)

    # create solver
    solver = nonlinear_solver(zeros(DT,N), function_stages_solver; nmax=nmax, atol=atol, rtol=rtol, stol=stol)

    # create initial guess
    iguess = InitialGuessPODE(interpolation, equation, Δt; periodicity=equation.periodicity)

    IntegratorVPRKpSecondary{DT, TT, ΑT, FT, GT, VT, ΩT, HT, typeof(params), typeof(solver), typeof(iguess.int)}(
                                        equation, tableau, Δt, params, solver, scache, pcache, iguess,
                                        q, p, qₑᵣᵣ, pₑᵣᵣ)
end


"Integrate ODE with variational partitioned Runge-Kutta integrator."
function integrate!(int::IntegratorVPRKpSecondary{DT,TT,ΑT,FT,GT,VT,ΩT,HT}, sol::SolutionPDAE{DT,TT,N}, m1::Int, m2::Int) where {DT,TT,ΑT,FT,GT,VT,N,ΩT,HT}
    @assert m1 ≥ 1
    @assert m2 ≤ sol.ni

    # loop over initial conditions
    for m in m1:m2
        # copy initial conditions from solution
        get_initial_conditions!(sol, int.q, int.p, m)

        # initialise initial guess
        initialize!(int.iguess, sol.t[0], int.q, int.p)

        for n in 1:sol.ntime
            # set time for nonlinear solver
            int.params.t = sol.t[0] + (n-1)*int.Δt

            # compute initial guess
            for i in 1:int.tableau.s
                evaluate!(int.iguess, int.scache.y, int.scache.z, int.scache.v, int.tableau.q.c[i], int.tableau.p.c[i])
                for k in 1:int.equation.d
                    int.solver.x[int.equation.d*(0*int.tableau.s+i-1)+k] = int.scache.v[k]
                    int.solver.x[int.equation.d*(1*int.tableau.s+i-1)+k] = 0
                end
            end

            # call nonlinear solver
            solve!(int.solver; refactorize=1)

            if !check_solver_status(int.solver.status, int.solver.params)
                println(int.solver.status, ", it=", n)
            end

            if isnan(int.solver.status.rₐ)
                break
            end

            # compute final update
            compute_stages_vprk!(int.solver.x,
                                 int.pcache.q̅, int.pcache.p̅, int.pcache.λ,
                                 int.scache.Q, int.scache.V, int.pcache.Λ, int.pcache.U,
                                 int.scache.P, int.scache.F, int.pcache.R, int.pcache.G,
                                 int.pcache.Φ, int.params)

            update_solution!(int, int.scache)
            project_solution!(int, int.pcache, int.params.R)

            # copy solution to initial guess for next time step
            update!(int.iguess, sol.t[0] + n*int.Δt, int.q, int.p)

            # take care of periodic solutions
            cut_periodic_solution!(int)

            # copy to solution
            copy_solution!(sol, int.q, int.p, int.pcache.λ, n, m)
        end
    end
end
