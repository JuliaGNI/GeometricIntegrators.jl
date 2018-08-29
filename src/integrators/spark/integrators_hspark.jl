
"Holds the tableau of an Hamiltonian Specialised Partitioned Additive Runge-Kutta method."
struct TableauHSPARK{T} <: AbstractTableau{T}
end

"Parameters for right-hand side function of Hamiltonian Specialised Partitioned Additive Runge-Kutta methods."
mutable struct ParametersHSPARK{DT,TT,VT,FT,ϕT,ψT} <: Parameters{DT,TT}
    f_v::VT
    f_f::FT
    f_ϕ::ϕT
    f_ψ::ψT

    Δt::TT

    d::Int
    m::Int
    s::Int
    σ::Int

    t_q::CoefficientsARK{TT}
    t_p::CoefficientsARK{TT}
    t_q̃::CoefficientsPRK{TT}
    t_p̃::CoefficientsPRK{TT}
    t_ω::Matrix{TT}

    t::TT

    q::Vector{DT}
    p::Vector{DT}

    function ParametersHSPARK{DT,TT,VT,FT,ϕT,ψT}(f_v, f_f, f_ϕ, f_ψ, Δt, d, m, s, σ, t_q, t_p, t_q̃, t_p̃, t_ω) where {DT,TT,VT,FT,ϕT,ψT}
        # create solution vectors
        q = zeros(DT,d)
        p = zeros(DT,d)

        new(f_v, f_f, f_ϕ, f_ψ,
            Δt, d, s, σ,
            t_q, t_p, t_q̃, t_p̃, t_ω,
            0, q, p)
    end
end


struct NonlinearFunctionCacheHSPARK{ST}

    q̅::Vector{ST}
    p̅::Vector{ST}

    y::Vector{ST}
    z::Vector{ST}

    Qi::Matrix{ST}
    Pi::Matrix{ST}
    Λi::Matrix{ST}

    Vi::Matrix{ST}
    Fi::Matrix{ST}
    Yi::Matrix{ST}
    Zi::Matrix{ST}
    Φi::Matrix{ST}

    Qp::Matrix{ST}
    Pp::Matrix{ST}
    Λp::Matrix{ST}
    Γp::Matrix{ST}

    Vp::Matrix{ST}
    Fp::Matrix{ST}
    Yp::Matrix{ST}
    Zp::Matrix{ST}
    Φp::Matrix{ST}
    Ψp::Matrix{ST}

    # Qt::Vector{ST}
    # Pt::Vector{ST}
    # Λt::Vector{ST}
    # Γt::Vector{ST}
    #
    # Vt::Vector{ST}
    # Ft::Vector{ST}
    # Φt::Vector{ST}
    # Ψt::Vector{ST}

    function NonlinearFunctionCacheHSPARK{ST}(d,m,s,σ) where {ST}
        # create solution vectors
        q̅ = zeros(ST,d)
        p̅ = zeros(ST,d)

        y = zeros(ST,d)
        z = zeros(ST,d)

        # create internal stage vectors
        Qi = zeros(DT,d,s)
        Pi = zeros(DT,d,s)
        Λi = zeros(DT,d,s)

        Vi = zeros(DT,d,s)
        Fi = zeros(DT,d,s)
        Yi = zeros(DT,d,s)
        Zi = zeros(DT,d,s)
        Φi = zeros(DT,m,s)

        Qp = zeros(DT,d,σ)
        Pp = zeros(DT,d,σ)
        Λp = zeros(DT,m,σ)
        Γp = zeros(DT,m,σ)

        Vp = zeros(DT,d,σ)
        Fp = zeros(DT,d,σ)
        Yp = zeros(DT,d,σ)
        Zp = zeros(DT,d,σ)
        Φp = zeros(DT,m,σ)
        Ψp = zeros(DT,m,σ)

        # create temporary vectors
        # Qt = zeros(DT,d)
        # Pt = zeros(DT,d)
        # Λt = zeros(DT,m)
        # Γt = zeros(DT,m)
        #
        # Vt = zeros(DT,d)
        # Ft = zeros(DT,d)
        # Φt = zeros(DT,m)
        # Ψt = zeros(DT,m)

        new(q̅, p̅, y, z, Qi, Pi, Λi, Vi, Fi, Yi, Zi, Φi, Qp, Pp, Λp, Γp, Vp, Fp, Yp, Zp, Φp, Ψp)
    end
end


# -> edited till here


"Compute stages of Hamiltonian Specialised Partitioned Additive Runge-Kutta methods."
function function_stages!(y::Vector{DT}, b::Vector{DT}, params::ParametersHSPARK{DT,TT,VT,FT,ϕT,ψT}) where {DT,TT,VT,FT,ϕT,ψT}
    local offset::Int
    local tpᵢ::TT
    local tλᵢ::TT

    for i in 1:params.s
        for k in 1:params.d
            # copy y to Y, Z
            offset = 3*(params.d*(i-1)+k-1)
            params.Yi[k,i] = y[offset+1]
            params.Zi[k,i] = y[offset+2]
            params.Vi[k,i] = y[offset+3]

            # compute Q and P
            params.Qi[k,i] = params.q[k] + params.Δt * params.Yi[k,i]
            params.Pi[k,i] = params.p[k] + params.Δt * params.Zi[k,i]
        end

        # compute f(X)
        tpᵢ = params.t + params.Δt * params.t_p.c[i]

        simd_copy_xy_first!(params.Qt, params.Qi, i)
        # simd_copy_xy_first!(params.Pt, params.Pi, i)
        simd_copy_xy_first!(params.Vt, params.Vi, i)
        params.f_f(tpᵢ, params.Qt, params.Vt, params.Ft)
        params.f_p(tpᵢ, params.Qt, params.Vt, params.Pt)
        # params.f_ϕ(tpᵢ, params.Qt, params.Pt, params.Φt)
        simd_copy_yx_first!(params.Ft, params.Fi, i)
        # simd_copy_yx_first!(params.Φt, params.Φi, i)
        simd_copy_yx_first!(params.Pt, params.Φi, i)
    end

    params.Φi .-= params.Pi

    for i in 1:params.r
        for k in 1:params.d
            # copy y to Y and Z
            offset = 3*params.d*params.s+3*(params.d*(i-1)+k-1)
            params.Yp[k,i] = y[offset+1]
            params.Zp[k,i] = y[offset+2]
            params.Vp[k,i] = y[offset+3]

            # compute Q and V
            params.Qp[k,i] = params.q[k] + params.Δt * params.Yp[k,i]
            params.Pp[k,i] = params.p[k] + params.Δt * params.Zp[k,i]
        end

        # compute f(X)
        tλᵢ = params.t + params.Δt * params.t_λ.c[i]

        simd_copy_xy_first!(params.Qt, params.Qp, i)
        simd_copy_xy_first!(params.Pt, params.Pp, i)
        simd_copy_xy_first!(params.Vt, params.Vp, i)
        params.f_u(tλᵢ, params.Qt, params.Pt, params.Vt, params.Ut)
        params.f_g(tλᵢ, params.Qt, params.Pt, params.Vt, params.Gt)
        params.f_ϕ(tλᵢ, params.Qt, params.Pt, params.Φt)
        simd_copy_yx_first!(params.Ut, params.Up, i)
        simd_copy_yx_first!(params.Gt, params.Gp, i)
        simd_copy_yx_first!(params.Φt, params.Φp, i)
    end

    for i in 1:params.ρ
        offset = 3*params.d*params.s+3*params.d*params.r+params.d*(i-1)
        for k in 1:params.d
            params.Λρ[k,i] = y[offset+k]
        end
    end

    if length(params.d_v) > 0
        offset = 3*params.d*params.s+3*params.d*params.r+params.d*params.ρ
        for k in 1:params.d
            params.μ[k] = y[offset+k]
        end
    end

    # compute b = - [(Y-AV-AU), (Z-AF-AG), Φ]
    for i in 1:params.s
        for k in 1:params.d
            offset = 3*(params.d*(i-1)+k-1)
            b[offset+1] = - params.Yi[k,i]
            b[offset+2] = - params.Zi[k,i]
            b[offset+3] = params.Φi[k,i]
            for j in 1:params.s
                b[offset+1] += params.t_q.a[i,j] * params.Vi[k,j]
                b[offset+2] += params.t_p.a[i,j] * params.Fi[k,j]
            end
            for j in 1:params.r
                b[offset+1] += params.t_q.α[i,j] * params.Up[k,j]
                b[offset+2] += params.t_p.α[i,j] * params.Gp[k,j]
            end
        end
    end

    # compute b = - [(Y-AV-AU), (Z-AF-AG), Φ]
    for i in 1:params.r
        for k in 1:params.d
            offset = 3*params.d*params.s+3*(params.d*(i-1)+k-1)
            b[offset+1] = - params.Yp[k,i]
            b[offset+2] = - params.Zp[k,i]
            b[offset+3] = - params.Vp[k,i]
            for j in 1:params.s
                b[offset+1] += params.t_q̃.a[i,j] * params.Vi[k,j]
                b[offset+2] += params.t_p̃.a[i,j] * params.Fi[k,j]
            end
            for j in 1:params.r
                b[offset+1] += params.t_q̃.α[i,j] * params.Up[k,j]
                b[offset+2] += params.t_p̃.α[i,j] * params.Gp[k,j]
            end
            for j in 1:params.ρ
                b[offset+3] += params.ω_λ[j,i] * params.Λρ[k,j]
            end
        end
    end

    for i in 1:params.ρ
        for k in 1:params.d
            offset = 3*params.d*params.s+3*params.d*params.r+params.d*(i-1)+k
            b[offset] = 0
            for j in 1:params.r
                b[offset] -= params.ω_λ[i,j] * params.Φp[k,j]
            end
        end
    end

    if length(params.d_v) > 0
        for i in 1:params.s
            for k in 1:params.d
                b[3*(params.d*(i-1)+k-1)+3] -= params.d_v[i] * params.μ[k]
            end
        end

        offset = 3*params.d*params.s+3*params.d*params.r+params.d*params.ρ
        for k in 1:params.d
            b[offset+k] = 0
            for i in 1:params.s
                b[offset+k] -= params.Vi[k,i] * params.d_v[i]
            end
        end
    end
end


"Hamiltonian Specialised Partitioned Additive Runge-Kutta integrator."
mutable struct IntegratorHSPARK{DT, TT, VT, FT, ϕT, ψT, SPT, ST, IT} <: Integrator{DT, TT}
    equation::IDAE{DT,TT,VT,FT,ϕT,ψT}
    tableau::TableauHSPARK{TT}
    Δt::TT

    params::SPT
    solver::ST
    iguess::InitialGuessPODE{DT, TT, VT, FT, IT}

    q::Array{DT,1}
    v::Array{DT,1}
    p::Array{DT,1}
    λ::Array{DT,1}
    y::Array{DT,1}
    z::Array{DT,1}

    Q::Array{DT,2}
    P::Array{DT,2}
    V::Array{DT,2}
    F::Array{DT,2}
    Λ::Array{DT,2}
    U::Array{DT,2}
    G::Array{DT,2}
end

function IntegratorHSPARK(equation::HDAE{DT,TT,VT,FT,ϕT,ψT}, tableau::TableauHSPARK{TT}, Δt::TT) where {DT,TT,VT,FT,ϕT,ψT}
    D = equation.d
    S = tableau.s
    R = tableau.r
    ρ = tableau.ρ

    N = 3*D*S + 3*D*R + D*ρ

    if isdefined(tableau, :d)
        N += D
        d_v = tableau.d
    else
        d_v = DT[]
    end

    # create solution vector for internal stages / nonlinear solver
    z = zeros(DT, N)

    # create params
    params = ParametersHSPARK{DT,TT,VT,FT,ϕT,ψT}(
                                                equation.f, equation.p, equation.u, equation.g, equation.ϕ,
                                                Δt, D, S, R, ρ,
                                                tableau.q, tableau.p, tableau.q̃, tableau.p̃, tableau.λ, tableau.ω, d_v)

    # create rhs function for nonlinear solver
    function_stages = (x,b) -> function_stages!(x, b, params)

    # create solver
    solver = get_config(:nls_solver)(z, function_stages)

    # create initial guess
    iguess = InitialGuessPODE(get_config(:ig_interpolation), equation, Δt)

    # create integrator
    IntegratorHSPARK{DT, TT, FT, PT, UT, GT, ϕT, VT, typeof(params), typeof(solver), typeof(iguess.int)}(
                                        equation, tableau, Δt, params, solver, iguess,
                                        params.q, params.v, params.p, params.λ,
                                        params.y, params.z,
                                        params.Qi, params.Pi, params.Vi, params.Fi,
                                        params.Vp, params.Up, params.Gp)
end


function initialize!(int::IntegratorHSPARK, sol::SolutionPDAE, m::Int)
    @assert m ≥ 1
    @assert m ≤ sol.ni

    # copy initial conditions from solution
    get_initial_conditions!(sol, int.q, int.p, int.λ, m)

    # initialise initial guess
    initialize!(int.iguess, m, sol.t[0], int.q, int.p)
end

"Integrate DAE with variational special partitioned additive Runge-Kutta integrator."
function integrate_step!(int::IntegratorHSPARK{DT,TT,VT,FT,ϕT,ψT}, sol::SolutionPDAE{DT,TT,N}, m::Int, n::Int) where {DT,TT,VT,FT,ϕT,ψT,N}
    local offset::Int

    # set time for nonlinear solver
    int.params.t = sol.t[n]

    # copy previous solution to initial guess
    update!(int.iguess, sol.t[n], int.q, int.p)

    # compute initial guess
    for i in 1:int.tableau.s
        evaluate!(int.iguess, int.y, int.z, int.v, int.tableau.q.c[i], int.tableau.p.c[i])
        for k in 1:int.equation.d
            offset = 3*(int.equation.d*(i-1)+k-1)
            int.solver.x[offset+1] = (int.y[k] - int.q[k])/(int.Δt)
            int.solver.x[offset+2] = (int.z[k] - int.p[k])/(int.Δt)
            int.solver.x[offset+3] = int.v[k]
        end
    end

    for i in 1:int.tableau.r
        evaluate!(int.iguess, int.y, int.z, int.v, int.tableau.q̃.c[i], int.tableau.p̃.c[i])
        for k in 1:int.equation.d
            offset = 3*int.equation.d*int.tableau.s+3*(int.equation.d*(i-1)+k-1)
            int.solver.x[offset+1] = (int.y[k] - int.q[k])/(int.Δt)
            int.solver.x[offset+2] = (int.z[k] - int.p[k])/(int.Δt)
            int.solver.x[offset+3] = 0
        end
    end

    for i in 1:int.tableau.ρ
        offset = 3*int.equation.d*int.tableau.s+3*int.equation.d*int.tableau.r+int.equation.d*(i-1)
        for k in 1:int.equation.d
            int.solver.x[offset+k] = 0
            # for j in 1:int.tableau.r
            #     int.solver.x[offset+k] += 0
            # end
        end
    end

    if isdefined(int.tableau, :d)
        offset = 3*int.equation.d*int.tableau.s+3*int.equation.d*int.tableau.r+int.equation.d*int.tableau.ρ
        for k in 1:int.equation.d
            int.solver.x[offset+k] = 0
        end
    end

    # println("ig:  ", join([@sprintf "%8.4e" x for x in int.solver.x], ", "))

    # call nonlinear solver
    solve!(int.solver)

    # println("sol: ", join([@sprintf "%8.4e" x for x in int.solver.x], ", "))
    #
    # println(int.solver.status)
    # println()

    if !check_solver_status(int.solver.status, int.solver.params)
        println(int.solver.status)
    end

    # if isnan(int.solver.status.rₐ)
    #     break
    # end

    # compute final update
    simd_mult!(int.y, int.V, int.tableau.q.b)
    simd_mult!(int.z, int.F, int.tableau.p.b)
    int.q .+= int.Δt .* int.y
    int.p .+= int.Δt .* int.z
    simd_mult!(int.y, int.U, int.tableau.q.β)
    simd_mult!(int.z, int.G, int.tableau.p.β)
    int.q .+= int.Δt .* int.y
    int.p .+= int.Δt .* int.z
    simd_mult!(int.λ, int.Λ, int.tableau.λ.b)

    # copy to solution
    copy_solution!(sol, int.q, int.p, int.λ, n, m)
end
