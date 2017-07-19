
"Parameters for right-hand side function of partitioned additive Runge-Kutta methods."
mutable struct NonlinearFunctionParametersPARK{DT,TT,FT,PT,UT,GT,ϕT} <: NonlinearFunctionParameters{DT}
    f_f::FT
    f_p::PT
    f_u::UT
    f_g::GT
    f_ϕ::ϕT

    Δt::TT

    d::Int
    s::Int
    r::Int

    t_q::CoefficientsARK{TT}
    t_p::CoefficientsARK{TT}
    t_q̃::CoefficientsPRK{TT}
    t_p̃::CoefficientsPRK{TT}
    t_λ::CoefficientsMRK{TT}

    t::TT

    q::Vector{DT}
    p::Vector{DT}
    λ::Vector{DT}

    y::Vector{DT}
    z::Vector{DT}

    Qi::Matrix{DT}
    Pi::Matrix{DT}
    Λi::Matrix{DT}
    Vi::Matrix{DT}
    Fi::Matrix{DT}
    Yi::Matrix{DT}
    Zi::Matrix{DT}
    Φi::Matrix{DT}

    Qp::Matrix{DT}
    Pp::Matrix{DT}
    Λp::Matrix{DT}
    Up::Matrix{DT}
    Gp::Matrix{DT}
    Yp::Matrix{DT}
    Zp::Matrix{DT}
    Φp::Matrix{DT}

    Qt::Vector{DT}
    Pt::Vector{DT}
    Λt::Vector{DT}
    Vt::Vector{DT}
    Ft::Vector{DT}
    Ut::Vector{DT}
    Gt::Vector{DT}
    Φt::Vector{DT}

    function NonlinearFunctionParametersPARK{DT,TT,FT,PT,UT,GT,ϕT}(f_f, f_p, f_u, f_g, f_ϕ, Δt, d, s, r, t_q, t_p, t_q̃, t_p̃, t_λ) where {DT,TT,FT,PT,UT,GT,ϕT}
        # create solution vectors
        q = zeros(DT,d)
        p = zeros(DT,d)
        λ = zeros(DT,d)
        y = zeros(DT,d)
        z = zeros(DT,d)

        # create internal stage vectors
        Qi = zeros(DT,d,s)
        Pi = zeros(DT,d,s)
        Λi = zeros(DT,d,s)
        Vi = zeros(DT,d,s)
        Fi = zeros(DT,d,s)
        Yi = zeros(DT,d,s)
        Zi = zeros(DT,d,s)
        Φi = zeros(DT,d,s)

        Qp = zeros(DT,d,r)
        Pp = zeros(DT,d,r)
        Λp = zeros(DT,d,r)
        Up = zeros(DT,d,r)
        Gp = zeros(DT,d,r)
        Yp = zeros(DT,d,r)
        Zp = zeros(DT,d,r)
        Φp = zeros(DT,d,r)

        # create temporary vectors
        Qt = zeros(DT,d)
        Pt = zeros(DT,d)
        Λt = zeros(DT,d)
        Vt = zeros(DT,d)
        Ft = zeros(DT,d)
        Ut = zeros(DT,d)
        Gt = zeros(DT,d)
        Φt = zeros(DT,d)

        new(f_f, f_p, f_u, f_g, f_ϕ, Δt, d, s, r,
            t_q, t_p, t_q̃, t_p̃, t_λ,
            0, q, p, λ, y, z,
            Qi, Pi, Λi, Vi, Fi, Yi, Zi, Φi,
            Qp, Pp, Λp, Up, Gp, Yp, Zp, Φp,
            Qt, Pt, Λt, Vt, Ft, Ut, Gt, Φt)
    end
end

"Compute stages of partitioned additive Runge-Kutta methods."
function function_stages!(y::Vector{DT}, b::Vector{DT}, params::NonlinearFunctionParametersPARK{DT,TT,FT,PT,UT,GT,ϕT}) where {DT,TT,FT,PT,UT,GT,ϕT}
    local tpᵢ::TT
    local tλᵢ::TT

    for i in 1:params.s
        for k in 1:params.d
            # copy y to Y, Z, Λ
            params.Yi[k,i] = y[2*(params.d*(i-1)+k-1)+1]
            params.Zi[k,i] = y[2*(params.d*(i-1)+k-1)+2]

            # compute Q and P
            params.Qi[k,i] = params.q[k] + params.Δt * params.Yi[k,i]
            params.Pi[k,i] = params.p[k] + params.Δt * params.Zi[k,i]
        end

        # compute f(X)
        tpᵢ = params.t + params.Δt * params.t_p.c[i]

        simd_copy_xy_first!(params.Qt, params.Qi, i)
        simd_copy_xy_first!(params.Pt, params.Pi, i)
        params.f_v(tpᵢ, params.Qt, params.Pt, params.Vt)
        params.f_f(tpᵢ, params.Qt, params.Pt, params.Ft)
        simd_copy_yx_first!(params.Vt, params.Vi, i)
        simd_copy_yx_first!(params.Ft, params.Fi, i)
    end

    for i in 1:params.r
        for k in 1:params.d
            # copy y to Y and Z
            params.Yp[k,i] = y[2*params.d*params.s+3*(params.d*(i-1)+k-1)+1]
            params.Zp[k,i] = y[2*params.d*params.s+3*(params.d*(i-1)+k-1)+2]
            params.Λp[k,i] = y[2*params.d*params.s+3*(params.d*(i-1)+k-1)+3]

            # compute Q and V
            params.Qp[k,i] = params.q[k] + params.Δt * params.Yp[k,i]
            params.Pp[k,i] = params.p[k] + params.Δt * params.Zp[k,i]
        end

        # compute f(X)
        tλᵢ = params.t + params.Δt * params.t_q.c[i]

        simd_copy_xy_first!(params.Qt, params.Qp, i)
        simd_copy_xy_first!(params.Pt, params.Pp, i)
        simd_copy_xy_first!(params.Λt, params.Λp, i)
        params.f_u(tλᵢ, params.Qt, params.Pt, params.Λt, params.Ut)
        params.f_g(tλᵢ, params.Qt, params.Pt, params.Λt, params.Gt)
        params.f_ϕ(tλᵢ, params.Qt, params.Pt, params.Φt)
        simd_copy_yx_first!(params.Ut, params.Up, i)
        simd_copy_yx_first!(params.Gt, params.Gp, i)
        simd_copy_yx_first!(params.Φt, params.Φp, i)
    end

    # compute b = - [(Y-AV-AU), (Z-AF-AG), Φ]
    for i in 1:params.s
        for k in 1:params.d
            b[2*(params.d*(i-1)+k-1)+1] = - params.Yi[k,i]
            b[2*(params.d*(i-1)+k-1)+2] = - params.Zi[k,i]
            for j in 1:params.s
                b[2*(params.d*(i-1)+k-1)+1] += params.t_q.a[i,j] * params.Vi[k,j]
                b[2*(params.d*(i-1)+k-1)+2] += params.t_p.a[i,j] * params.Fi[k,j]
            end
            for j in 1:params.r
                b[2*(params.d*(i-1)+k-1)+1] += params.t_q̃.a[i,j] * params.Up[k,j]
                b[2*(params.d*(i-1)+k-1)+2] += params.t_p̃.a[i,j] * params.Gp[k,j]
            end
        end
    end

    # compute b = - [(Y-AV-AU), (Z-AF-AG), Φ]
    for i in 1:params.r
        for k in 1:params.d
            b[2*params.d*params.s+3*(params.d*(i-1)+k-1)+1] = - params.Yp[k,i]
            b[2*params.d*params.s+3*(params.d*(i-1)+k-1)+2] = - params.Zp[k,i]
            b[2*params.d*params.s+3*(params.d*(i-1)+k-1)+3] = - params.Φp[k,i]
            for j in 1:params.s
                b[2*params.d*params.s+3*(params.d*(i-1)+k-1)+1] += params.t_q.α[i,j] * params.Vi[k,j]
                b[2*params.d*params.s+3*(params.d*(i-1)+k-1)+2] += params.t_p.α[i,j] * params.Fi[k,j]
            end
            for j in 1:params.r
                b[2*params.d*params.s+3*(params.d*(i-1)+k-1)+1] += params.t_q̃.α[i,j] * params.Up[k,j]
                b[2*params.d*params.s+3*(params.d*(i-1)+k-1)+2] += params.t_p̃.α[i,j] * params.Gp[k,j]
            end
        end
    end

    # compute b = - [Λ₁-λ]
    if params.t_λ.c[1] == 0
        for k in 1:params.d
            b[2*params.d*params.s+3*(k-1)+3] = - params.Λp[k,1] + params.λ[k]
        end
    end
end


"Implicit partitioned additive Runge-Kutta integrator."
struct IntegratorPARK{DT, TT, FT, PT, UT, GT, ϕT, ST} <: Integrator{DT, TT}
    equation::IDAE{DT,TT,FT,PT,UT,GT,ϕT}
    tableau::TableauPARK{TT}
    Δt::TT

    solver::ST

    q::Array{DT,1}
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

function IntegratorPARK(equation::IDAE{DT,TT,FT,PT,UT,GT,ϕT}, tableau::TableauPARK{TT}, Δt::TT;
                        nonlinear_solver=DEFAULT_NonlinearSolver,
                        nmax=DEFAULT_nmax, atol=DEFAULT_atol, rtol=DEFAULT_rtol, stol=DEFAULT_stol,
                        interpolation=HermiteInterpolation{DT}) where {DT,TT,FT,PT,UT,GT,ϕT}
    D = equation.d
    S = tableau.s
    R = tableau.r

    # create solution vector for internal stages / nonlinear solver
    z = zeros(DT, 2*D*S + 3*D*R)

    # create params
    params = NonlinearFunctionParametersPARK{DT,TT,FT,PT,UT,GT,ϕT}(
                                                equation.f, equation.p, equation.u, equation.g, equation.ϕ,
                                                Δt, D, S, R, tableau.q, tableau.p, tableau.q̃, tableau.p̃, tableau.λ)

    # create solver
    solver = nonlinear_solver(z, params; nmax=nmax, atol=atol, rtol=rtol, stol=stol)

    # create integrator
    IntegratorPARK{DT, TT, FT, PT, UT, GT, ϕT, typeof(solver)}(
                                        equation, tableau, Δt, solver,
                                        params.q, params.p, params.λ, params.y, params.z,
                                        params.Qi, params.Pi, params.Vi, params.Fi,
                                        params.Λp, params.Up, params.Gp)
end


function initialize!(int::IntegratorPARK, sol::Union{SolutionPDAE, PSolutionPDAE}, m::Int)
    @assert m ≥ 1
    @assert m ≤ sol.ni

    # copy initial conditions from solution
    get_initial_conditions!(sol, int.q, int.p, int.λ, m)
end

"Integrate DAE with partitioned additive Runge-Kutta integrator."
function integrate_step!(int::IntegratorPARK{DT,TT,FT,PT,UT,GT,ϕT}, sol::SolutionPDAE{DT,TT,N}, m::Int, n::Int) where {DT,TT,FT,PT,UT,GT,ϕT,N}
    # set time for nonlinear solver
    int.solver.Fparams.t = sol.t[n]

    # compute initial guess
    for i in 1:int.tableau.s
        for k in 1:int.equation.d
            # TODO initial guess for y and z
            int.solver.x[2*(int.equation.d*(i-1)+k-1)+1] = 0
            int.solver.x[2*(int.equation.d*(i-1)+k-1)+2] = 0
        end
    end

    for i in 1:int.tableau.r
        for k in 1:int.equation.d
            # TODO initial guess for y and z
            int.solver.x[2*int.equation.d*int.tableau.s+3*(int.equation.d*(i-1)+k-1)+1] = 0
            int.solver.x[2*int.equation.d*int.tableau.s+3*(int.equation.d*(i-1)+k-1)+2] = 0
            int.solver.x[2*int.equation.d*int.tableau.s+3*(int.equation.d*(i-1)+k-1)+3] = 0
        end
    end

    # call nonlinear solver
    solve!(int.solver)

    if !check_solver_status(int.solver.status, int.solver.params)
        println(int.solver.status)
    end

    # compute final update
    simd_mult!(int.y, int.V, int.tableau.q.b)
    simd_mult!(int.z, int.F, int.tableau.p.b)
    simd_axpy!(int.Δt, int.y, int.q)
    simd_axpy!(int.Δt, int.z, int.p)
    simd_mult!(int.y, int.U, int.tableau.q.β)
    simd_mult!(int.z, int.G, int.tableau.p.β)
    simd_axpy!(int.Δt, int.y, int.q)
    simd_axpy!(int.Δt, int.z, int.p)
    simd_mult!(int.λ, int.Λ, int.tableau.λ.b)

    # copy to solution
    copy_solution!(sol, int.q, int.p, int.λ, n, m)
end
