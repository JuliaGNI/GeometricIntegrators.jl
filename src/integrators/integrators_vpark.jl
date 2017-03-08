
"Parameters for right-hand side function of variational partitioned additive Runge-Kutta methods."
type NonlinearFunctionParametersVPARK{DT,TT,FT,PT,UT,GT,ϕT} <: NonlinearFunctionParameters{DT}
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
    d_v::Vector{TT}

    t::TT

    q::Vector{DT}
    v::Vector{DT}
    p::Vector{DT}
    λ::Vector{DT}
    μ::Vector{DT}

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

    function NonlinearFunctionParametersVPARK(f_f, f_p, f_u, f_g, f_ϕ, Δt, d, s, r, t_q, t_p, t_q̃, t_p̃, t_λ, d_v)
        # create solution vectors
        q = zeros(DT,d)
        v = zeros(DT,d)
        p = zeros(DT,d)
        λ = zeros(DT,d)
        μ = zeros(DT,d)

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
            t_q, t_p, t_q̃, t_p̃, t_λ, d_v,
            0, q, v, p, λ, μ, y, z,
            Qi, Pi, Λi, Vi, Fi, Yi, Zi, Φi,
            Qp, Pp, Λp, Up, Gp, Yp, Zp, Φp,
            Qt, Pt, Λt, Vt, Ft, Ut, Gt, Φt)
    end
end

"Compute stages of variational partitioned additive Runge-Kutta methods."
function function_stages!{DT,TT,FT,PT,UT,GT,ϕT}(y::Vector{DT}, b::Vector{DT}, params::NonlinearFunctionParametersVPARK{DT,TT,FT,PT,UT,GT,ϕT})
    local tpᵢ::TT
    local tλᵢ::TT

    for i in 1:params.s
        for k in 1:params.d
            # copy y to Y, Z
            params.Yi[k,i] = y[2*(params.d*(i-1)+k-1)+1]
            params.Zi[k,i] = y[2*(params.d*(i-1)+k-1)+2]

            # compute Q and V
            params.Qi[k,i] = params.q[k] + params.Δt * params.Yi[k,i]
            params.Vi[k,i] = params.Zi[k,i]
        end

        # compute f(X)
        tpᵢ = params.t + params.Δt * params.t_p.c[i]

        simd_copy_xy_first!(params.Qt, params.Qi, i)
        simd_copy_xy_first!(params.Vt, params.Vi, i)
        params.f_f(tpᵢ, params.Qt, params.Vt, params.Ft)
        params.f_p(tpᵢ, params.Qt, params.Vt, params.Pt)
        simd_copy_yx_first!(params.Ft, params.Fi, i)
        simd_copy_yx_first!(params.Pt, params.Pi, i)
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
        tλᵢ = params.t + params.Δt * params.t_λ.c[i]

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

    if length(params.d_v) > 0
        for k in 1:params.d
            params.μ[k] = y[2*params.d*params.s+3*params.d*params.r+k]
        end
    end

    # compute b = - [(Y-AV-AU), (Z-AF-AG), Φ]
    for i in 1:params.s
        for k in 1:params.d
            b[2*(params.d*(i-1)+k-1)+1] = - params.Yi[k,i]
            b[2*(params.d*(i-1)+k-1)+2] = - params.Pi[k,i] + params.p[k]
            for j in 1:params.s
                b[2*(params.d*(i-1)+k-1)+1] += params.t_q.a[i,j] * params.Vi[k,j]
                b[2*(params.d*(i-1)+k-1)+2] += params.t_p.a[i,j] * params.Fi[k,j] * params.Δt
            end
            for j in 1:params.r
                b[2*(params.d*(i-1)+k-1)+1] += params.t_q.α[i,j] * params.Up[k,j]
                b[2*(params.d*(i-1)+k-1)+2] += params.t_p.α[i,j] * params.Gp[k,j] * params.Δt
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
                b[2*params.d*params.s+3*(params.d*(i-1)+k-1)+1] += params.t_q̃.a[i,j] * params.Vi[k,j]
                b[2*params.d*params.s+3*(params.d*(i-1)+k-1)+2] += params.t_p̃.a[i,j] * params.Fi[k,j]
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

    if length(params.d_v) > 0
        for i in 1:params.s
            for k in 1:params.d
                b[2*(params.d*(i-1)+k-1)+2] -= params.d_v[i] * params.μ[k]
            end
        end

        for k in 1:params.d
            b[2*params.d*params.s+3*params.d*params.r+k] = 0
            for i in 1:params.s
                b[2*params.d*params.s+3*params.d*params.r+k] -= params.Vi[k,i] * params.d_v[i]
            end
        end
    end
end


"Variational partitioned additive Runge-Kutta integrator."
immutable IntegratorVPARK{DT, TT, FT, PT, UT, GT, ϕT, VT, SPT, ST, IT} <: Integrator{DT, TT}
    equation::IDAE{DT,TT,FT,PT,UT,GT,ϕT}
    tableau::TableauVPARK{TT}
    Δt::TT

    params::SPT
    solver::ST
    iguess::InitialGuessIODE{DT, TT, VT, FT, IT}

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

function IntegratorVPARK{DT,TT,FT,PT,UT,GT,ϕT,VT}(equation::IDAE{DT,TT,FT,PT,UT,GT,ϕT,VT},
                                               tableau::TableauVPARK{TT}, Δt::TT;
                                               nonlinear_solver=DEFAULT_NonlinearSolver,
                                               nmax=DEFAULT_nmax, atol=DEFAULT_atol, rtol=DEFAULT_rtol, stol=DEFAULT_stol,
                                               interpolation=HermiteInterpolation{DT})
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

    # create solution vector for internal stages / nonlinear solver
    x = zeros(DT, N)

    # create params
    params = NonlinearFunctionParametersVPARK{DT,TT,FT,PT,UT,GT,ϕT}(
                                                equation.f, equation.p, equation.u, equation.g, equation.ϕ,
                                                Δt, D, S, R,
                                                tableau.q, tableau.p, tableau.q̃, tableau.p̃, tableau.λ, d_v)

    # create rhs function for nonlinear solver
    function_stages = (x,b) -> function_stages!(x, b, params)

    # create solver
    solver = nonlinear_solver(x, function_stages; nmax=nmax, atol=atol, rtol=rtol, stol=stol, autodiff=false)

    # create initial guess
    iguess = InitialGuessIODE(interpolation, equation, Δt)

    # create integrator
    IntegratorVPARK{DT, TT, FT, PT, UT, GT, ϕT, VT, typeof(params), typeof(solver), typeof(iguess.int)}(
                                        equation, tableau, Δt, params, solver, iguess,
                                        params.q, params.v, params.p, params.λ,
                                        params.y, params.z,
                                        params.Qi, params.Pi, params.Vi, params.Fi,
                                        params.Λp, params.Up, params.Gp)
end


"Integrate DAE with variational partitioned additive Runge-Kutta integrator."
function integrate!{DT,TT,FT,PT,UT,GT,ϕT,VT,N}(int::IntegratorVPARK{DT,TT,FT,PT,UT,GT,ϕT,VT}, sol::SolutionPDAE{DT,TT,N}, m1::Int, m2::Int)
    @assert m1 ≥ 1
    @assert m2 ≤ sol.ni

    # loop over initial conditions
    for m in m1:m2
        # copy initial conditions from solution
        get_initial_conditions!(sol, int.q, int.p, int.λ, m)

        # initialise initial guess
        initialize!(int.iguess, sol.t[0], int.q, int.p)

        for n in 1:sol.ntime
            # set time for nonlinear solver
            int.params.t = sol.t[n]

            # copy previous solution to initial guess
            update!(int.iguess, sol.t[n], int.q, int.p)

            # compute initial guess
            for i in 1:int.tableau.s
                evaluate!(int.iguess, int.y, int.z, int.v, int.tableau.q.c[i], int.tableau.p.c[i])
                for k in 1:int.equation.d
                    int.solver.x[2*(int.equation.d*(i-1)+k-1)+1] = (int.y[k] - int.q[k])/(int.Δt)
                    int.solver.x[2*(int.equation.d*(i-1)+k-1)+2] = int.v[k]
                end
            end

            for i in 1:int.tableau.r
                evaluate!(int.iguess, int.y, int.z, int.v, int.tableau.q̃.c[i], int.tableau.p̃.c[i])
                for k in 1:int.equation.d
                    int.solver.x[2*int.equation.d*int.tableau.s+3*(int.equation.d*(i-1)+k-1)+1] = (int.y[k] - int.q[k])/(int.Δt)
                    int.solver.x[2*int.equation.d*int.tableau.s+3*(int.equation.d*(i-1)+k-1)+2] = (int.z[k] - int.p[k])/(int.Δt)
                    int.solver.x[2*int.equation.d*int.tableau.s+3*(int.equation.d*(i-1)+k-1)+3] = 0
                end
            end

            if isdefined(int.tableau, :d)
                for k in 1:int.equation.d
                    int.solver.x[2*int.equation.d*int.tableau.s+3*int.equation.d*int.tableau.r+k] = 0
                end
            end

            # call nonlinear solver
            solve!(int.solver)

            if !solverStatusOK(int.solver.status, int.solver.params)
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
    end
    nothing
end
