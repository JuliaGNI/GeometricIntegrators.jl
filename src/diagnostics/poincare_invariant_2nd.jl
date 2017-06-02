
using ApproxFun
using StaticArrays


immutable PoincareInvariant2nd{ET,DT,TT,ΩT}
    equ::ET
    ω::ΩT
    Δt::TT
    nx::Int
    ny::Int
    ntime::Int
    nsave::Int
    nplot::Int
    odir::String
    atol::DT
    rtol::DT
    nmax::Int
end

function PoincareInvariant2nd{DT,TT,ΩT}(f_equ::Function, f_surface::Function, ω::ΩT, Δt::TT, d::Int, nx::Int, ny::Int, ntime::Int, nsave::Int, nplot::Int, odir::String;
                                        atol::DT=2*eps(), rtol::DT=2*eps(), nmax::Int=100)

    println()
    println("Second Euler-Poincaré Integral Invariant")
    println("========================================")
    println()
    println(" nx    = ", nx)
    println(" ny    = ", ny)
    println(" ntime = ", ntime)
    println(" Δt    = ", Δt)
    println()

    if !isdir(odir)
        mkpath(odir)
    end

    # compute Chebyshev points
    c = points(Chebyshev(0..1)^2, nx*ny)

    # compute initial conditions
    q₀ = zeros(d, length(c))

    for i in 1:length(c)
        q₀[:,i] .= f_surface(c[i][1], c[i][2])
    end

    # initialise euation
    equ = f_equ(q₀)

    # initialise Poincare invariant
    PoincareInvariant2nd{typeof(equ),DT,TT,ΩT,}(equ, ω, Δt, nx, ny, ntime, nsave, nplot, odir, atol, rtol, nmax)
end



function evaluate_poincare_invariant{ET,DT,TT}(pinv::PoincareInvariant2nd{ET,DT,TT}, integrator, tableau, runid)

    int = integrator(pinv.equ, tableau, pinv.Δt; atol=pinv.atol, rtol=pinv.rtol, nmax=pinv.nmax)
    sol = Solution(pinv.equ, pinv.Δt, pinv.ntime, pinv.nsave)

    println("Running ", tableau.name, " (", runid, ")...")

    integrate!(int, sol)
    # try
    #     integrate!(int, sol)
    # catch DomainError
    #     println("DOMAIN ERROR")
    # end

    println("   Computing surface integrals...")

    local SC  = Chebyshev(0..1)^2
    local SU  = Ultraspherical(1, 0..1)^2
    local SCV = ApproxFun.ArraySpace(SC, sol.nd)
    local SUV = ApproxFun.ArraySpace(SU, sol.nd)
    local Dx  = Derivative(SC, [1,0])
    local Dy  = Derivative(SC, [0,1])
    local Q   = DefiniteIntegral(SU[1]) ⊗ DefiniteIntegral(SU[2])
    # local Q   = DefiniteIntegral(SU)

    local CSU = Conversion(SC, SU)
    local CXU = Conversion(Ultraspherical(1, 0..1) ⊗ Chebyshev(0..1), SU)
    local CYU = Conversion(Chebyshev(0..1) ⊗ Ultraspherical(1, 0..1), SU)

    fq = Fun(Derivative(SC, [1,0]) * Fun(SC, ApproxFun.transform(SC, pinv.equ.q₀[1,:])), SU)

    local nc = ncoefficients(fq)
    local nv = length(values(fq))

    local I = zeros(DT, sol.t.n+1)
    local J = zeros(DT, sol.t.n+1)

    local Iᵢⱼ = zeros(DT, nv)
    local Jᵢⱼ = zeros(DT, nv)

    local qᵢⱼ = zeros(DT, sol.nd, nv)
    local vᵢ  = zeros(DT, sol.nd, nv)
    local vⱼ  = zeros(DT, sol.nd, nv)
    local aᵢ  = zeros(DT, sol.nd, nv)
    local aⱼ  = zeros(DT, sol.nd, nv)

    local ω::Matrix{DT} = zeros(DT, sol.nd, sol.nd)

    local q = permutedims(sol.q.d, [3,1,2])
    local p = permutedims(sol.p.d, [3,1,2])

    for i in 1:size(sol.q.d,2)
        println("      it = ", i-1)

        # create funs for q and p
        fq = Fun(SCV, ApproxFun.transform(SCV, q[:,:,i]))
        fp = Fun(SCV, ApproxFun.transform(SCV, p[:,:,i]))

        for k in 1:sol.nd
            # compute derivatives of q and p and store values
            vᵢ[k,:] .= values(CXU * (Dx * fq[k]))
            vⱼ[k,:] .= values(CYU * (Dy * fq[k]))
            aᵢ[k,:] .= values(CXU * (Dx * fp[k]))
            aⱼ[k,:] .= values(CYU * (Dy * fp[k]))

            # obtain values of q at the same points as the derivatives
            qᵢⱼ[k,:] .= itransform(SU, resize!((CSU * fq[k]).coefficients, nc))
        end

        # compute integrands of integral invariants at all points
        for j in 1:nv
            pinv.ω(sol.t.t[i], qᵢⱼ[:,j], ω)
            Iᵢⱼ[j] = vector_matrix_vector_product(vⱼ[:,j], ω, vᵢ[:,j])
            Jᵢⱼ[j] = dot(aᵢ[:,j], vⱼ[:,j]) - dot(vᵢ[:,j], aⱼ[:,j])
        end

        # compute noncanonical integral invariant
        I[i] = Q * Fun(SU, ApproxFun.transform(SU, Iᵢⱼ))

        # compute canonical integral invariant
        J[i] = Q * Fun(SU, ApproxFun.transform(SU, Jᵢⱼ))
    end


    println("   Plotting results...")

    plot_integral_error(sol.t.t, I, pinv.odir * "/" * runid * "_poincare_2nd_q.png";
                        plot_title=L"$\Delta \int_S \omega_{ij} (q) \, dq^i \wedge dq^j$")
    plot_integral_error(sol.t.t, J, pinv.odir * "/" * runid * "_poincare_2nd_p.png";
                        plot_title=L"$\Delta \int_S dp_i \, dq^i$")
    plot_surface(sol, pinv.nplot, pinv.odir * "/" * runid * "_area.png")

    return I, J
end
