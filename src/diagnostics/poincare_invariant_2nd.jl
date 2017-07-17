
using ApproxFun
using StaticArrays


struct PoincareInvariant2nd{ET,DT,TT,ΩT}
    equ::ET
    ω::ΩT
    Δt::TT
    nx::Int
    ny::Int
    ntime::Int
    nsave::Int
    nt::Int
    I::Vector{DT}
    J::Vector{DT}
    L::Vector{DT}
end

function PoincareInvariant2nd(f_equ::Function, f_surface::Function, ω::ΩT, Δt::TT, d::Int, nx::Int, ny::Int, ntime::Int, nsave::Int, DT=Float64) where {TT,ΩT}

    println()
    println("Second Euler-Poincaré Integral Invariant")
    println("========================================")
    println()
    println(" nx    = ", nx)
    println(" ny    = ", ny)
    println(" ntime = ", ntime)
    println(" nsave = ", nsave)
    println(" Δt    = ", Δt)
    println()

    # compute Chebyshev points
    c = points(Chebyshev(0..1)^2, nx*ny)

    # compute initial conditions
    q₀ = zeros(DT, (d, length(c)))

    for i in 1:length(c)
        q₀[:,i] .= f_surface(c[i][1], c[i][2])
    end

    # initialise euation
    equ = f_equ(q₀)

    # create arrays for results
    nt = div(ntime, nsave)

    I = zeros(DT, nt+1)
    J = zeros(DT, nt+1)
    L = zeros(DT, nt+1)

    # initialise Poincare invariant
    PoincareInvariant2nd{typeof(equ),DT,TT,ΩT}(equ, ω, Δt, nx, ny, ntime, nsave, nt, I, J, L)
end


function evaluate_poincare_invariant(pinv::PoincareInvariant2nd{ET,DT,TT}, sol::Solution) where {ET,DT,TT}

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

    local Iᵢⱼ = zeros(DT, nv)
    local Jᵢⱼ = zeros(DT, nv)
    local Kᵢⱼ = zeros(DT, nv)
    local Lᵢⱼ = zeros(DT, nv)

    local qᵢⱼ = zeros(DT, sol.nd, nv)
    local vᵢ  = zeros(DT, sol.nd, nv)
    local vⱼ  = zeros(DT, sol.nd, nv)
    local aᵢ  = zeros(DT, sol.nd, nv)
    local aⱼ  = zeros(DT, sol.nd, nv)
    local γᵢ  = zeros(DT, sol.nd, nv)
    local γⱼ  = zeros(DT, sol.nd, nv)

    local ω::Matrix{DT} = zeros(DT, sol.nd, sol.nd)

    local q = permutedims(sol.q.d, [3,1,2])
    if isdefined(sol, :p)
        local p = permutedims(sol.p.d, [3,1,2])
    end
    if isdefined(sol, :λ)
        local λ = permutedims(sol.λ.d, [3,1,2])
    end

    for i in 1:size(sol.q.d,2)
        println("      it = ", i-1)

        # create funs for q and p
        fq = Fun(SCV, ApproxFun.transform(SCV, q[:,:,i]))
        if isdefined(sol, :p)
            fp = Fun(SCV, ApproxFun.transform(SCV, p[:,:,i]))
        end
        if isdefined(sol, :λ)
            fλ = Fun(SCV, ApproxFun.transform(SCV, λ[:,:,i]))
        end

        for k in 1:sol.nd
            # compute derivatives of q and p and store values
            vᵢ[k,:] .= values(CXU * (Dx * fq[k]))
            vⱼ[k,:] .= values(CYU * (Dy * fq[k]))
            if isdefined(sol, :p)
                aᵢ[k,:] .= values(CXU * (Dx * fp[k]))
                aⱼ[k,:] .= values(CYU * (Dy * fp[k]))
            end
            if isdefined(sol, :λ)
                γᵢ[k,:] .= values(CXU * (Dx * fλ[k]))
                γⱼ[k,:] .= values(CYU * (Dy * fλ[k]))
            end

            # obtain values of q at the same points as the derivatives
            qᵢⱼ[k,:] .= itransform(SU, resize!((CSU * fq[k]).coefficients, nc))
        end

        # compute integrands of integral invariants at all points
        for j in 1:nv
            pinv.ω(sol.t.t[i], qᵢⱼ[:,j], ω)
            Iᵢⱼ[j] = vector_matrix_vector_product(vⱼ[:,j], ω, vᵢ[:,j])
            if isdefined(sol, :p)
                Jᵢⱼ[j] = dot(aᵢ[:,j], vⱼ[:,j]) - dot(vᵢ[:,j], aⱼ[:,j])
            end
            if isdefined(sol, :λ)
                Kᵢⱼ[j] = vector_matrix_vector_product(γⱼ[:,j], ω, γᵢ[:,j])
                Lᵢⱼ[j] = Iᵢⱼ[j] - pinv.Δt^2 * Kᵢⱼ[j]
            end
        end

        # compute noncanonical integral invariant
        pinv.I[i] = Q * Fun(SU, ApproxFun.transform(SU, Iᵢⱼ))

        # compute canonical integral invariant
        if isdefined(sol, :p)
            pinv.J[i] = Q * Fun(SU, ApproxFun.transform(SU, Jᵢⱼ))
        end

        # compute corrected integral invariant
        if isdefined(sol, :λ)
            pinv.L[i] = Q * Fun(SU, ApproxFun.transform(SU, Lᵢⱼ))
        end
    end

    return (pinv.I, pinv.J, pinv.L)
end


function CommonFunctions.write_to_hdf5(pinv::PoincareInvariant2nd, sol::Solution, output_file::String)
    # h5open(output_file, isfile(output_file) ? "r+" : "w") do h5
    h5open(output_file, "w") do h5

        write(h5, "t", sol.t.t)
        write(h5, "I", pinv.I)

        isdefined(sol, :p) ? write(h5, "J", pinv.J) : nothing
        isdefined(sol, :λ) ? write(h5, "L", pinv.L) : nothing

    end
end
