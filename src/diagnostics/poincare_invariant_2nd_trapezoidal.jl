
immutable PoincareInvariant2ndTrapezoidal{ET,DT,TT,ΩT}
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

function PoincareInvariant2ndTrapezoidal{DT,TT,ΩT}(f_equ::Function, f_surface::Function, ω::ΩT, Δt::TT, d::Int, nx::Int, ny::Int, ntime::Int, nsave::Int, nplot::Int, odir::String;
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

    q₀ = zeros(d, nx*ny)

    for j in 1:ny
        for i in 1:nx
            q₀[:,nx*(j-1)+i] = f_surface(i/nx, j/nx)
        end
    end

    equ = f_equ(q₀)

    PoincareInvariant2ndTrapezoidal{typeof(equ),DT,TT,ΩT}(equ, ω, Δt, nx, ny, ntime, nsave, nplot, odir, atol, rtol, nmax)
end



function evaluate_poincare_invariant(pinv::PoincareInvariant2ndTrapezoidal, integrator, tableau, runid)

    int = integrator(pinv.equ, tableau, pinv.Δt; atol=pinv.atol, rtol=pinv.rtol, nmax=pinv.nmax)
    sol = Solution(pinv.equ, pinv.Δt, pinv.ntime, pinv.nsave)

    println("Running ", tableau.name, " (", runid, ")...")

    integrate!(int, sol)
    # try
    #     integrate!(int, sol)
    # catch DomainError
    #     println("DOMAIN ERROR")
    # end

    I = zeros(sol.t.n+1)
    J = zeros(sol.t.n+1)

    for i in 1:size(sol.q.d,2)
        I[i] = surface_integral(sol.t.t[i], sol.q.d[:,i,:], pinv.ω, pinv.nx, pinv.ny)
        J[i] = surface_integral_canonical(sol.q.d[:,i,:], sol.p.d[:,i,:], pinv.nx, pinv.ny)
    end

    plot_integral_error(sol.t.t, I, pinv.odir * "/" * runid * "_poincare_2nd_q.png";
                        plot_title=L"$\Delta \int_S \omega_{ij} (q) \, dq^i \wedge dq^j$")
    plot_integral_error(sol.t.t, J, pinv.odir * "/" * runid * "_poincare_2nd_p.png";
                        plot_title=L"$\Delta \int_S dp_i \, dq^i$")
    plot_surface(sol, pinv.nplot, pinv.odir * "/" * runid * "_area.png")

    return I, J
end


function interpolate_trajectory(x, i1, j1, λ, μ, γ, nx, ny)
    @assert length(γ) == size(x, 1)
    @assert λ ≥ 0.
    @assert λ ≤ 1.
    @assert μ ≥ 0.
    @assert μ ≤ 1.

    @assert i1 > 0
    @assert i1 < nx
    @assert j1 > 0
    @assert j1 < ny

    i2 = i1 + 1
    j2 = j1 + 1

    for k in 1:length(γ)
        γ[k] = x[k, nx*(j1-1)+i1] * (1-λ) * (1-μ) +
               x[k, nx*(j1-1)+i2] *    λ  * (1-μ) +
               x[k, nx*(j2-1)+i1] * (1-λ) *    μ  +
               x[k, nx*(j2-1)+i2] *    λ  *    μ
    end
end


function interpolate_derivative_i(x, i1, j1, λ, μ, γ̇, nx, ny)
    @assert length(γ̇) == size(x, 1)
    @assert λ ≥ 0.
    @assert λ ≤ 1.
    @assert μ ≥ 0.
    @assert μ ≤ 1.

    @assert i1 > 0
    @assert i1 < nx
    @assert j1 > 0
    @assert j1 < ny

    i2 = i1 + 1
    j2 = j1 + 1

    for k in 1:length(γ̇)
        γ̇[k] = (x[k, nx*(j1-1)+i2] - x[k, nx*(j1-1)+i1]) * (1-μ) +
               (x[k, nx*(j2-1)+i2] - x[k, nx*(j2-1)+i1]) *    μ
    end
end


function interpolate_derivative_j(x, i1, j1, λ, μ, γ̇, nx, ny)
    @assert length(γ̇) == size(x, 1)
    @assert λ ≥ 0.
    @assert λ ≤ 1.
    @assert μ ≥ 0.
    @assert μ ≤ 1.

    @assert i1 > 0
    @assert i1 < nx
    @assert j1 > 0
    @assert j1 < ny

    i2 = i1 + 1
    j2 = j1 + 1

    for k in 1:length(γ̇)
        γ̇[k] = (x[k, nx*(j2-1)+i1] - x[k, nx*(j1-1)+i1]) * (1-λ) +
               (x[k, nx*(j2-1)+i2] - x[k, nx*(j1-1)+i2]) *    λ
    end
end


function integrate{DT,TT}(t, γ, γ̇ᵢ, γ̇ⱼ, ω, b::Vector{TT}, c::Vector{TT}, q::Vector{DT}, vᵢ::Vector{DT}, vⱼ::Vector{DT}, B::Matrix{DT})
    @assert length(b) == length(c)

    local result = zero(DT)

    for i in 1:length(b)
        for j in 1:length(b)
            γ(c[i], c[j], q)
            γ̇ᵢ(c[i], c[j], vᵢ)
            γ̇ⱼ(c[i], c[j], vⱼ)
            ω(t, q, B)

            result += b[i] * b[j] * vector_matrix_vector_product(vⱼ, B, vᵢ)
        end
    end

    return result
end

function surface_integral{DT}(t, x::Matrix{DT}, ω, nx, ny)
    const b = [0.5, 0.5]
    const c = [0.0, 1.0]

    local q  = zeros(DT, size(x,1))
    local vᵢ = zeros(DT, size(x,1))
    local vⱼ = zeros(DT, size(x,1))
    local B  = zeros(DT, size(x,1), size(x,1))
    local I  = zero(DT)

    integrate_trapezoidal = (γ, γ̇ᵢ, γ̇ⱼ) -> integrate(t, γ, γ̇ᵢ, γ̇ⱼ, ω, b, c, q, vᵢ, vⱼ, B)

    for j in 1:ny-1
        for i in 1:nx-1
            γ  = (λ, μ, y) -> interpolate_trajectory(x, i, j, λ, μ, y, nx, ny)
            γ̇ᵢ = (λ, μ, y) -> interpolate_derivative_i(x, i, j, λ, μ, y, nx, ny)
            γ̇ⱼ = (λ, μ, y) -> interpolate_derivative_j(x, i, j, λ, μ, y, nx, ny)
            I += integrate_trapezoidal(γ, γ̇ᵢ, γ̇ⱼ)
        end
    end

    return I
end



function integrate_canonical{DT,TT}(γ̇ᵢ, γ̇ⱼ, Θ̇ᵢ, Θ̇ⱼ, b::Vector{TT}, c::Vector{TT}, vᵢ::Vector{DT}, vⱼ::Vector{DT})
    @assert length(b) == length(c)

    local result = zero(DT)

    for i in 1:length(b)
        for j in 1:length(b)
            Θ̇ᵢ(c[i], c[j], vᵢ)
            γ̇ⱼ(c[i], c[j], vⱼ)
            result += b[i] * b[j] * dot(vᵢ,vⱼ)

            γ̇ᵢ(c[i], c[j], vᵢ)
            Θ̇ⱼ(c[i], c[j], vⱼ)
            result -= b[i] * b[j] * dot(vᵢ,vⱼ)
        end
    end

    return result
end

function surface_integral_canonical{DT}(q::Matrix{DT}, p::Matrix{DT}, nx, ny)
    const b = [0.5, 0.5]
    const c = [0.0, 1.0]

    local vᵢ = zeros(DT, size(q,1))
    local vⱼ = zeros(DT, size(q,1))
    local I  = zero(DT)

    integrate_trapezoidal = (γ, γ̇ᵢ, γ̇ⱼ, Θ̇ᵢ, Θ̇ⱼ) -> integrate_canonical(γ̇ᵢ, γ̇ⱼ, Θ̇ᵢ, Θ̇ⱼ, b, c, vᵢ, vⱼ)

    for j in 1:ny-1
        for i in 1:nx-1
            γ̇ᵢ = (λ, μ, y) -> interpolate_derivative_i(q, i, j, λ, μ, y, nx, ny)
            γ̇ⱼ = (λ, μ, y) -> interpolate_derivative_j(q, i, j, λ, μ, y, nx, ny)
            Θ̇ᵢ = (λ, μ, y) -> interpolate_derivative_i(p, i, j, λ, μ, y, nx, ny)
            Θ̇ⱼ = (λ, μ, y) -> interpolate_derivative_j(p, i, j, λ, μ, y, nx, ny)
            I += integrate_trapezoidal(γ, γ̇ᵢ, γ̇ⱼ, Θ̇ᵢ, Θ̇ⱼ)
        end
    end

    return I
end
