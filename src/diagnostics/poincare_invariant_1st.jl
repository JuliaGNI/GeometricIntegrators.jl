
immutable PoincareInvariant1st{ET,DT,TT,ΘT}
    equ::ET
    Θ::ΘT
    Δt::TT
    nloop::Int
    ntime::Int
    nsave::Int
    nplot::Int
    odir::String
    atol::DT
    rtol::DT
    nmax::Int
end

function PoincareInvariant1st{DT,TT,ΘT}(equ::Equation{DT,TT}, Θ::ΘT, Δt::TT, nloop::Int, ntime::Int, nsave::Int, nplot::Int, odir::String;
                                        atol::DT=2*eps(), rtol::DT=2*eps(), nmax::Int=100)

    println()
    println("First Euler-Poincaré Integral Invariant")
    println("=======================================")
    println()
    println(" nloop = ", nloop)
    println(" ntime = ", ntime)
    println(" Δt    = ", Δt)
    println()

    if !isdir(odir)
        mkpath(odir)
    end

    PoincareInvariant1st{typeof(equ),DT,TT,ΘT}(equ, Θ, Δt, nloop, ntime, nsave, nplot, odir, atol, rtol, nmax)
end



function evaluate_poincare_invariant(pinv::PoincareInvariant1st, integrator, tableau, runid)

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
    p = zeros(sol.q.d)
    v = zeros(size(sol.q.d,1), size(sol.q.d,3))

    compute_one_form(sol.t.t, sol.q.d, p, pinv.Θ)

    for i in 1:size(sol.q.d,2)
        compute_velocity(sol.q.d[:,i,:], v)
        I[i] = compute_loop_integral(p[:,i,:], v)
        J[i] = compute_loop_integral(sol.p.d[:,i,:], v)
    end

    plot_integral_error(sol.t.t, I, pinv.odir * "/" * runid * "_poincare_1st_q.png";
                        plot_title=L"$\Delta \int_\gamma \theta_i (q) \, dq^i$")
    plot_integral_error(sol.t.t, J, pinv.odir * "/" * runid * "_poincare_1st_p.png";
                        plot_title=L"$\Delta \int_\gamma p_i \, dq^i$")
    plot_loop(sol, pinv.nplot, pinv.odir * "/" * runid * "_loop.png")

end



function compute_one_form(t, x, p, Θ)
    @assert size(x) == size(p)
    @assert length(t) == size(x,2) == size(p,2)

    tp = zeros(eltype(x), size(x,1))

    for k in 1:size(x,3)
        for j in 1:size(x,2)
            Θ(t[j], x[:,j,k], tp)
            p[:,j,k] .= tp
        end
    end
end


function compute_velocity(γ, γ̇)
    @assert size(γ) == size(γ̇)

    n = size(γ,2)
    k = 1im*2π*collect(0:div(n-1,2)+1)

    for l in 1:size(γ,1)
        γ̂ = rfft(γ[l,:])
        γ̇[l,:] .= irfft(k .* γ̂, n) / n
    end
end


function compute_loop_integral{T}(p::Array{T,2}, v::Array{T,2})
    local result = zero(T)
    local error  = zero(T)

    for i in 1:size(v,2)
        for l in 1:size(v,1)
            result, error = compensated_summation(p[l,i] * v[l,i], result, error)
        end
    end

    return result
end
