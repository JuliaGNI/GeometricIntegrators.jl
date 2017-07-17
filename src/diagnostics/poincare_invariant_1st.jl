
struct PoincareInvariant1st{ET,DT,TT,ΘT}
    equ::ET
    Θ::ΘT
    Δt::TT
    nloop::Int
    ntime::Int
    nsave::Int
    nplot::Int
    nt::Int
    I::Vector{DT}
    J::Vector{DT}
    L::Vector{DT}
end

function PoincareInvariant1st(equ::Equation{DT,TT}, Θ::ΘT, Δt::TT, nloop::Int, ntime::Int, nsave::Int, nplot::Int) where {DT,TT,ΘT}

    println()
    println("First Euler-Poincaré Integral Invariant")
    println("=======================================")
    println()
    println(" nloop = ", nloop)
    println(" ntime = ", ntime)
    println(" nsave = ", nsave)
    println(" Δt    = ", Δt)
    println()

    nt = div(ntime, nsave)

    I = zeros(DT, nt+1)
    J = zeros(DT, nt+1)
    L = zeros(DT, nt+1)

    PoincareInvariant1st{typeof(equ),DT,TT,ΘT}(equ, Θ, Δt, nloop, ntime, nsave, nplot, nt, I, J, L)
end


function evaluate_poincare_invariant(pinv::PoincareInvariant1st, sol::Solution)
    p = zeros(sol.q.d)
    g = zeros(sol.q.d)
    v = zeros(size(sol.q.d,1), size(sol.q.d,3))
    γ = zeros(size(sol.q.d,1), size(sol.q.d,3))

    compute_one_form(sol.t.t, sol.q.d, p, pinv.Θ)

    if isdefined(sol, :λ)
        compute_correction(sol.t.t, sol.q.d, sol.λ.d, g, pinv.equ.g)
    end

    for i in 1:size(sol.q.d,2)
        compute_velocity(sol.q.d[:,i,:], v)
        pinv.I[i] = compute_loop_integral(p[:,i,:], v)

        if isdefined(sol, :p)
            pinv.J[i] = compute_loop_integral(sol.p.d[:,i,:], v)
        end

        if isdefined(sol, :λ)
            compute_velocity(sol.λ.d[:,i,:], γ)
            pinv.L[i] = compute_loop_integral(p[:,i,:] .- pinv.Δt .* g[:,i,:], v .- pinv.Δt .* γ)
        end
    end

    (pinv.I, pinv.J, pinv.L)
end


function CommonFunctions.write_to_hdf5(pinv::PoincareInvariant1st, sol::Solution, output_file::String)
    # h5open(output_file, isfile(output_file) ? "r+" : "w") do h5
    h5open(output_file, "w") do h5

        write(h5, "t", sol.t.t)
        write(h5, "I", pinv.I)

        isdefined(sol, :p) ? write(h5, "J", pinv.J) : nothing
        isdefined(sol, :λ) ? write(h5, "L", pinv.L) : nothing

    end
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


function compute_correction(t, x, λ, σ, g)
    @assert size(x) == size(λ) == size(σ)
    @assert length(t) == size(x,2) == size(σ,2)

    tσ = zeros(eltype(x), size(x,1))

    for k in 1:size(x,3)
        for j in 1:size(x,2)
            g(t[j], x[:,j,k], λ[:,j,k], tσ)
            σ[:,j,k] .= tσ
        end
    end
end


function compute_loop_integral(p::Array{T,2}, v::Array{T,2}) where {T}
    local result = zero(T)
    local error  = zero(T)

    for i in 1:size(v,2)
        for l in 1:size(v,1)
            result, error = compensated_summation(p[l,i] * v[l,i], result, error)
        end
    end

    return result
end
