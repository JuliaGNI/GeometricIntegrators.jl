
struct PoincareInvariant1st{ET,DT,TT,ΘT}
    equ::ET
    Θ::ΘT
    Δt::TT
    nloop::Int
    ntime::Int
    nsave::Int
    nt::Int
    I::Vector{DT}
    J::Vector{DT}
    L::Vector{DT}
end

function PoincareInvariant1st(f_equ::Function, f_loop::Function, Θ::ΘT, Δt::TT, d::Int, nloop::Int, ntime::Int, nsave::Int=1, DT=Float64) where {TT,ΘT}

    if get_config(:verbosity) > 1
        println()
        println("First Euler-Poincaré Integral Invariant")
        println("=======================================")
        println()
        println(" nloop = ", nloop)
        println(" ntime = ", ntime)
        println(" nsave = ", nsave)
        println(" Δt    = ", Δt)
        println()
    end

    # compute initial conditions
    q₀ = zeros(DT, (d, nloop))

    for i in 1:nloop
        q₀[:,i] .= f_loop(i/nloop)
    end

    # initialise euation
    equ = f_equ(q₀)

    # create arrays for results
    nt = div(ntime, nsave)

    I = zeros(DT, nt+1)
    J = zeros(DT, nt+1)
    L = zeros(DT, nt+1)

    PoincareInvariant1st{typeof(equ),DT,TT,ΘT}(equ, Θ, Δt, nloop, ntime, nsave, nt, I, J, L)
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
