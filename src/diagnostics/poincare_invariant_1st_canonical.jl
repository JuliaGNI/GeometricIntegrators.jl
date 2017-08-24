
struct PoincareInvariant1stCanonical{ET,DT,TT}
    equ::ET
    Δt::TT
    nloop::Int
    ntime::Int
    nsave::Int
    nt::Int
    I::Vector{DT}
end

function PoincareInvariant1stCanonical(f_equ::Function, f_loop_q::Function, f_loop_p::Function, Δt::TT, d::Int, nloop::Int, ntime::Int, nsave::Int, DT=Float64) where {TT}

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
    p₀ = zeros(DT, (d, nloop))

    for i in 1:nloop
        q₀[:,i] .= f_loop_q(i/nloop)
        p₀[:,i] .= f_loop_p(i/nloop)
    end

    # initialise euation
    equ = f_equ(q₀, p₀)

    # create arrays for results
    nt = div(ntime, nsave)

    I = zeros(DT, nt+1)

    PoincareInvariant1stCanonical{typeof(equ),DT,TT}(equ, Δt, nloop, ntime, nsave, nt, I)
end


function evaluate_poincare_invariant(pinv::PoincareInvariant1stCanonical, sol::SolutionPODE)
    v = zeros(size(sol.q.d,1), size(sol.q.d,3))

    for i in 1:size(sol.q.d,2)
        compute_velocity(sol.q.d[:,i,:], v)
        pinv.I[i] = compute_loop_integral(sol.p.d[:,i,:], v)
    end

    pinv.I
end


function CommonFunctions.write_to_hdf5(pinv::PoincareInvariant1stCanonical, sol::Solution, output_file::String)
    # h5open(output_file, isfile(output_file) ? "r+" : "w") do h5
    h5open(output_file, "w") do h5
        write(h5, "t", sol.t.t)
        write(h5, "I", pinv.I)
    end
end
