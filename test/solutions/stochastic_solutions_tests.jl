
using GeometricIntegrators.CommonFunctions
using GeometricIntegrators.Solutions
using GeometricIntegrators.Solutions: createHDF5
using GeometricProblems.KuboOscillator
using Test


nd = 3
ns = 5
nt = 10
Δt = .1
dt = Float64

t0    = 0.
x0    = rand(2)
y0    = rand(2)
q0    = rand(1)
p0    = rand(1)

t1    = 1.
x1    = rand(2, ns)
y1    = rand(2, ns)
q1    = rand(1, ns)
p1    = rand(1, ns)

t2    = t1 + (t1-t0)
x2    = rand(2, ns)
y2    = rand(2, ns)
q2    = rand(1, ns)
p2    = rand(1, ns)

tx    = zero(x0)
ty    = zero(y0)
tq    = zero(q0)
tp    = zero(p0)

xs    = rand(2,nt)
ys    = rand(2,nt)
qs    = rand(1,nt)
ps    = rand(1,nt)

Xs    = rand(2,nt,ns)
Ys    = rand(2,nt,ns)
Qs    = rand(1,nt,ns)
Ps    = rand(1,nt,ns)

sde  = kubo_oscillator_sde_2()
psde = kubo_oscillator_psde_2()

x100  = rand(ndims(sde),  100)
q100  = rand(ndims(psde), 100)
p100  = rand(ndims(psde), 100)

h5file = "test.hdf5"


@testset "$(rpad("Wiener Process",80))" begin

    wp = WienerProcess(dt, 1, nt, 1,  Δt, :strong)
    @test wp == WienerProcess(Δt, wp.ΔW[:], wp.ΔZ[:], :strong)
    @test ndims(wp) == 2

    wp = WienerProcess(dt, 1, nt, 1,  Δt, :weak)
    @test wp == WienerProcess(Δt, wp.ΔW[:], wp.ΔZ[:], :weak)
    @test ndims(wp) == 2

    wp = WienerProcess(dt, nd, nt, 1,  Δt, :strong)
    @test wp == WienerProcess(Δt, wp.ΔW, wp.ΔZ, :strong)
    @test ndims(wp) == 2

    wp = WienerProcess(dt, nd, nt, 1,  Δt, :weak)
    @test wp == WienerProcess(Δt, wp.ΔW, wp.ΔZ, :weak)
    @test ndims(wp) == 2

    wp = WienerProcess(dt, nd, nt, ns, Δt, :strong)
    @test wp == WienerProcess(Δt, wp.ΔW, wp.ΔZ, :strong)
    @test ndims(wp) == 3

    wp = WienerProcess(dt, nd, nt, ns, Δt, :weak)
    @test wp == WienerProcess(Δt, wp.ΔW, wp.ΔZ, :weak)
    @test ndims(wp) == 3

end


@testset "$(rpad("SDE Solution",80))" begin
    asol = AtomicSolution(sde)

    # test constructors and general functionality
    sol = Solution(kubo_oscillator_sde_1(), Δt, nt)
    @test typeof(sol) <: SolutionSDE
    @test sol.nd == 2
    @test sol.nm == 1
    @test sol.ns == 1
    @test sol.nt == nt

    sol = Solution(kubo_oscillator_sde_2(), Δt, nt)
    @test typeof(sol) <: SolutionSDE
    @test sol.nd == 2
    @test sol.nm == 1
    @test sol.ns == 3
    @test sol.nt == nt

    sol = Solution(kubo_oscillator_sde_3(), Δt, nt)
    @test typeof(sol) <: SolutionSDE
    @test sol.nd == 2
    @test sol.nm == 1
    @test sol.ns == 3
    @test sol.nt == nt

    sol  = Solution(sde, Δt, nt)
    sol0 = Solution(similar(sde, x0), Δt, nt)
    sol1 = Solution(similar(sde, x1), Δt, nt)

    @test typeof(sol0) <: SolutionSDE
    @test typeof(sol1) <: SolutionSDE

    @test sol != sol0
    @test sol != sol1

    @test sol == SSolutionSDE(sde, Δt, sol.W.ΔW, sol.W.ΔZ, sol.ntime)

    @test nsave(sol) == 1
    @test ntime(sol) == nt
    @test timesteps(sol) == Δt .* collect(0:nt)

    # test initial conditions
    set_initial_conditions!(sol, t0, x0)
    get_initial_conditions!(sol, tx, 1)
    get_initial_conditions!(sol, asol, 1)

    @test tx == x0
    @test asol.t == t0
    @test asol.q == x0

    δt, δx = get_initial_conditions(sol, 1)
    @test δt == t0
    @test δx == x0

    set_initial_conditions!(sol1, similar(sde, t1, x2))
    for i in 1:ns
        get_initial_conditions!(sol1, tx, i)
        @test tx == x2[:,i]
    end


    # test set/get solution
    sol1 = Solution(similar(kubo_oscillator_sde_3(), x0), Δt, nt)
    sol2 = Solution(similar(kubo_oscillator_sde_3(), x0), Δt, nt)
    for i in 1:nt
        tx .= xs[:,i]
        asol.q .= xs[:,i]
        set_solution!(sol1, tx, i)
        set_solution!(sol2, asol, i)
    end
    @test sol1.q[:,1:nt] == xs
    @test sol2.q[:,1:nt] == xs

    sol1 = Solution(similar(sde, x1), Δt, nt)
    sol2 = Solution(similar(sde, x1), Δt, nt)
    for i in 1:nt
        for k in 1:ns
            tx .= Xs[:,i,k]
            asol.q .= Xs[:,i,k]
            set_solution!(sol1, tx, i, k)
            set_solution!(sol2, asol, i, k)
        end
    end
    @test sol1.q[:,1:nt,:] == Xs
    @test sol2.q[:,1:nt,:] == Xs

    # test hdf5 in- and output
    sol1 = Solution(kubo_oscillator_sde_3(x0), Δt, nt)
    create_hdf5!(sol1, h5file)
    write_to_hdf5(sol1)
    close(sol1)
    @test isfile(h5file)
    sol2 = SSolutionSDE(h5file)
    @test sol1   != sol2
    @test sol1.t == sol2.t
    @test sol1.q == sol2.q
    @test sol1.ntime == sol2.ntime
    @test sol1.nsave == sol2.nsave
    @test sol1.nd == sol2.nd
    @test sol1.nm == sol2.nm
    @test sol1.ns == sol2.ns
    @test sol1.K  == sol2.K
    @test conv(sol1) == conv(sol2)
    rm(h5file)

    sol1 = Solution(kubo_oscillator_sde_3(x1), Δt, nt)
    create_hdf5!(sol1, h5file)
    write_to_hdf5(sol1)
    close(sol1)
    @test isfile(h5file)
    sol2 = SSolutionSDE(h5file)
    @test sol1   != sol2
    @test sol1.t == sol2.t
    @test sol1.q == sol2.q
    @test sol1.ntime == sol2.ntime
    @test sol1.nsave == sol2.nsave
    @test sol1.nd == sol2.nd
    @test sol1.nm == sol2.nm
    @test sol1.ns == sol2.ns
    @test sol1.K  == sol2.K
    @test conv(sol1) == conv(sol2)
    rm(h5file)

    sde1 = kubo_oscillator_sde_3(x1)
    sol1 = Solution(sde1, Δt, nt)
    create_hdf5!(sol1, h5file; save_W=false)
    write_to_hdf5(sol1)
    close(sol1)
    @test isfile(h5file)
    sol2 = SSolutionSDE(h5file)
    @test sol2.W.ΔW == zeros(eltype(q1), 0, 0, 0)
    @test sol2.W.ΔZ == zeros(eltype(q1), 0, 0, 0)
    rm(h5file)

    # test nsave and nwrite parameters
    sol = Solution(sde, Δt, 20, nsave=2)
    @test sol.nt == 10

    sol = Solution(sde, Δt, 20, nsave=2, nwrite=10)
    @test sol.nt == 5

    # test reset
    sol = Solution(sde, Δt, nt)
    reset!(sol)
    @test sol.t[0]   == t1
    @test sol.t[end] == t2
    @test offset(sol) == nt

    sde1 = kubo_oscillator_sde_1()
    sol1 = Solution(sde1, Δt, 100, nwrite=10)
    create_hdf5!(sol1, h5file)
    ΔW = zeros(sde1.m, 100)
    ΔZ = zeros(sde1.m, 100)
    for i in 1:10
        ΔW[:,(i-1)*10+1:i*10] .= sol1.W.ΔW
        ΔZ[:,(i-1)*10+1:i*10] .= sol1.W.ΔZ
        for j in 1:10
            set_solution!(sol1, x100[:,(i-1)*10+j], j)
        end
        write_to_hdf5(sol1)
        reset!(sol1)
    end
    close(sol1)
    sol2 = SSolutionSDE(h5file)
    @test sol2.t.t ≈ Δt .* collect(0:100) atol=1E-14
    @test sol2.q.d == hcat(reshape(sde1.q₀, (ndims(sde1),1)), x100)
    @test sol2.W.ΔW == ΔW
    @test sol2.W.ΔZ == ΔZ
    @test sol2.ntime == 100
    @test sol2.nsave == 1
    rm(h5file)

    sol1 = Solution(sde1, Δt, 100, nsave=2, nwrite=10)
    create_hdf5!(sol1, h5file)
    ΔW = zeros(sde1.m, 100)
    ΔZ = zeros(sde1.m, 100)
    for i in 1:10
        ΔW[:,(i-1)*10+1:i*10] .= sol1.W.ΔW
        ΔZ[:,(i-1)*10+1:i*10] .= sol1.W.ΔZ
        for j in 1:10
            set_solution!(sol1, x100[:,(i-1)*10+j], j)
        end
        write_to_hdf5(sol1)
        reset!(sol1)
    end
    close(sol1)
    sol2 = SSolutionSDE(h5file)
    @test sol2.t.t ≈ Δt .* collect(0:2:100) atol=1E-14
    @test sol2.q.d == hcat(reshape(sde1.q₀, (ndims(sde1),1)), x100)[:,1:2:end]
    @test sol2.W.ΔW == ΔW
    @test sol2.W.ΔZ == ΔZ
    @test sol2.ntime == 100
    @test sol2.nsave == 2
    rm(h5file)
end


@testset "$(rpad("PSDE Solution",80))" begin
    asol = AtomicSolution(psde)

    # test constructors and general functionality
    sol = Solution(kubo_oscillator_psde_1(), Δt, nt)
    @test typeof(sol) <: SolutionPSDE
    @test sol.nd == 1
    @test sol.nm == 1
    @test sol.ns == 1
    @test sol.nt == nt

    sol = Solution(kubo_oscillator_psde_2(), Δt, nt)
    @test typeof(sol) <: SolutionPSDE
    @test sol.nd == 1
    @test sol.nm == 1
    @test sol.ns == 3
    @test sol.nt == nt

    sol = Solution(kubo_oscillator_psde_3(), Δt, nt)
    @test typeof(sol) <: SolutionPSDE
    @test sol.nd == 1
    @test sol.nm == 1
    @test sol.ns == 3
    @test sol.nt == nt

    sol  = Solution(psde, Δt, nt)
    sol0 = Solution(similar(psde, q0, p0), Δt, nt)
    sol1 = Solution(similar(psde, q1, p1), Δt, nt)

    @test typeof(sol0) <: SolutionPSDE
    @test typeof(sol1) <: SolutionPSDE

    @test sol != sol0
    @test sol != sol1

    @test sol == SSolutionPSDE(psde, Δt, sol.W.ΔW, sol.W.ΔZ, sol.ntime)

    @test nsave(sol) == 1
    @test ntime(sol) == nt
    @test timesteps(sol) == Δt .* collect(0:nt)

    # test initial conditions
    set_initial_conditions!(sol, t0, q0, p0)
    get_initial_conditions!(sol, tq, tp, 1)
    get_initial_conditions!(sol, asol, 1)

    @test tq == q0
    @test tp == p0
    @test asol.t == t0
    @test asol.q == q0
    @test asol.p == p0

    δt, δq, δp = get_initial_conditions(sol, 1)
    @test δt == t0
    @test δq == q0
    @test δp == p0

    set_initial_conditions!(sol1, similar(psde, t1, q2, p2))
    for i in 1:ns
        get_initial_conditions!(sol1, tq, tp, i)
        @test tq == q2[:,i]
        @test tp == p2[:,i]
    end

    # test set/get solution
    sol1 = Solution(similar(kubo_oscillator_psde_3(), q0, p0), Δt, nt)
    sol2 = Solution(similar(kubo_oscillator_psde_3(), q0, p0), Δt, nt)
    for i in 1:nt
        tq .= qs[:,i]
        tp .= ps[:,i]
        asol.q .= qs[:,i]
        asol.p .= ps[:,i]
        set_solution!(sol1, tq, tp, i)
        set_solution!(sol2, asol, i)
    end
    @test sol1.q[:,1:nt] == qs
    @test sol1.p[:,1:nt] == ps
    @test sol2.q[:,1:nt] == qs
    @test sol2.p[:,1:nt] == ps

    sol1 = Solution(similar(kubo_oscillator_psde_3(), q1, p1), Δt, nt)
    sol2 = Solution(similar(kubo_oscillator_psde_3(), q1, p1), Δt, nt)
    for i in 1:nt
        for k in 1:ns
            tq .= Qs[:,i,k]
            tp .= Ps[:,i,k]
            asol.q .= Qs[:,i,k]
            asol.p .= Ps[:,i,k]
            set_solution!(sol1, tq, tp, i, k)
            set_solution!(sol2, asol, i, k)
        end
    end
    @test sol1.q[:,1:nt,:] == Qs
    @test sol1.p[:,1:nt,:] == Ps
    @test sol2.q[:,1:nt,:] == Qs
    @test sol2.p[:,1:nt,:] == Ps

    # test hdf5 in- and output
    sol1 = Solution(kubo_oscillator_psde_3(q0, p0), Δt, nt)
    create_hdf5!(sol1, h5file)
    write_to_hdf5(sol1)
    close(sol1)
    @test isfile(h5file)
    sol2 = SSolutionPSDE(h5file)
    @test sol1   != sol2
    @test sol1.t == sol2.t
    @test sol1.q == sol2.q
    @test sol1.p == sol2.p
    @test sol1.ntime == sol2.ntime
    @test sol1.nsave == sol2.nsave
    @test sol1.nd == sol2.nd
    @test sol1.nm == sol2.nm
    @test sol1.ns == sol2.ns
    @test sol1.K  == sol2.K
    @test conv(sol1) == conv(sol2)
    rm(h5file)

    sol1 = Solution(kubo_oscillator_psde_3(q1, p1), Δt, nt)
    create_hdf5!(sol1, h5file)
    write_to_hdf5(sol1)
    close(sol1)
    @test isfile(h5file)
    sol2 = SSolutionPSDE(h5file)
    @test sol1   != sol2
    @test sol1.t == sol2.t
    @test sol1.q == sol2.q
    @test sol1.p == sol2.p
    @test sol1.ntime == sol2.ntime
    @test sol1.nsave == sol2.nsave
    @test sol1.nd == sol2.nd
    @test sol1.nm == sol2.nm
    @test sol1.ns == sol2.ns
    @test sol1.K  == sol2.K
    @test conv(sol1) == conv(sol2)
    rm(h5file)

    psde1 = kubo_oscillator_psde_3(q1, p1)
    sol1 = Solution(psde1, Δt, nt)
    create_hdf5!(sol1, h5file; save_W=false)
    write_to_hdf5(sol1)
    close(sol1)
    @test isfile(h5file)
    sol2 = SSolutionPSDE(h5file)
    @test sol2.W.ΔW == zeros(eltype(q1), 0, 0, 0)
    @test sol2.W.ΔZ == zeros(eltype(q1), 0, 0, 0)
    rm(h5file)

    # test nsave and nwrite parameters
    sol = Solution(psde, Δt, 20, nsave=2)
    @test sol.nt == 10

    sol = Solution(psde, Δt, 20, nsave=2, nwrite=10)
    @test sol.nt == 5

    # test reset
    sol = Solution(psde, Δt, nt)
    reset!(sol)
    @test sol.t[0]   == t1
    @test sol.t[end] == t2
    @test offset(sol) == nt

    psde1 = kubo_oscillator_psde_1()
    sol1 = Solution(psde1, Δt, 100, nwrite=10)
    create_hdf5!(sol1, h5file)
    ΔW = zeros(psde1.m, 100)
    ΔZ = zeros(psde1.m, 100)
    for i in 1:10
        ΔW[:,(i-1)*10+1:i*10] .= sol1.W.ΔW
        ΔZ[:,(i-1)*10+1:i*10] .= sol1.W.ΔZ
        for j in 1:10
            oj = (i-1)*10+j
            set_solution!(sol1, q100[:,oj], p100[:,oj], j)
        end
        write_to_hdf5(sol1)
        reset!(sol1)
    end
    close(sol1)
    sol2 = SSolutionPSDE(h5file)
    @test sol2.t.t ≈ Δt .* collect(0:100) atol=1E-14
    @test sol2.q.d == hcat(reshape(psde1.q₀, (ndims(psde1),1)), q100)
    @test sol2.p.d == hcat(reshape(psde1.p₀, (ndims(psde1),1)), p100)
    @test sol2.W.ΔW == ΔW
    @test sol2.W.ΔZ == ΔZ
    @test sol2.ntime == 100
    @test sol2.nsave == 1
    rm(h5file)

    sol1 = Solution(psde1, Δt, 100, nsave=2, nwrite=10)
    create_hdf5!(sol1, h5file)
    ΔW = zeros(psde1.m, 100)
    ΔZ = zeros(psde1.m, 100)
    for i in 1:10
        ΔW[:,(i-1)*10+1:i*10] .= sol1.W.ΔW
        ΔZ[:,(i-1)*10+1:i*10] .= sol1.W.ΔZ
        for j in 1:10
            oj = (i-1)*10+j
            set_solution!(sol1, q100[:,oj], p100[:,oj], j)
        end
        write_to_hdf5(sol1)
        reset!(sol1)
    end
    close(sol1)
    sol2 = SSolutionPSDE(h5file)
    @test sol2.t.t ≈ Δt .* collect(0:2:100) atol=1E-14
    @test sol2.q.d == hcat(reshape(psde1.q₀, (ndims(psde1),1)), q100)[:,1:2:end]
    @test sol2.p.d == hcat(reshape(psde1.p₀, (ndims(psde1),1)), p100)[:,1:2:end]
    @test sol2.W.ΔW == ΔW
    @test sol2.W.ΔZ == ΔZ
    @test sol2.ntime == 100
    @test sol2.nsave == 2
    rm(h5file)
end


@testset "$(rpad("SPSDE Solution",80))" begin

    sol = Solution(kubo_oscillator_spsde_1(), Δt, nt)
    @test typeof(sol) <: SolutionPSDE
    @test sol.nd == 1
    @test sol.nm == 1
    @test sol.ns == 1
    @test sol.nt == nt

    sol = Solution(kubo_oscillator_spsde_2(), Δt, nt)
    @test typeof(sol) <: SolutionPSDE
    @test sol.nd == 1
    @test sol.nm == 1
    @test sol.ns == 3
    @test sol.nt == nt

    sol = Solution(kubo_oscillator_spsde_3(), Δt, nt)
    @test typeof(sol) <: SolutionPSDE
    @test sol.nd == 1
    @test sol.nm == 1
    @test sol.ns == 3
    @test sol.nt == nt

end
