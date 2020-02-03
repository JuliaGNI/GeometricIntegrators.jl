
using GeometricIntegrators.CommonFunctions
using GeometricIntegrators.Solutions
using GeometricIntegrators.Solutions: createHDF5
using GeometricIntegrators.TestProblems.KuboOscillatorProblem
using HDF5: HDF5File
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

sde  = kubo_oscillator_sde_1()
psde = kubo_oscillator_psde_1()

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
    sol = Solution(sde, Δt, nt)
    @test typeof(sol) <: SolutionSDE

    sol0 = Solution(similar(sde, x0), Δt, nt)
    @test typeof(sol0) <: SolutionSDE

    sol1 = Solution(similar(sde, x1), Δt, nt)
    @test typeof(sol1) <: SolutionSDE

    @test sol != sol0
    @test sol != sol1

    @test sol == SolutionSDE(sde, Δt, sol.W.ΔW, sol.W.ΔZ, sol.ntime)

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
    sol1 = Solution(similar(kubo_oscillator_sde_2(), x0), Δt, nt)
    sol2 = Solution(similar(kubo_oscillator_sde_2(), x0), Δt, nt)
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

    # test reset
    reset!(sol)
    @test sol.t[0]   == t1
    @test sol.t[end] == t2
    @test offset(sol) == nt

    # test hdf5 in- and output
    sol1 = Solution(kubo_oscillator_sde_2(x0), Δt, nt)
    create_hdf5(sol1, h5file)
    write_to_hdf5(sol1)
    close(sol1)
    @test isfile(h5file)
    sol2 = SolutionSDE(h5file)
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

    sol1 = Solution(kubo_oscillator_sde_2(x1), Δt, nt)
    create_hdf5(sol1, h5file)
    write_to_hdf5(sol1)
    close(sol1)
    @test isfile(h5file)
    sol2 = SolutionSDE(h5file)
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

    # test nsave and nwrite parameters
    sol = Solution(sde, Δt, 20, 2)
    @test sol.nt == 10

    sol = Solution(sde, Δt, 20, 2, 10)
    @test sol.nt == 5
end


@testset "$(rpad("PSDE Solution",80))" begin
    psde  = kubo_oscillator_psde_1()
    sol = Solution(psde, Δt, nt)
    @test typeof(sol) <: SolutionPSDE

    # test hdf5 in- and output
    h5 = createHDF5(sol, h5file)
    @test typeof(h5) == HDF5File
    close(h5)
    @test isfile(h5file)

    h5 = createHDF5(sol, h5file, overwrite=false)
    @test typeof(h5) == HDF5File
    close(h5)
    @test isfile(h5file)
    rm(h5file)
end


@testset "$(rpad("SPSDE Solution",80))" begin
    spsde  = kubo_oscillator_spsde_1()
    sol = Solution(spsde, Δt, nt)
    @test typeof(sol) <: SolutionPSDE
end
