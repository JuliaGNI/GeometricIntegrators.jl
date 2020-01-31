
using GeometricIntegrators.CommonFunctions
using GeometricIntegrators.Solutions
using GeometricIntegrators.Solutions: createHDF5
using GeometricIntegrators.TestProblems.HarmonicOscillatorProblem
using HDF5: HDF5File
using Test


nt    = 10
Δt    = .1
ni    = 5

t0    = 0.
x0    = rand(2)
q0    = rand(1)
z0    = [x0[1], x0[2], x0[1] + x0[2]]
p0    = q0.^2
λ0    = [0.]

t1    = 1.
x1    = rand(2, ni)
q1    = rand(1, ni)
p1    = q1.^2
λ1    = zeros(1, ni)

z1    = rand(3, ni)
z1[3,:] .= z0[1,:] .+ z0[2,:]

tx    = zero(x0)
tq    = zero(q0)
tp    = zero(p0)

t2    = t1 + (t1-t0)
x2    = rand(2, ni)
q2    = rand(1, ni)
p2    = rand(1, ni)


ode   = harmonic_oscillator_ode()
dae   = harmonic_oscillator_dae()
pode  = harmonic_oscillator_pode()
pdae  = harmonic_oscillator_pdae()

h5file = "test.hdf5"


@testset "$(rpad("ODE Solution",80))" begin
    sol = Solution(ode, Δt, nt)
    @test typeof(sol) <: SolutionODE

    sol0 = Solution(similar(ode, x0), Δt, nt)
    @test typeof(sol0) <: SolutionODE

    sol1 = Solution(similar(ode, x1), Δt, nt)
    @test typeof(sol1) <: SolutionODE

    @test sol != sol0
    @test sol != sol1

    @test nsave(sol) == 1
    @test ntime(sol) == nt
    @test timesteps(sol) == Δt .* collect(0:nt)

    set_initial_conditions!(sol, t0, x0)
    get_initial_conditions!(sol, tx, 1)

    @test tx == x0

    δt, δx = get_initial_conditions(sol, 1)
    @test δt == t0
    @test δx == x0

    set_initial_conditions!(sol1, similar(ode, t1, x2))
    get_initial_conditions!(sol1, tx, 1)
    @test tx == x2[:,1]

    reset!(sol)
    @test sol.t[0]   == t1
    @test sol.t[end] == t2
    @test offset(sol) == nt

    # test hdf5 in- and output
    sol1 = Solution(ode, Δt, nt)
    h5 = createHDF5(sol1, h5file)
    @test typeof(h5) == HDF5File
    close(h5)
    @test isfile(h5file)

    sol1 = Solution(ode, Δt, nt)
    h5 = createHDF5(sol1, h5file; overwrite=false)
    @test typeof(h5) == HDF5File
    close(h5)
    @test isfile(h5file)
    rm(h5file)

    sol1 = Solution(harmonic_oscillator_ode(x0), Δt, nt)
    create_hdf5(sol1, h5file)
    write_to_hdf5(sol1)
    close(sol1)
    sol2 = SolutionODE(h5file)
    @test sol1   != sol2
    @test sol1.t == sol2.t
    @test sol1.q == sol2.q
    rm(h5file)

    sol1 = Solution(harmonic_oscillator_ode(x1), Δt, nt)
    create_hdf5(sol1, h5file)
    write_to_hdf5(sol1)
    close(sol1)
    sol2 = SolutionODE(h5file)
    @test sol1   != sol2
    @test sol1.t == sol2.t
    @test sol1.q == sol2.q
    rm(h5file)

    # test nsave and nwrite parameters
    sol = Solution(ode, Δt, 20, 2)
    @test sol.nt == 10

    sol = Solution(ode, Δt, 20, 2, 10)
    @test sol.nt == 5
end


@testset "$(rpad("PODE Solution",80))" begin
    sol = Solution(pode, Δt, nt)
    @test typeof(sol) <: SolutionPODE

    sol0 = Solution(similar(pode, q0, p0), Δt, nt)
    @test typeof(sol0) <: SolutionPODE

    sol1 = Solution(similar(pode, q1, p1), Δt, nt)
    @test typeof(sol1) <: SolutionPODE

    @test sol != sol0
    @test sol != sol1

    @test nsave(sol) == 1
    @test ntime(sol) == nt
    @test timesteps(sol) == Δt .* collect(0:nt)

    set_initial_conditions!(sol, t0, q0, p0)
    get_initial_conditions!(sol, tq, tp, 1)

    @test tq == q0
    @test tp == p0

    δt, δq, δp = get_initial_conditions(sol, 1)
    @test δt == t0
    @test δq == q0
    @test δp == p0

    set_initial_conditions!(sol1, similar(pode, t1, q2, p2))
    get_initial_conditions!(sol1, tq, tp, 1)
    @test tq == q2[:,1]
    @test tp == p2[:,1]

    reset!(sol)
    @test sol.t[0]   == t1
    @test sol.t[end] == t2
    @test offset(sol) == nt

    # test hdf5 in- and output
    sol1 = Solution(pode, Δt, nt)
    h5 = createHDF5(sol1, h5file)
    @test typeof(h5) == HDF5File
    close(h5)
    @test isfile(h5file)

    sol1 = Solution(pode, Δt, nt)
    h5 = createHDF5(sol1, h5file; overwrite=false)
    @test typeof(h5) == HDF5File
    close(h5)
    @test isfile(h5file)
    rm(h5file)

    sol1 = Solution(harmonic_oscillator_pode(q0, p0), Δt, nt)
    create_hdf5(sol1, h5file)
    write_to_hdf5(sol1)
    close(sol1)
    sol2 = SolutionPODE(h5file)
    @test sol1   != sol2
    @test sol1.t == sol2.t
    @test sol1.q == sol2.q
    @test sol1.p == sol2.p
    rm(h5file)

    sol1 = Solution(harmonic_oscillator_pode(q1, p1), Δt, nt)
    create_hdf5(sol1, h5file)
    write_to_hdf5(sol1)
    close(sol1)
    sol2 = SolutionPODE(h5file)
    @test sol1   != sol2
    @test sol1.t == sol2.t
    @test sol1.q == sol2.q
    @test sol1.p == sol2.p
    rm(h5file)
end


@testset "$(rpad("DAE Solution",80))" begin
    sol = Solution(dae, Δt, nt)
    @test typeof(sol) <: SolutionDAE

    sol0 = Solution(similar(dae, z0, λ0), Δt, nt)
    @test typeof(sol0) <: SolutionDAE

    sol1 = Solution(similar(dae, z1, λ1), Δt, nt)
    @test typeof(sol1) <: SolutionDAE

    @test sol != sol0
    @test sol != sol1

    @test nsave(sol) == 1
    @test ntime(sol) == nt
    @test timesteps(sol) == Δt .* collect(0:nt)

    # test hdf5 in- and output
    create_hdf5(sol, "test.hdf5")
    write_to_hdf5(sol)
    close(sol)

    sol2 = SolutionDAE("test.hdf5")
    @test sol != sol2
    rm("test.hdf5")
end


@testset "$(rpad("PDAE Solution",80))" begin
    sol = Solution(pdae, Δt, nt)
    @test typeof(sol) <: SolutionPDAE

    sol0 = Solution(similar(pdae, x0, x0.^2, zero(x0)), Δt, nt)
    @test typeof(sol0) <: SolutionPDAE

    sol1 = Solution(similar(pdae, x1, x1.^2, zero(x1)), Δt, nt)
    @test typeof(sol1) <: SolutionPDAE

    @test sol != sol0
    @test sol != sol1

    @test nsave(sol) == 1
    @test ntime(sol) == nt
    @test timesteps(sol) == Δt .* collect(0:nt)

    # test hdf5 in- and output
    create_hdf5(sol, "test.hdf5")
    write_to_hdf5(sol)
    close(sol)

    sol2 = SSolutionPDAE("test.hdf5")
    @test sol != sol2
    rm("test.hdf5")
end
