
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
y0    = rand(2)
q0    = rand(1)
p0    = rand(1)
z0    = rand(3)
λ0    = rand(1)
μ0    = rand(2)

t1    = 1.
x1    = rand(2, ni)
y1    = rand(2, ni)
q1    = rand(1, ni)
p1    = rand(1, ni)
z1    = rand(3, ni)
λ1    = rand(1, ni)
μ1    = rand(2, ni)

t2    = t1 + (t1-t0)
x2    = rand(2, ni)
y2    = rand(2, ni)
q2    = rand(1, ni)
p2    = rand(1, ni)
z2    = rand(3, ni)
λ2    = rand(1, ni)
μ2    = rand(2, ni)

tx    = zero(x0)
ty    = zero(y0)
tq    = zero(q0)
tp    = zero(p0)
tz    = zero(z0)
tλ    = zero(λ0)
tμ    = zero(μ0)


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
    for i in 1:ni
        get_initial_conditions!(sol1, tx, i)
        @test tx == x2[:,i]
    end

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
    @test isfile(h5file)
    sol2 = SolutionODE(h5file)
    @test sol1   != sol2
    @test sol1.t == sol2.t
    @test sol1.q == sol2.q
    @test sol1.ntime == sol2.ntime
    @test sol1.nsave == sol2.nsave
    rm(h5file)

    sol1 = Solution(harmonic_oscillator_ode(x1), Δt, nt)
    create_hdf5(sol1, h5file)
    write_to_hdf5(sol1)
    close(sol1)
    @test isfile(h5file)
    sol2 = SolutionODE(h5file)
    @test sol1   != sol2
    @test sol1.t == sol2.t
    @test sol1.q == sol2.q
    @test sol1.ntime == sol2.ntime
    @test sol1.nsave == sol2.nsave
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
    for i in 1:ni
        get_initial_conditions!(sol1, tq, tp, i)
        @test tq == q2[:,i]
        @test tp == p2[:,i]
    end

    reset!(sol)
    @test sol.t[0]   == t1
    @test sol.t[end] == t2
    @test offset(sol) == nt

    # test hdf5 in- and output
    sol1 = Solution(harmonic_oscillator_pode(q0, p0), Δt, nt)
    create_hdf5(sol1, h5file)
    write_to_hdf5(sol1)
    close(sol1)
    @test isfile(h5file)
    sol2 = SolutionPODE(h5file)
    @test sol1   != sol2
    @test sol1.t == sol2.t
    @test sol1.q == sol2.q
    @test sol1.p == sol2.p
    @test sol1.ntime == sol2.ntime
    @test sol1.nsave == sol2.nsave
    rm(h5file)

    sol1 = Solution(harmonic_oscillator_pode(q1, p1), Δt, nt)
    create_hdf5(sol1, h5file)
    write_to_hdf5(sol1)
    close(sol1)
    @test isfile(h5file)
    sol2 = SolutionPODE(h5file)
    @test sol1   != sol2
    @test sol1.t == sol2.t
    @test sol1.q == sol2.q
    @test sol1.p == sol2.p
    @test sol1.ntime == sol2.ntime
    @test sol1.nsave == sol2.nsave
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

    set_initial_conditions!(sol, t0, z0, λ0)
    get_initial_conditions!(sol, tz, tλ, 1)

    @test tz == z0
    @test tλ == λ0

    δt, δz, δλ = get_initial_conditions(sol, 1)
    @test δt == t0
    @test δz == z0
    @test δλ == λ0

    set_initial_conditions!(sol1, similar(dae, t1, z2, λ2))
    for i in 1:ni
        get_initial_conditions!(sol1, tz, tλ, i)
        @test tz == z2[:,i]
        @test tλ == λ2[:,i]
    end

    reset!(sol)
    @test sol.t[0]   == t1
    @test sol.t[end] == t2
    @test offset(sol) == nt

    # test hdf5 in- and output
    sol1 = Solution(harmonic_oscillator_dae(z0, λ0), Δt, nt)
    create_hdf5(sol1, h5file)
    write_to_hdf5(sol1)
    close(sol1)
    @test isfile(h5file)
    sol2 = SolutionDAE(h5file)
    @test sol1   != sol2
    @test sol1.t == sol2.t
    @test sol1.q == sol2.q
    @test sol1.λ == sol2.λ
    @test sol1.ntime == sol2.ntime
    @test sol1.nsave == sol2.nsave
    rm(h5file)

    sol1 = Solution(harmonic_oscillator_dae(z1, λ1), Δt, nt)
    create_hdf5(sol1, h5file)
    write_to_hdf5(sol1)
    close(sol1)
    @test isfile(h5file)
    sol2 = SolutionDAE(h5file)
    @test sol1   != sol2
    @test sol1.t == sol2.t
    @test sol1.q == sol2.q
    @test sol1.λ == sol2.λ
    @test sol1.ntime == sol2.ntime
    @test sol1.nsave == sol2.nsave
    rm(h5file)
end


@testset "$(rpad("PDAE Solution",80))" begin
    sol = Solution(pdae, Δt, nt)
    @test typeof(sol) <: SolutionPDAE

    sol0 = Solution(similar(pdae, x0, y0, μ0), Δt, nt)
    @test typeof(sol0) <: SolutionPDAE

    sol1 = Solution(similar(pdae, x1, y1, μ1), Δt, nt)
    @test typeof(sol1) <: SolutionPDAE

    @test sol != sol0
    @test sol != sol1

    @test nsave(sol) == 1
    @test ntime(sol) == nt
    @test timesteps(sol) == Δt .* collect(0:nt)

    set_initial_conditions!(sol, t0, x0, y0, μ0)
    get_initial_conditions!(sol, tx, ty, tμ, 1)

    @test tx == x0
    @test ty == y0
    @test tμ == μ0

    δt, δx, δy, δμ = get_initial_conditions(sol, 1)
    @test δt == t0
    @test δx == x0
    @test δy == y0
    @test δμ == μ0

    set_initial_conditions!(sol1, similar(pdae, t1, x2, y2, μ2))
    for i in 1:ni
        get_initial_conditions!(sol1, tx, ty, tμ, i)
        @test tx == x2[:,i]
        @test ty == y2[:,i]
        @test tμ == μ2[:,i]#
    end

    reset!(sol)
    @test sol.t[0]   == t1
    @test sol.t[end] == t2
    @test offset(sol) == nt

    # test hdf5 in- and output
    sol1 = Solution(harmonic_oscillator_pdae(x0, y0, μ0), Δt, nt)
    create_hdf5(sol1, h5file)
    write_to_hdf5(sol1)
    close(sol1)
    @test isfile(h5file)
    sol2 = SSolutionPDAE(h5file)
    @test sol1   != sol2
    @test sol1.t == sol2.t
    @test sol1.q == sol2.q
    @test sol1.p == sol2.p
    @test sol1.λ == sol2.λ
    @test sol1.ntime == sol2.ntime
    @test sol1.nsave == sol2.nsave
    rm(h5file)

    sol1 = Solution(harmonic_oscillator_pdae(x1, y1, μ1), Δt, nt)
    create_hdf5(sol1, h5file)
    write_to_hdf5(sol1)
    close(sol1)
    @test isfile(h5file)
    sol2 = SSolutionPDAE(h5file)
    @test sol1   != sol2
    @test sol1.t == sol2.t
    @test sol1.q == sol2.q
    @test sol1.p == sol2.p
    @test sol1.λ == sol2.λ
    @test sol1.ntime == sol2.ntime
    @test sol1.nsave == sol2.nsave
    rm(h5file)
end
