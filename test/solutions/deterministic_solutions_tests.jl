
using GeometricIntegrators.CommonFunctions
using GeometricIntegrators.Solutions
using GeometricIntegrators.Solutions: createHDF5
using GeometricIntegrators.TestProblems.HarmonicOscillator
using HDF5: HDF5File
using Test


ntime = 10
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
x2    = rand(2, ni)

ode   = oscillator_ode()
dae   = oscillator_dae()
pode  = oscillator_pode()
pdae  = oscillator_pdae()

h5file = "test.hdf5"


@testset "$(rpad("ODE Solution",80))" begin
    sol = Solution(ode, Δt, ntime)
    @test typeof(sol) <: SolutionODE

    sol0 = Solution(similar(ode, x0), Δt, ntime)
    @test typeof(sol0) <: SolutionODE

    sol1 = Solution(similar(ode, x1), Δt, ntime)
    @test typeof(sol1) <: SolutionODE

    @test sol != sol0
    @test sol != sol1

    set_initial_conditions!(sol, t0, x0)
    get_initial_conditions!(sol, tx, 1)

    # @test sol != sol0
    @test tx == x0

    set_initial_conditions!(sol1, similar(ode, t1, x2))
    get_initial_conditions!(sol1, tx, 1)
    @test tx == x2[:,1]

    # test hdf5 in- and output
    h5 = createHDF5(sol, h5file)
    @test typeof(h5) == HDF5File
    close(h5)
    @test isfile(h5file)

    h5 = createHDF5(sol, h5file; overwrite=false)
    @test typeof(h5) == HDF5File
    close(h5)
    @test isfile(h5file)
    rm(h5file)

    write_to_hdf5(sol, h5file)
    @test isfile(h5file)
    rm(h5file)

    create_hdf5(sol, h5file)
    write_to_hdf5(sol)
    close(sol)

    sol2 = SolutionODE(h5file)
    @test sol != sol2
    rm(h5file)

    # test nsave and nwrite parameters
    sol = Solution(ode, Δt, 20, 2)
    @test sol.nt == 10

    sol = Solution(ode, Δt, 20, 2, 10)
    @test sol.nt == 5
end


@testset "$(rpad("PODE Solution",80))" begin
    psol = Solution(pode, Δt, ntime)
    @test typeof(psol) <: SolutionPODE

    psol0 = Solution(similar(pode, q0, p0), Δt, ntime)
    @test typeof(psol0) <: SolutionPODE

    psol1 = Solution(similar(pode, q1, p1), Δt, ntime)
    @test typeof(psol1) <: SolutionPODE

    @test psol != psol0
    @test psol != psol1

    # test hdf5 in- and output
    create_hdf5(psol, "test.hdf5")
    write_to_hdf5(psol)
    close(psol)

    psol2 = SolutionPODE("test.hdf5")
    @test psol != psol2
    rm("test.hdf5")
end


@testset "$(rpad("DAE Solution",80))" begin
    sol = Solution(dae, Δt, ntime)
    @test typeof(sol) <: SolutionDAE

    sol0 = Solution(similar(dae, z0, λ0), Δt, ntime)
    @test typeof(sol0) <: SolutionDAE

    sol1 = Solution(similar(dae, z1, λ1), Δt, ntime)
    @test typeof(sol1) <: SolutionDAE

    @test sol != sol0
    @test sol != sol1

    # test hdf5 in- and output
    create_hdf5(sol, "test.hdf5")
    write_to_hdf5(sol)
    close(sol)

    sol2 = SolutionDAE("test.hdf5")
    @test sol != sol2
    rm("test.hdf5")
end


@testset "$(rpad("PDAE Solution",80))" begin
    psol = Solution(pdae, Δt, ntime)
    @test typeof(psol) <: SolutionPDAE

    psol0 = Solution(similar(pdae, x0, x0.^2, zero(x0)), Δt, ntime)
    @test typeof(psol0) <: SolutionPDAE

    psol1 = Solution(similar(pdae, x1, x1.^2, zero(x1)), Δt, ntime)
    @test typeof(psol1) <: SolutionPDAE

    @test psol != psol0
    @test psol != psol1

    # test hdf5 in- and output
    create_hdf5(psol, "test.hdf5")
    write_to_hdf5(psol)
    close(psol)

    psol2 = SSolutionPDAE("test.hdf5")
    @test psol != psol2
    rm("test.hdf5")
end
