using GeometricIntegrators
using GeometricIntegrators.Solutions: save#load, 
using GeometricEquations.Tests.HarmonicOscillator
using Test
import HDF5

ode  = odeproblem()
dae  = daeproblem()
pode = podeproblem()
pdae = pdaeproblem()

nt = 10
Δt = 0.1
ni = 5

t0 = 0.0
x0 = rand(2)
y0 = rand(2)
q0 = rand(1)
p0 = rand(1)
z0 = rand(3)
λ0 = rand(1)
μ0 = rand(2)

t1 = 1.0
x1 = [rand(2) for j = 1:ni]
y1 = [rand(2) for j = 1:ni]
q1 = [rand(1) for j = 1:ni]
p1 = [rand(1) for j = 1:ni]
z1 = [rand(3) for j = 1:ni]
λ1 = [rand(1) for j = 1:ni]
μ1 = [rand(2) for j = 1:ni]

x100 = rand(length(vec(ode.ics.q)), 100)
y100 = rand(length(vec(ode.ics.q)), 100)
z100 = rand(length(vec(dae.ics.q)), 100)
λ100 = rand(length(vec(dae.ics.λ)), 100)
q100 = rand(length(vec(pode.ics.q)), 100)
p100 = rand(length(vec(pode.ics.p)), 100)
μ100 = rand(length(vec(pdae.ics.λ)), 100)

h5file = "test.hdf5"


@testset "$(rpad("HDF5 I/O for ODE Solution",80))" begin
    # test general hdf5 functions
    sol1 = Solution(ode)
    h5io = SolutionHDF5(h5file, sol1; overwrite = true)
    @test typeof(hdf5(h5io)) == HDF5.File
    save(h5io, sol1)
    close(h5io)
    @test isfile(h5file)

    sol1 = Solution(ode)
    h5io = SolutionHDF5(h5file; overwrite = false)
    @test typeof(hdf5(h5io)) == HDF5.File
    save(h5io, sol1)
    close(h5io)
    @test isfile(h5file)

    rm(h5file)

    # test hdf5 in- and output
    sol1 = Solution(harmonic_oscillator_ode(x0))
    h5io = SolutionHDF5(h5file, sol1)
    save(h5io, sol1)
    close(h5io)
    @test isfile(h5file)
    # h5io = SolutionHDF5(h5file) # TODO
    sol2 = SolutionODE(h5file)
    @test sol1 != sol2
    @test sol1.t == sol2.t
    @test sol1.q == sol2.q
    @test sol1.ntime == sol2.ntime
    @test sol1.nsave == sol2.nsave
    rm(h5file)

    # sol1 = Solution(harmonic_oscillator_ode(x1))
    # h5io = SolutionHDF5(h5file, sol1)
    # save(h5io, sol1)
    # close(h5io)
    # @test isfile(h5file)
    # sol2 = SolutionODE(h5file)
    # @test sol1 != sol2
    # @test sol1.t == sol2.t
    # @test sol1.q == sol2.q
    # @test sol1.ntime == sol2.ntime
    # @test sol1.nsave == sol2.nsave
    # rm(h5file)

    sol1 = Solution(similar(ode; tspan = 10 .* tspan(ode)), nwrite = 10)
    h5io = SolutionHDF5(h5file, sol1)
    for i = 1:10
        for j = 1:10
            set_solution!(sol1, x100[:, (i-1)*10+j], j)
        end
        save(h5io, sol1)
        reset!(sol1)
    end
    close(h5io)
    sol2 = SolutionODE(h5file)
    @test sol2.t.t ≈ Δt .* collect(0:100) atol = 1E-14
    @test sol2.q.d == vcat([ode.ics.q], [x100[:, j] for j in axes(x100, 2)])
    @test sol2.ntime == 100
    @test sol2.nsave == 1
    rm(h5file)

    sol1 = Solution(similar(ode; tspan = 10 .* tspan(ode)), nsave = 2, nwrite = 10)
    h5io = SolutionHDF5(h5file, sol1)
    for i = 1:10
        for j = 1:10
            set_solution!(sol1, x100[:, (i-1)*10+j], j)
        end
        save(h5io, sol1)
        reset!(sol1)
    end
    close(h5io)
    sol2 = SolutionODE(h5file)
    @test sol2.t.t ≈ Δt .* collect(0:2:100) atol = 1E-14
    @test sol2.q.d == vcat([ode.ics.q], [x100[:, j] for j in axes(x100, 2)[2:2:end]])
    @test sol2.ntime == 100
    @test sol2.nsave == 2
    rm(h5file)

end


@testset "$(rpad("HDF5 I/O for PODE Solution",80))" begin
    sol1 = Solution(harmonic_oscillator_pode(q0, p0))
    h5io = SolutionHDF5(h5file, sol1)
    save(h5io, sol1)
    close(h5io)
    @test isfile(h5file)
    sol2 = SolutionPODE(h5file)
    @test sol1 != sol2
    @test sol1.t == sol2.t
    @test sol1.q == sol2.q
    @test sol1.p == sol2.p
    @test sol1.ntime == sol2.ntime
    @test sol1.nsave == sol2.nsave
    rm(h5file)

    # sol1 = Solution(harmonic_oscillator_pode(q1, p1))
    # h5io = SolutionHDF5(h5file, sol1)
    # save(h5io, sol1)
    # close(h5io)
    # @test isfile(h5file)
    # sol2 = SolutionPODE(h5file)
    # @test sol1 != sol2
    # @test sol1.t == sol2.t
    # @test sol1.q == sol2.q
    # @test sol1.p == sol2.p
    # @test sol1.ntime == sol2.ntime
    # @test sol1.nsave == sol2.nsave
    # rm(h5file)

    sol1 = Solution(similar(pode; tspan = 10 .* tspan(pode)), nwrite = 10)
    h5io = SolutionHDF5(h5file, sol1)
    for i = 1:10
        for j = 1:10
            oj = (i - 1) * 10 + j
            set_solution!(sol1, q100[:, oj], p100[:, oj], j)
        end
        save(h5io, sol1)
        reset!(sol1)
    end
    close(h5io)
    sol2 = SolutionPODE(h5file)
    @test sol2.t.t ≈ Δt .* collect(0:100) atol = 1E-14
    @test sol2.q.d == vcat([pode.ics.q], [q100[:, j] for j in axes(q100, 2)])
    @test sol2.p.d == vcat([pode.ics.p], [p100[:, j] for j in axes(p100, 2)])
    @test sol2.ntime == 100
    @test sol2.nsave == 1
    rm(h5file)

    sol1 = Solution(similar(pode; tspan = 10 .* tspan(pode)), nsave = 2, nwrite = 10)
    h5io = SolutionHDF5(h5file, sol1)
    for i = 1:10
        for j = 1:10
            oj = (i - 1) * 10 + j
            set_solution!(sol1, q100[:, oj], p100[:, oj], j)
        end
        save(h5io, sol1)
        reset!(sol1)
    end
    close(h5io)
    sol2 = SolutionPODE(h5file)
    @test sol2.t.t ≈ Δt .* collect(0:2:100) atol = 1E-14
    @test sol2.q.d == vcat([pode.ics.q], [q100[:, j] for j in axes(q100, 2)[2:2:end]])
    @test sol2.p.d == vcat([pode.ics.p], [p100[:, j] for j in axes(p100, 2)[2:2:end]])
    @test sol2.ntime == 100
    @test sol2.nsave == 2
    rm(h5file)
end


@testset "$(rpad("HDF5 I/O for DAE Solution",80))" begin
    sol1 = Solution(harmonic_oscillator_dae(z0, λ0))
    h5io = SolutionHDF5(h5file, sol1)
    save(h5io, sol1)
    close(h5io)
    @test isfile(h5file)
    sol2 = SolutionDAE(h5file)
    @test sol1 != sol2
    @test sol1.t == sol2.t
    @test sol1.q == sol2.q
    @test sol1.λ == sol2.λ
    @test sol1.ntime == sol2.ntime
    @test sol1.nsave == sol2.nsave
    rm(h5file)

    # sol1 = Solution(harmonic_oscillator_dae(z1, λ1))
    # h5io = SolutionHDF5(h5file, sol1)
    # save(h5io, sol1)
    # close(h5io)
    # @test isfile(h5file)
    # sol2 = SolutionDAE(h5file)
    # @test sol1 != sol2
    # @test sol1.t == sol2.t
    # @test sol1.q == sol2.q
    # @test sol1.λ == sol2.λ
    # @test sol1.ntime == sol2.ntime
    # @test sol1.nsave == sol2.nsave
    # rm(h5file)

    sol1 = Solution(similar(dae; tspan = 10 .* tspan(dae)), nwrite = 10)
    h5io = SolutionHDF5(h5file, sol1)
    for i = 1:10
        for j = 1:10
            oj = (i - 1) * 10 + j
            set_solution!(sol1, z100[:, oj], λ100[:, oj], j)
        end
        save(h5io, sol1)
        reset!(sol1)
    end
    close(h5io)
    sol2 = SolutionDAE(h5file)
    @test sol2.t.t ≈ Δt .* collect(0:100) atol = 1E-14
    @test sol2.q.d == vcat([dae.ics.q], [z100[:, j] for j in axes(z100, 2)])
    @test sol2.λ.d == vcat([dae.ics.λ], [λ100[:, j] for j in axes(λ100, 2)])
    @test sol2.ntime == 100
    @test sol2.nsave == 1
    rm(h5file)

    sol1 = Solution(similar(dae; tspan = 10 .* tspan(dae)), nsave = 2, nwrite = 10)
    h5io = SolutionHDF5(h5file, sol1)
    for i = 1:10
        for j = 1:10
            oj = (i - 1) * 10 + j
            set_solution!(sol1, z100[:, oj], λ100[:, oj], j)
        end
        save(h5io, sol1)
        reset!(sol1)
    end
    close(h5io)
    sol2 = SolutionDAE(h5file)
    @test sol2.t.t ≈ Δt .* collect(0:2:100) atol = 1E-14
    @test sol2.q.d == vcat([dae.ics.q], [z100[:, j] for j in axes(z100, 2)[2:2:end]])
    @test sol2.λ.d == vcat([dae.ics.λ], [λ100[:, j] for j in axes(λ100, 2)[2:2:end]])
    @test sol2.ntime == 100
    @test sol2.nsave == 2
    rm(h5file)
end


@testset "$(rpad("HDF5 I/O for PDAE Solution",80))" begin
    sol1 = Solution(harmonic_oscillator_pdae(x0, y0, μ0))
    h5io = SolutionHDF5(h5file, sol1)
    save(h5io, sol1)
    close(h5io)
    @test isfile(h5file)
    sol2 = SolutionPDAE(h5file)
    @test sol1 != sol2
    @test sol1.t == sol2.t
    @test sol1.q == sol2.q
    @test sol1.p == sol2.p
    @test sol1.λ == sol2.λ
    @test sol1.ntime == sol2.ntime
    @test sol1.nsave == sol2.nsave
    rm(h5file)

    # sol1 = Solution(harmonic_oscillator_pdae(x1, y1, μ1))
    # h5io = SolutionHDF5(h5file, sol1)
    # save(h5io, sol1)
    # close(h5io)
    # @test isfile(h5file)
    # sol2 = SolutionPDAE(h5file)
    # @test sol1 != sol2
    # @test sol1.t == sol2.t
    # @test sol1.q == sol2.q
    # @test sol1.p == sol2.p
    # @test sol1.λ == sol2.λ
    # @test sol1.ntime == sol2.ntime
    # @test sol1.nsave == sol2.nsave
    # rm(h5file)

    sol1 = Solution(similar(pdae; tspan = 10 .* tspan(pdae)), nwrite = 10)
    h5io = SolutionHDF5(h5file, sol1)
    for i = 1:10
        for j = 1:10
            oj = (i - 1) * 10 + j
            set_solution!(sol1, x100[:, oj], y100[:, oj], μ100[:, oj], j)
        end
        save(h5io, sol1)
        reset!(sol1)
    end
    close(h5io)
    sol2 = SolutionPDAE(h5file)
    @test sol2.t.t ≈ Δt .* collect(0:100) atol = 1E-14
    @test sol2.q.d == vcat([pdae.ics.q], [x100[:, j] for j in axes(x100, 2)])
    @test sol2.p.d == vcat([pdae.ics.p], [y100[:, j] for j in axes(y100, 2)])
    @test sol2.λ.d == vcat([pdae.ics.λ], [μ100[:, j] for j in axes(μ100, 2)])
    # @test sol2.q.d == hcat(reshape(pdae.q₀, (ndims(pdae),1)), x100)
    # @test sol2.p.d == hcat(reshape(pdae.p₀, (ndims(pdae),1)), y100)
    # @test sol2.λ.d == hcat(reshape(pdae.λ₀, (pdae.m,1)), μ100)
    @test sol2.ntime == 100
    @test sol2.nsave == 1
    rm(h5file)

    sol1 = Solution(similar(pdae; tspan = 10 .* tspan(pdae)), nsave = 2, nwrite = 10)
    h5io = SolutionHDF5(h5file, sol1)
    for i = 1:10
        for j = 1:10
            oj = (i - 1) * 10 + j
            set_solution!(sol1, x100[:, oj], y100[:, oj], μ100[:, oj], j)
        end
        save(h5io, sol1)
        reset!(sol1)
    end
    close(h5io)
    sol2 = SolutionPDAE(h5file)
    @test sol2.t.t ≈ Δt .* collect(0:2:100) atol = 1E-14
    @test sol2.q.d == vcat([pdae.ics.q], [x100[:, j] for j in axes(x100, 2)[2:2:end]])
    @test sol2.p.d == vcat([pdae.ics.p], [y100[:, j] for j in axes(y100, 2)[2:2:end]])
    @test sol2.λ.d == vcat([pdae.ics.λ], [μ100[:, j] for j in axes(μ100, 2)[2:2:end]])
    # @test sol2.q.d == hcat(reshape(pdae.q₀, (ndims(pdae),1)), x100)[:,1:2:end]
    # @test sol2.p.d == hcat(reshape(pdae.p₀, (ndims(pdae),1)), y100)[:,1:2:end]
    # @test sol2.λ.d == hcat(reshape(pdae.λ₀, (pdae.m,1)), μ100)[:,1:2:end]
    @test sol2.ntime == 100
    @test sol2.nsave == 2
    rm(h5file)
end
