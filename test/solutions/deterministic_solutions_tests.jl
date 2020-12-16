
using GeometricIntegrators.CommonFunctions
using GeometricIntegrators.Solutions
using GeometricIntegrators.Solutions: createHDF5
using GeometricProblems.HarmonicOscillator
using Test
import HDF5


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

xs    = rand(2,nt)
ys    = rand(2,nt)
qs    = rand(1,nt)
ps    = rand(1,nt)
zs    = rand(3,nt)
λs    = rand(1,nt)
μs    = rand(2,nt)

Xs    = rand(2,nt,ni)
Ys    = rand(2,nt,ni)
Qs    = rand(1,nt,ni)
Ps    = rand(1,nt,ni)
Zs    = rand(3,nt,ni)
Λs    = rand(1,nt,ni)
Ms    = rand(2,nt,ni)

ode   = harmonic_oscillator_ode()
dae   = harmonic_oscillator_dae()
pode  = harmonic_oscillator_pode()
pdae  = harmonic_oscillator_pdae()

x100  = rand(ndims(ode),  100)
y100  = rand(ndims(ode),  100)
z100  = rand(ndims(dae),  100)
q100  = rand(ndims(pode), 100)
p100  = rand(ndims(pode), 100)
λ100  = rand(dae.m,  100)
μ100  = rand(pdae.m, 100)

h5file = "test.hdf5"


@testset "$(rpad("ODE Solution",80))" begin
    asol = AtomicSolution(ode)

    # test constructors and general functionality
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

    set_initial_conditions!(sol1, similar(ode, t1, x2))
    for i in 1:ni
        get_initial_conditions!(sol1, tx, i)
        @test tx == x2[:,i]
    end

    # test set/get solution
    sol1 = Solution(similar(ode, x0), Δt, nt)
    sol2 = Solution(similar(ode, x0), Δt, nt)
    for i in 1:nt
        tx .= xs[:,i]
        asol.q .= xs[:,i]
        set_solution!(sol1, tx, i)
        set_solution!(sol2, asol, i)
    end
    @test sol1.q[:,1:nt] == xs
    @test sol2.q[:,1:nt] == xs

    sol1 = Solution(similar(ode, x1), Δt, nt)
    sol2 = Solution(similar(ode, x1), Δt, nt)
    for i in 1:nt
        for k in 1:ni
            tx .= Xs[:,i,k]
            asol.q .= Xs[:,i,k]
            set_solution!(sol1, tx, i, k)
            set_solution!(sol2, asol, i, k)
        end
    end
    @test sol1.q[:,1:nt,:] == Xs
    @test sol2.q[:,1:nt,:] == Xs

    # test general hdf5 functions
    sol1 = Solution(ode, Δt, nt)
    h5 = createHDF5(sol1, h5file)
    @test typeof(h5) == HDF5.File
    close(h5)
    @test isfile(h5file)

    sol1 = Solution(ode, Δt, nt)
    h5 = createHDF5(sol1, h5file; overwrite=false)
    @test typeof(h5) == HDF5.File
    close(h5)
    @test isfile(h5file)
    rm(h5file)

    sol1 = Solution(ode, Δt, nt)
    write_to_hdf5(sol1, h5file)
    @test isfile(h5file)
    rm(h5file)

    # test hdf5 in- and output
    sol1 = Solution(harmonic_oscillator_ode(x0), Δt, nt)
    create_hdf5!(sol1, h5file)
    write_to_hdf5(sol1)
    close(sol1)
    @test isfile(h5file)
    sol2 = SSolutionODE(h5file)
    @test sol1   != sol2
    @test sol1.t == sol2.t
    @test sol1.q == sol2.q
    @test sol1.ntime == sol2.ntime
    @test sol1.nsave == sol2.nsave
    rm(h5file)

    sol1 = Solution(harmonic_oscillator_ode(x1), Δt, nt)
    create_hdf5!(sol1, h5file)
    write_to_hdf5(sol1)
    close(sol1)
    @test isfile(h5file)
    sol2 = SSolutionODE(h5file)
    @test sol1   != sol2
    @test sol1.t == sol2.t
    @test sol1.q == sol2.q
    @test sol1.ntime == sol2.ntime
    @test sol1.nsave == sol2.nsave
    rm(h5file)

    # test nsave and nwrite parameters
    sol = Solution(ode, Δt, 20, nsave=2)
    @test sol.nt == 10

    sol = Solution(ode, Δt, 20, nsave=2, nwrite=10)
    @test sol.nt == 5

    # test reset
    sol = Solution(ode, Δt, nt)
    reset!(sol)
    @test sol.t[0]   == t1
    @test sol.t[end] == t2
    @test offset(sol) == nt

    sol1 = Solution(ode, Δt, 100, nwrite=10)
    create_hdf5!(sol1, h5file)
    for i in 1:10
        for j in 1:10
            set_solution!(sol1, x100[:,(i-1)*10+j], j)
        end
        write_to_hdf5(sol1)
        reset!(sol1)
    end
    close(sol1)
    sol2 = SSolutionODE(h5file)
    @test sol2.t.t ≈ Δt .* collect(0:100) atol=1E-14
    @test sol2.q.d == hcat(reshape(ode.q₀, (ndims(ode),1)), x100)
    @test sol2.ntime == 100
    @test sol2.nsave == 1
    rm(h5file)

    sol1 = Solution(ode, Δt, 100, nsave=2, nwrite=10)
    create_hdf5!(sol1, h5file)
    for i in 1:10
        for j in 1:10
            set_solution!(sol1, x100[:,(i-1)*10+j], j)
        end
        write_to_hdf5(sol1)
        reset!(sol1)
    end
    close(sol1)
    sol2 = SSolutionODE(h5file)
    @test sol2.t.t ≈ Δt .* collect(0:2:100) atol=1E-14
    @test sol2.q.d == hcat(reshape(ode.q₀, (ndims(ode),1)), x100)[:,1:2:end]
    @test sol2.ntime == 100
    @test sol2.nsave == 2
    rm(h5file)
end


@testset "$(rpad("PODE Solution",80))" begin
    asol = AtomicSolution(pode)

    # test constructors and general functionality
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

    set_initial_conditions!(sol1, similar(pode, t1, q2, p2))
    for i in 1:ni
        get_initial_conditions!(sol1, tq, tp, i)
        @test tq == q2[:,i]
        @test tp == p2[:,i]
    end

    # test set/get solution
    sol1 = Solution(similar(pode, q0, p0), Δt, nt)
    sol2 = Solution(similar(pode, q0, p0), Δt, nt)
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

    sol1 = Solution(similar(pode, q1, p1), Δt, nt)
    sol2 = Solution(similar(pode, q1, p1), Δt, nt)
    for i in 1:nt
        for k in 1:ni
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
    sol1 = Solution(harmonic_oscillator_pode(q0, p0), Δt, nt)
    create_hdf5!(sol1, h5file)
    write_to_hdf5(sol1)
    close(sol1)
    @test isfile(h5file)
    sol2 = SSolutionPODE(h5file)
    @test sol1   != sol2
    @test sol1.t == sol2.t
    @test sol1.q == sol2.q
    @test sol1.p == sol2.p
    @test sol1.ntime == sol2.ntime
    @test sol1.nsave == sol2.nsave
    rm(h5file)

    sol1 = Solution(harmonic_oscillator_pode(q1, p1), Δt, nt)
    create_hdf5!(sol1, h5file)
    write_to_hdf5(sol1)
    close(sol1)
    @test isfile(h5file)
    sol2 = SSolutionPODE(h5file)
    @test sol1   != sol2
    @test sol1.t == sol2.t
    @test sol1.q == sol2.q
    @test sol1.p == sol2.p
    @test sol1.ntime == sol2.ntime
    @test sol1.nsave == sol2.nsave
    rm(h5file)

    # test nsave and nwrite parameters
    sol = Solution(pode, Δt, 20, nsave=2)
    @test sol.nt == 10

    sol = Solution(pode, Δt, 20, nsave=2, nwrite=10)
    @test sol.nt == 5

    # test reset
    sol = Solution(pode, Δt, nt)
    reset!(sol)
    @test sol.t[0]   == t1
    @test sol.t[end] == t2
    @test offset(sol) == nt

    sol1 = Solution(pode, Δt, 100, nwrite=10)
    create_hdf5!(sol1, h5file)
    for i in 1:10
        for j in 1:10
            oj = (i-1)*10+j
            set_solution!(sol1, q100[:,oj], p100[:,oj], j)
        end
        write_to_hdf5(sol1)
        reset!(sol1)
    end
    close(sol1)
    sol2 = SSolutionPODE(h5file)
    @test sol2.t.t ≈ Δt .* collect(0:100) atol=1E-14
    @test sol2.q.d == hcat(reshape(pode.q₀, (ndims(pode),1)), q100)
    @test sol2.p.d == hcat(reshape(pode.p₀, (ndims(pode),1)), p100)
    @test sol2.ntime == 100
    @test sol2.nsave == 1
    rm(h5file)

    sol1 = Solution(pode, Δt, 100, nsave=2, nwrite=10)
    create_hdf5!(sol1, h5file)
    for i in 1:10
        for j in 1:10
            oj = (i-1)*10+j
            set_solution!(sol1, q100[:,oj], p100[:,oj], j)
        end
        write_to_hdf5(sol1)
        reset!(sol1)
    end
    close(sol1)
    sol2 = SSolutionPODE(h5file)
    @test sol2.t.t ≈ Δt .* collect(0:2:100) atol=1E-14
    @test sol2.q.d == hcat(reshape(pode.q₀, (ndims(pode),1)), q100)[:,1:2:end]
    @test sol2.p.d == hcat(reshape(pode.p₀, (ndims(pode),1)), p100)[:,1:2:end]
    @test sol2.ntime == 100
    @test sol2.nsave == 2
    rm(h5file)
end


@testset "$(rpad("DAE Solution",80))" begin
    asol = AtomicSolution(dae)

    # test constructors and general functionality
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

    # test initial conditions
    set_initial_conditions!(sol, t0, z0, λ0)
    get_initial_conditions!(sol, tz, tλ, 1)
    get_initial_conditions!(sol, asol, 1)

    @test tz == z0
    @test tλ == λ0
    @test asol.t == t0
    @test asol.q == z0
    @test asol.λ == λ0

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

    # test set/get solution
    sol1 = Solution(similar(dae, z0, λ0), Δt, nt)
    sol2 = Solution(similar(dae, z0, λ0), Δt, nt)
    for i in 1:nt
        tz .= zs[:,i]
        tλ .= λs[:,i]
        asol.q .= zs[:,i]
        asol.λ .= λs[:,i]
        set_solution!(sol1, tz, tλ, i)
        set_solution!(sol2, asol, i)
    end
    @test sol1.q[:,1:nt] == zs
    @test sol1.λ[:,1:nt] == λs
    @test sol2.q[:,1:nt] == zs
    @test sol2.λ[:,1:nt] == λs

    sol1 = Solution(similar(dae, z1, λ1), Δt, nt)
    sol2 = Solution(similar(dae, z1, λ1), Δt, nt)
    for i in 1:nt
        for k in 1:ni
            tz .= Zs[:,i,k]
            tλ .= Λs[:,i,k]
            asol.q .= Zs[:,i,k]
            asol.λ .= Λs[:,i,k]
            set_solution!(sol1, tz, tλ, i, k)
            set_solution!(sol2, asol, i, k)
        end
    end
    @test sol1.q[:,1:nt,:] == Zs
    @test sol1.λ[:,1:nt,:] == Λs
    @test sol2.q[:,1:nt,:] == Zs
    @test sol2.λ[:,1:nt,:] == Λs

    # test hdf5 in- and output
    sol1 = Solution(harmonic_oscillator_dae(z0, λ0), Δt, nt)
    create_hdf5!(sol1, h5file)
    write_to_hdf5(sol1)
    close(sol1)
    @test isfile(h5file)
    sol2 = SSolutionDAE(h5file)
    @test sol1   != sol2
    @test sol1.t == sol2.t
    @test sol1.q == sol2.q
    @test sol1.λ == sol2.λ
    @test sol1.ntime == sol2.ntime
    @test sol1.nsave == sol2.nsave
    rm(h5file)

    sol1 = Solution(harmonic_oscillator_dae(z1, λ1), Δt, nt)
    create_hdf5!(sol1, h5file)
    write_to_hdf5(sol1)
    close(sol1)
    @test isfile(h5file)
    sol2 = SSolutionDAE(h5file)
    @test sol1   != sol2
    @test sol1.t == sol2.t
    @test sol1.q == sol2.q
    @test sol1.λ == sol2.λ
    @test sol1.ntime == sol2.ntime
    @test sol1.nsave == sol2.nsave
    rm(h5file)

    # test nsave and nwrite parameters
    sol = Solution(dae, Δt, 20, nsave=2)
    @test sol.nt == 10

    sol = Solution(dae, Δt, 20, nsave=2, nwrite=10)
    @test sol.nt == 5

    # test reset
    sol = Solution(dae, Δt, nt)
    reset!(sol)
    @test sol.t[0]   == t1
    @test sol.t[end] == t2
    @test offset(sol) == nt

    sol1 = Solution(dae, Δt, 100, nwrite=10)
    create_hdf5!(sol1, h5file)
    for i in 1:10
        for j in 1:10
            oj = (i-1)*10+j
            set_solution!(sol1, z100[:,oj], λ100[:,oj], j)
        end
        write_to_hdf5(sol1)
        reset!(sol1)
    end
    close(sol1)
    sol2 = SSolutionDAE(h5file)
    @test sol2.t.t ≈ Δt .* collect(0:100) atol=1E-14
    @test sol2.q.d == hcat(reshape(dae.q₀, (ndims(dae),1)), z100)
    @test sol2.λ.d == hcat(reshape(dae.λ₀, (dae.m,1)), λ100)
    @test sol2.ntime == 100
    @test sol2.nsave == 1
    rm(h5file)

    sol1 = Solution(dae, Δt, 100, nsave=2, nwrite=10)
    create_hdf5!(sol1, h5file)
    for i in 1:10
        for j in 1:10
            oj = (i-1)*10+j
            set_solution!(sol1, z100[:,oj], λ100[:,oj], j)
        end
        write_to_hdf5(sol1)
        reset!(sol1)
    end
    close(sol1)
    sol2 = SSolutionDAE(h5file)
    @test sol2.t.t ≈ Δt .* collect(0:2:100) atol=1E-14
    @test sol2.q.d == hcat(reshape(dae.q₀, (ndims(dae),1)), z100)[:,1:2:end]
    @test sol2.λ.d == hcat(reshape(dae.λ₀, (dae.m,1)), λ100)[:,1:2:end]
    @test sol2.ntime == 100
    @test sol2.nsave == 2
    rm(h5file)
end


@testset "$(rpad("PDAE Solution",80))" begin
    asol = AtomicSolution(pdae)

    # test constructors and general functionality
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

    # test initial conditions
    set_initial_conditions!(sol, t0, x0, y0, μ0)
    get_initial_conditions!(sol, tx, ty, tμ, 1)
    get_initial_conditions!(sol, asol, 1)

    @test tx == x0
    @test ty == y0
    @test tμ == μ0
    @test asol.t == t0
    @test asol.q == x0
    @test asol.p == y0
    @test asol.λ == μ0

    tx .= 0
    ty .= 0
    get_initial_conditions!(sol, tx, ty, 1)
    @test tx == x0
    @test ty == y0

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
        @test tμ == μ2[:,i]
    end

    # test set/get solution
    sol1 = Solution(similar(pdae, x0, y0, μ0), Δt, nt)
    sol2 = Solution(similar(pdae, x0, y0, μ0), Δt, nt)
    sol3 = Solution(similar(pdae, x0, y0, μ0), Δt, nt)
    for i in 1:nt
        tx .= xs[:,i]
        ty .= ys[:,i]
        tμ .= μs[:,i]
        asol.q .= xs[:,i]
        asol.p .= ys[:,i]
        asol.λ .= μs[:,i]
        set_solution!(sol1, tx, ty, tμ, i)
        set_solution!(sol2, tx, ty, i)
        set_solution!(sol3, asol, i)
    end
    @test sol1.q[:,1:nt] == xs
    @test sol1.p[:,1:nt] == ys
    @test sol1.λ[:,1:nt] == μs
    @test sol2.q[:,1:nt] == xs
    @test sol2.p[:,1:nt] == ys
    @test sol3.q[:,1:nt] == xs
    @test sol3.p[:,1:nt] == ys
    @test sol3.λ[:,1:nt] == μs

    sol1 = Solution(similar(pdae, x1, y1, μ1), Δt, nt)
    sol2 = Solution(similar(pdae, x1, y1, μ1), Δt, nt)
    sol3 = Solution(similar(pdae, x1, y1, μ1), Δt, nt)
    for i in 1:nt
        for k in 1:ni
            tx .= Xs[:,i,k]
            ty .= Ys[:,i,k]
            tμ .= Ms[:,i,k]
            asol.q .= Xs[:,i,k]
            asol.p .= Ys[:,i,k]
            asol.λ .= Ms[:,i,k]
            set_solution!(sol1, tx, ty, tμ, i, k)
            set_solution!(sol2, tx, ty, i, k)
            set_solution!(sol3, asol, i, k)
        end
    end
    @test sol1.q[:,1:nt,:] == Xs
    @test sol1.p[:,1:nt,:] == Ys
    @test sol1.λ[:,1:nt,:] == Ms
    @test sol2.q[:,1:nt,:] == Xs
    @test sol2.p[:,1:nt,:] == Ys
    @test sol3.q[:,1:nt,:] == Xs
    @test sol3.p[:,1:nt,:] == Ys
    @test sol3.λ[:,1:nt,:] == Ms

    # test hdf5 in- and output
    sol1 = Solution(harmonic_oscillator_pdae(x0, y0, μ0), Δt, nt)
    create_hdf5!(sol1, h5file)
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
    create_hdf5!(sol1, h5file)
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

    # test nsave and nwrite parameters
    sol = Solution(pdae, Δt, 20, nsave=2)
    @test sol.nt == 10

    sol = Solution(pdae, Δt, 20, nsave=2, nwrite=10)
    @test sol.nt == 5

    # test reset
    sol = Solution(pdae, Δt, nt)
    reset!(sol)
    @test sol.t[0]   == t1
    @test sol.t[end] == t2
    @test offset(sol) == nt

    sol1 = Solution(pdae, Δt, 100, nwrite=10)
    create_hdf5!(sol1, h5file)
    for i in 1:10
        for j in 1:10
            oj = (i-1)*10+j
            set_solution!(sol1, x100[:,oj], y100[:,oj], μ100[:,oj], j)
        end
        write_to_hdf5(sol1)
        reset!(sol1)
    end
    close(sol1)
    sol2 = SSolutionPDAE(h5file)
    @test sol2.t.t ≈ Δt .* collect(0:100) atol=1E-14
    @test sol2.q.d == hcat(reshape(pdae.q₀, (ndims(pdae),1)), x100)
    @test sol2.p.d == hcat(reshape(pdae.p₀, (ndims(pdae),1)), y100)
    @test sol2.λ.d == hcat(reshape(pdae.λ₀, (pdae.m,1)), μ100)
    @test sol2.ntime == 100
    @test sol2.nsave == 1
    rm(h5file)

    sol1 = Solution(pdae, Δt, 100, nsave=2, nwrite=10)
    create_hdf5!(sol1, h5file)
    for i in 1:10
        for j in 1:10
            oj = (i-1)*10+j
            set_solution!(sol1, x100[:,oj], y100[:,oj], μ100[:,oj], j)
        end
        write_to_hdf5(sol1)
        reset!(sol1)
    end
    close(sol1)
    sol2 = SSolutionPDAE(h5file)
    @test sol2.t.t ≈ Δt .* collect(0:2:100) atol=1E-14
    @test sol2.q.d == hcat(reshape(pdae.q₀, (ndims(pdae),1)), x100)[:,1:2:end]
    @test sol2.p.d == hcat(reshape(pdae.p₀, (ndims(pdae),1)), y100)[:,1:2:end]
    @test sol2.λ.d == hcat(reshape(pdae.λ₀, (pdae.m,1)), μ100)[:,1:2:end]
    @test sol2.ntime == 100
    @test sol2.nsave == 2
    rm(h5file)
end
