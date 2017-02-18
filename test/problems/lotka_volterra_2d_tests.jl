
module LotkaVolterra2dTests

    using GeometricIntegrators
    using GeometricIntegrators.Problems.LotkaVolterra2d

    const Δt = 0.1
    const nt = 1

    export test_lotka_volterra_2d

    function test_lotka_volterra_2d()
        ode  = lotka_volterra_2d_ode()
        iode = lotka_volterra_2d_iode()
        idae = lotka_volterra_2d_idae()

        int = Integrator(ode, getTableauGLRK(1), Δt)
        sol = integrate(int, nt)

        vint = IntegratorVPRKpStandard(iode, getTableauVPGLRK(1), Δt)
        isol = integrate(vint, nt)

        dint = Integrator(idae, getTableauGLRKpSymplectic(1), Δt)
        dsol = integrate(dint, nt)

        dint = Integrator(idae, getTableauGLRKpSymmetric(1), Δt)
        dsol = integrate(dint, nt)
    end

end


using LotkaVolterra2dTests

test_lotka_volterra_2d()
