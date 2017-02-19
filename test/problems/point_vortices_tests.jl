
module PointVorticesTests

    using GeometricIntegrators
    using GeometricIntegrators.Problems.PointVortices

    const Δt = 0.1
    const nt = 1

    export test_point_vortices

    function test_point_vortices()
        int = Integrator(point_vortices_ode(), getTableauGLRK(1), Δt)
        sol = integrate(int, nt)

        vint = IntegratorVPRKpStandard(point_vortices_iode(), getTableauVPGLRK(1), Δt)
        isol = integrate(vint, nt)
    end

end


using PointVorticesTests

test_point_vortices()
