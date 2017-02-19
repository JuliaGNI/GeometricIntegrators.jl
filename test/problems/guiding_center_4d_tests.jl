

module GuidingCenter4dTests

    using GeometricIntegrators

    const Δt = 400.
    const nt = 1

    export test_guiding_center_4d_glrk, test_guiding_center_4d_vpglrk

    function test_guiding_center_4d_glrk(ode)
        int = Integrator(ode, getTableauGLRK(1), Δt)
        sol = integrate(int, nt)
    end

    function test_guiding_center_4d_vpglrk(ode)
        int = Integrator(ode, getTableauVPGLRK(1), Δt)
        sol = integrate(int, nt)
    end

end


using GuidingCenter4dTests

import GeometricIntegrators.Problems.GuidingCenter4dTokamakPassing
import GeometricIntegrators.Problems.GuidingCenter4dTokamakTrapped
import GeometricIntegrators.Problems.GuidingCenter4dTokamakBarelyTrapped

test_guiding_center_4d_glrk(GeometricIntegrators.Problems.GuidingCenter4dTokamakPassing.guiding_center_4d_ode())
test_guiding_center_4d_glrk(GeometricIntegrators.Problems.GuidingCenter4dTokamakTrapped.guiding_center_4d_ode())
test_guiding_center_4d_glrk(GeometricIntegrators.Problems.GuidingCenter4dTokamakBarelyTrapped.guiding_center_4d_ode())

test_guiding_center_4d_vpglrk(GeometricIntegrators.Problems.GuidingCenter4dTokamakPassing.guiding_center_4d_iode())
test_guiding_center_4d_vpglrk(GeometricIntegrators.Problems.GuidingCenter4dTokamakTrapped.guiding_center_4d_iode())
test_guiding_center_4d_vpglrk(GeometricIntegrators.Problems.GuidingCenter4dTokamakBarelyTrapped.guiding_center_4d_iode())
