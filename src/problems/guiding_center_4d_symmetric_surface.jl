module GuidingCenter4dSymmetricSurface

    export guiding_center_4d_surface_ode, guiding_center_4d_surface_iode,
           hamiltonian, toroidal_momentum, α, α1, α2, α3, α4, β, β1, β2, β3, b1, b2, b3

    include("magnetic_field_symmetric.jl")

    const μ  = 2.5E-6


    function f_surface(s,t)
        φ0 = 0.0
        u0 = 4E-03
        r0 = 0.1

        Rt = 2r0*(s-0.5)
        Zt = 2r0*(t-0.5)

        qt = [Rt, Zt, φ0, u0]

        return qt
    end

    function f_surface(i, j, nx, ny)
        f_surface(i/nx, j/ny)
    end

    function get_initial_conditions(nx, ny)
        q₀ = zeros(4, nx*ny)

        for j in 1:ny
            for i in 1:nx
                q₀[:,nx*(j-1)+i] = f_surface(i, j, nx, ny)
            end
        end

        return q₀
    end


    function guiding_center_4d_surface_ode(nx, ny)
        q₀ = get_initial_conditions(nx, ny)
        guiding_center_4d_ode(q₀; periodic=false)
    end

    function guiding_center_4d_surface_iode(nx, ny)
        q₀ = get_initial_conditions(nx, ny)
        guiding_center_4d_iode(q₀; periodic=false)
    end


    include("guiding_center_4d.jl")

end
