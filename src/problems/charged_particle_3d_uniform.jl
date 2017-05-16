module ChargedParticle3dUniform

    export charged_particle_3d_iode, hamiltonian, angular_momentum

    const A₀ = 1
    const E₀ = 2
    const q₀ = [1., 0., 0., 0., 1., 1.]


    function A1(t,q)
       return - 0.5 * A₀ * q[2]
    end

    function A2(t,q)
       return + 0.5 * A₀ * q[1]
    end

    function A3(t,q)
       return zero(eltype(q))
    end


    function dA1d1(t,q)
       return zero(eltype(q))
    end

    function dA1d2(t,q)
       return - 0.5 * A₀
    end

    function dA1d3(t,q)
       return zero(eltype(q))
    end


    function dA2d1(t,q)
       return + 0.5 * A₀
    end

    function dA2d2(t,q)
       return zero(eltype(q))
    end

    function dA2d3(t,q)
       return zero(eltype(q))
    end


    function dA3d1(t,q)
       return zero(eltype(q))
    end

    function dA3d2(t,q)
       return zero(eltype(q))
    end

    function dA3d3(t,q)
       return zero(eltype(q))
    end


    function B(t,q)
       return A₀
    end


    function dBd1(t,q)
       return zero(eltype(q))
    end

    function dBd2(t,q)
       return zero(eltype(q))
    end

    function dBd3(t,q)
       return zero(eltype(q))
    end


    function B1(t,q)
       return zero(eltype(q))
    end

    function B2(t,q)
       return zero(eltype(q))
    end

    function B3(t,q)
       return B(t,q)
    end


    function b1(t,q)
       return zero(eltype(q))
    end

    function b2(t,q)
       return zero(eltype(q))
    end

    function b3(t,q)
       return one(eltype(q))
    end


    function db1d1(t,q)
       return zero(eltype(q))
    end

    function db1d2(t,q)
       return zero(eltype(q))
    end

    function db1d3(t,q)
       return zero(eltype(q))
    end


    function db2d1(t,q)
       return zero(eltype(q))
    end

    function db2d2(t,q)
       return zero(eltype(q))
    end

    function db2d3(t,q)
       return zero(eltype(q))
    end


    function db3d1(t,q)
       return zero(eltype(q))
    end

    function db3d2(t,q)
       return zero(eltype(q))
    end

    function db3d3(t,q)
       return zero(eltype(q))
    end


    function R(t,x)
        return one(eltype(x))
    end

    function dRd1(t,x)
        return zero(eltype(x))
    end

    function dRd2(t,x)
        return zero(eltype(x))
    end

    function dRd3(t,x)
        return zero(eltype(x))
    end


    include("charged_particle_3d.jl")

end
