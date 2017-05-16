module ChargedParticle3dSymmetric

    export charged_particle_3d_iode, hamiltonian, angular_momentum

    const E₀ = 2
    const q₀ = [1., 0., 0., 0., 1., 1.]


    function A1(t,x)
       return - x[2] * (2 + x[1]^2 + x[2]^2) / 4
    end

    function A2(t,x)
       return + x[1] * (2 + x[1]^2 + x[2]^2) / 4
    end

    function A3(t,x)
       return zero(eltype(x))
    end


    function dA1d1(t,x)
       return - x[1] * x[2] / 2
    end

    function dA1d2(t,x)
       return - (2 + x[1]^2 + 3 * x[2]^2) / 4
    end

    function dA1d3(t,x)
       return zero(eltype(x))
    end


    function dA2d1(t,x)
       return + (2 + 3 * x[1]^2 + x[2]^2) / 4
    end

    function dA2d2(t,x)
       return + x[1] * x[2] / 2
    end

    function dA2d3(t,x)
       return zero(eltype(x))
    end


    function dA3d1(t,x)
       return zero(eltype(x))
    end

    function dA3d2(t,x)
       return zero(eltype(x))
    end

    function dA3d3(t,x)
       return zero(eltype(x))
    end


    function B1(t,x)
        zero(eltype(x))
    end

    function B2(t,x)
        zero(eltype(x))
    end

    function B3(t,x)
        x[1]^2 + x[2]^2 + 1
    end


    function B(t,x)
        x[1]^2 + x[2]^2 + 1
    end


    function b1(t,x)
        zero(eltype(x))
    end

    function b2(t,x)
        zero(eltype(x))
    end

    function b3(t,x)
        one(eltype(x))
    end


    function R(t,x)
        return one(eltype(x))
    end

    function Z(t,x)
        return zero(eltype(x))
    end

    function phi(t,x)
        return zero(eltype(x))
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


    function db3d2(t,x)
        zero(eltype(x))
    end

    function db3d1(t,x)
        zero(eltype(x))
    end

    function dB2d1(t,x)
        zero(eltype(x))
    end

    function dB1d2(t,x)
        zero(eltype(x))
    end

    function dB1d3(t,x)
        zero(eltype(x))
    end

    function db2d2(t,x)
        zero(eltype(x))
    end

    function dB2d3(t,x)
        zero(eltype(x))
    end

    function dBd1(t,x)
        2*x[1]
    end

    function dB1d1(t,x)
        zero(eltype(x))
    end

    function db3d3(t,x)
        zero(eltype(x))
    end

    function dBd2(t,x)
        2*x[2]
    end

    function dB2d2(t,x)
        zero(eltype(x))
    end

    function dB3d3(t,x)
        zero(eltype(x))
    end

    function db1d1(t,x)
        zero(eltype(x))
    end

    function dB3d1(t,x)
        2*x[1]
    end

    function db1d2(t,x)
        zero(eltype(x))
    end

    function db2d3(t,x)
        zero(eltype(x))
    end

    function db1d3(t,x)
        zero(eltype(x))
    end

    function db2d1(t,x)
        zero(eltype(x))
    end

    function dB3d2(t,x)
        2*x[2]
    end

    function dBd3(t,x)
        zero(eltype(x))
    end


    include("charged_particle_3d.jl")

end
