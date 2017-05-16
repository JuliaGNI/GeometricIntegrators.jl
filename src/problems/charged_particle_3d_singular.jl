module ChargedParticle3dSingular

    export charged_particle_3d_iode, hamiltonian, angular_momentum

    const E₀ = 0
    const q₀ = [1., 0., 0., 0., -1., 0.]


    function denominator(x,y)
        sqrt(x^2 + y^2)
        # sqrt(sqrt(eps()) + x^2 + y^2)
        # sqrt(1 + x^2 + y^2)
    end


    function A1(t,q)
       return + q[2] / denominator(q[1], q[2])^3
    end

    function A2(t,q)
       return - q[1] / denominator(q[1], q[2])^3
    end

    function A3(t,q)
       return zero(eltype(q))
    end


    function dA1d1(t,q)
       return - 3 * q[1] * q[2] / denominator(q[1], q[2])^5
    end

    function dA1d2(t,q)
       return - 3 * q[2] * q[2] / denominator(q[1], q[2])^5 + 1 / denominator(q[1], q[2])^3
    end

    function dA1d3(t,q)
       return zero(eltype(q))
    end


    function dA2d1(t,q)
       return + 3 * q[1] * q[1] / denominator(q[1], q[2])^5 - 1 / denominator(q[1], q[2])^3
    end

    function dA2d2(t,q)
       return + 3 * q[1] * q[2] / denominator(q[1], q[2])^5
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


    function B₀(x,y)
       return 1 / denominator(x, y)^3
    end

    function B1(t,q)
       return zero(eltype(q))
    end

    function B2(t,q)
       return zero(eltype(q))
    end

    function B3(t,q)
       return B₀(q[1], q[2])
    end



    function B(t,q)
       return B₀(q[1], q[2])
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


    include("charged_particle_3d.jl")

end
