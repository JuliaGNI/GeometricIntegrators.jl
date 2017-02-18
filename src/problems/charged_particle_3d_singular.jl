module ChargedParticle3dSingular

    export charged_particle_3d_iode, hamiltonian, angular_momentum

    const E₀ = 0
    const q₀ = [1., 0., 0., 0., -1., 0.]

    include("magnetic_field_singular.jl")
    include("charged_particle_3d.jl")

end
