module GuidingCenter4dTokamakTrapped

    export guiding_center_4d_ode, guiding_center_4d_iode,
           hamiltonian, toroidal_momentum, u, α, α1, α2, α3, α4, dα, β, β1, β2, β3, b1, b2, b3, dH

    const ν  = 1.
    const μ  = 2.250E-6
    const q₀ = [1.05, 0., 0., 4.306E-4]

    include("magnetic_field_tokamak.jl")
    include("guiding_center_4d.jl")

end
