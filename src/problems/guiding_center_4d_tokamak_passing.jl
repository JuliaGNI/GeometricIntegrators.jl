module GuidingCenter4dTokamakPassing

    export guiding_center_4d_ode, guiding_center_4d_iode,
           hamiltonian, toroidal_momentum, α, α1, α2, α3, α4, β1, β2, β3, b1, b2, b3

    const μ  = 2.448E-6
    const q₀ = [1.05, 0., 0., 8.117E-04]

    include("magnetic_field_tokamak.jl")
    include("guiding_center_4d.jl")

end
