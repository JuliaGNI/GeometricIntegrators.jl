__precompile__()

module Problems

    export ChargedParticle3dSingular, ChargedParticle3dSymmetric, ChargedParticle3dUniform,
           ExponentialGrowth, LotkaVolterra2d, Pendulum,
           GuidingCenter4dTokamakPassing, GuidingCenter4dTokamakTrapped,
           GuidingCenter4dTokamakBarelyTrapped


    include("problems/charged_particle_3d_singular.jl")
    include("problems/charged_particle_3d_symmetric.jl")
    include("problems/charged_particle_3d_uniform.jl")
    include("problems/exponential_growth.jl")
    include("problems/guiding_center_4d_tokamak_barely_trapped.jl")
    include("problems/guiding_center_4d_tokamak_passing.jl")
    include("problems/guiding_center_4d_tokamak_trapped.jl")
    include("problems/pendulum.jl")
    include("problems/lotka_volterra_2d.jl")



end
