__precompile__()

module Problems

    export ChargedParticle3dSingular, ChargedParticle3dSymmetric,
           ChargedParticle3dUniform, GuidingCenter4dTokamakBarelyTrapped,
           GuidingCenter4dTokamakPassing, GuidingCenter4dTokamakTrapped,
           ExponentialGrowth, LotkaVolterra2d, Pendulum, PointVortices

    include("problems/charged_particle_3d_singular.jl")
    include("problems/charged_particle_3d_symmetric.jl")
    include("problems/charged_particle_3d_uniform.jl")
    include("problems/guiding_center_4d_tokamak_barely_trapped.jl")
    include("problems/guiding_center_4d_tokamak_passing.jl")
    include("problems/guiding_center_4d_tokamak_trapped.jl")
    include("problems/exponential_growth.jl")
    include("problems/lotka_volterra_2d.jl")
    include("problems/pendulum.jl")
    include("problems/point_vortices.jl")

end
