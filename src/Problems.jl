__precompile__()

module Problems

    export ChargedParticle3dSingular, ChargedParticle3dSymmetric,
           ChargedParticle3dUniform, GuidingCenter4dTokamakBarelyTrapped,
           GuidingCenter4dTokamakPassing, GuidingCenter4dTokamakTrapped,
           GuidingCenter4dTokamakLoop, GuidingCenter4dTokamakSurface,
           GuidingCenter4dSymmetricLoop, GuidingCenter4dSymmetricSurface,
           GuidingCenter4dUniformLoop,
           ExponentialGrowth, LotkaVolterra2d, Pendulum, PointVortices

    include("problems/charged_particle_3d_singular.jl")
    include("problems/charged_particle_3d_symmetric.jl")
    include("problems/charged_particle_3d_uniform.jl")
    include("problems/guiding_center_4d_tokamak_barely_trapped.jl")
    include("problems/guiding_center_4d_tokamak_passing.jl")
    include("problems/guiding_center_4d_tokamak_trapped.jl")
    include("problems/guiding_center_4d_tokamak_loop.jl")
    include("problems/guiding_center_4d_tokamak_surface.jl")
    include("problems/guiding_center_4d_uniform_loop.jl")
    include("problems/guiding_center_4d_symmetric_loop.jl")
    include("problems/guiding_center_4d_symmetric_surface.jl")
    include("problems/exponential_growth.jl")
    include("problems/lotka_volterra_2d.jl")
    include("problems/pendulum.jl")
    include("problems/point_vortices.jl")

end
