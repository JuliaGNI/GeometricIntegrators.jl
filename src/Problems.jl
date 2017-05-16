__precompile__()

module Problems

    export ChargedParticle3dSingular, ChargedParticle3dSymmetric, ChargedParticle3dUniform,
           ExponentialGrowth, LotkaVolterra2d, Pendulum, PointVortices

    include("problems/charged_particle_3d_singular.jl")
    include("problems/charged_particle_3d_symmetric.jl")
    include("problems/charged_particle_3d_uniform.jl")
    include("problems/exponential_growth.jl")
    include("problems/lotka_volterra_2d.jl")
    include("problems/pendulum.jl")
    include("problems/point_vortices.jl")

end
