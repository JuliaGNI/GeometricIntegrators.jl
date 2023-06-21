```@meta
CurrentModule = GeometricIntegrators.SPARK
```

# Special Partitioned Additive Runge-Kutta Integrators

SPARK or *Special Partitioned Additive Runge-Kutta Integrators* are a family of integrators that have been introduced by [Laurent O. Jay](http://homepage.divms.uiowa.edu/~ljay/) for the integration of differential algebraic equations and in particular systems subject to holonomic and nonholonomic constraints [[Jay:1998](@cite), [Jay:2003](@cite), [Jay:2006](@cite), [Jay:2007a](@cite), [Jay:2007b](@cite)].
Recently, the idea of SPARK methods has been generalized and adapted to facilitate the integration of degenerate Lagrangian systems as well as Hamiltonian systems subject to Dirac constraints [[Kraus:2020](@cite)].

GeometricIntegrators.jl provides several flavours of such SPARK methods (*some are still experimental*):

| Integrator                          | Description                                                                                          |
|:------------------------------------|:-----------------------------------------------------------------------------------------------------|
| [`IntegratorHPARK`](@ref)           | Partitioned additive methods for Hamiltonian system subject to a general constraint $\phi(q,p) = 0$  |
| [`IntegratorVPARK`](@ref)           | Partitioned additive methods for Lagrangian system subject to a general constraint $\phi(q,p) = 0$   |
| [`IntegratorSPARK`](@ref)           | SPARK methods for general index-two differential algebraic equations                                 |
| [`IntegratorHSPARK`](@ref)          | Hamiltonian system subject to a general constraint $\phi(q,p) = 0$                                   |
| [`IntegratorHSPARKprimary`](@ref)   | Hamiltonian system subject primary constraint in the sense of Dirac                                  |
| [`IntegratorHSPARKsecondary`](@ref) | Hamiltonian system enforcing primary & secondary Dirac constraint                                    |
| [`IntegratorVSPARK`](@ref)          | Lagrangian system in implicit form subject to a general constraint $\phi(q,p) = 0$                   |
| [`IntegratorVSPARKprimary`](@ref)   | Degenerate Lagrangian system subject primary constraint in the sense of Dirac                        |
| [`IntegratorVSPARKsecondary`](@ref) | Degenerate Lagrangian system enforcing primary & secondary Dirac constraint                          |

These integrators are applied to either an `IDAE`, `HDAE` or `LDAE`.
