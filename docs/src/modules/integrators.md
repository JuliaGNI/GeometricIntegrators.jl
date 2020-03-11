
# Integrators


## Common

```@autodocs
Modules = [GeometricIntegrators.Integrators]
Pages   = ["integrators/abstract_coefficients.jl",
           "integrators/abstract_integrator.jl",
           "integrators/abstract_tableau.jl",
           "integrators/integrator_cache.jl",
           "integrators/integrators_common.jl",
           "integrators/integrators.jl"]
```

## Initial Guesses

```@autodocs
Modules = [GeometricIntegrators.Integrators]
Pages   = ["integrators/initial_guess/extrapolation.jl",
           "integrators/initial_guess/initial_guess_ode.jl",
           "integrators/initial_guess/initial_guess_pode.jl"]
```


## Splitting Methods

```@autodocs
Modules = [GeometricIntegrators.Integrators]
Pages   = ["integrators/splitting/integrators_splitting.jl"]
```


## Runge-Kutta Methods

```@autodocs
Modules = [GeometricIntegrators.Integrators]
Pages   = ["integrators/rk/bstract_integrator_rk.jl",
           "integrators/rk/coefficients.jl",
           "integrators/rk/tableaus.jl",
           "integrators/rk/integrators_erk.jl",
           "integrators/rk/integrators_dirk.jl",
           "integrators/rk/integrators_firk.jl",
           "integrators/rk/integrators_eprk.jl",
           "integrators/rk/integrators_iprk.jl",
           "integrators/rk/integrators_flrk.jl",
           "integrators/rk/integrators_pglrk.jl"]
```


## SPARK Methods

```@autodocs
Modules = [GeometricIntegrators.Integrators]
Pages   = ["integrators/spark/abstract_integrator_spark.jl",
           "integrators/spark/coefficients.jl",
           "integrators/spark/integrators_hpark.jl",
           "integrators/spark/integrators_hspark_common.jl",
           "integrators/spark/integrators_hspark.jl",
           "integrators/spark/integrators_hspark_primary.jl",
           "integrators/spark/integrators_hspark_secondary.jl",
           "integrators/spark/integrators_spark_cache.jl",
           "integrators/spark/integrators_spark_common.jl",
           "integrators/spark/integrators_spark_tableau.jl",
           "integrators/spark/integrators_vpark.jl",
           "integrators/spark/integrators_vspark_common.jl",
           "integrators/spark/integrators_vspark.jl",
           "integrators/spark/integrators_vspark_primary.jl",
           "integrators/spark/integrators_vspark_secondary.jl",
           "integrators/spark/integrators_slrk.jl"]
```


## Galerkin Variational Integrators

```@autodocs
Modules = [GeometricIntegrators.Integrators]
Pages   = ["integrators/cgvi/integrators_cgvi.jl",
           "integrators/dgvi/integrators_dgvi.jl"]
```
