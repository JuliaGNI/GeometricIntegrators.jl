
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
Pages   = ["integrators/splitting/integrators_composition.jl",
           "integrators/splitting/integrators_splitting.jl",
           "integrators/splitting/integrators_exact_ode.jl",
           "integrators/splitting/splitting_tableau.jl"]
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
Modules = [GeometricIntegrators.Integrators.SPARK]
```


## VPRK Methods

```@autodocs
Modules = [GeometricIntegrators.Integrators.VPRK]
```


## Galerkin Variational Integrators

```@autodocs
Modules = [GeometricIntegrators.Integrators]
Pages   = ["integrators/cgvi/integrators_cgvi.jl",
           "integrators/dgvi/integrators_dgvi.jl"]
```
