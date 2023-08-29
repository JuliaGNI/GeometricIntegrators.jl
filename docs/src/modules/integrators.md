
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
Pages   = [
            "integrators/initial_guess/initial_guess.jl",
            "integrators/initial_guess/hermite.jl",
            "integrators/initial_guess/midpoint.jl",
          ]
```


## Extrapolation Methods

The extrapolation routines are exclusively used for computing
initial guesses and are usually not called directly by the user.

```@autodocs
Modules = [GeometricIntegrators.Integrators]
Order   = [:constant, :type, :macro, :function]
Pages   = ["integrators/extrapolation/extrapolation.jl",
           "integrators/extrapolation/aitken_neville.jl",
           "integrators/extrapolation/euler.jl",
           "integrators/extrapolation/hermite.jl",
           "integrators/extrapolation/midpoint.jl"]
```


# Euler Integrators

```@docs
GeometricIntegrators.Integrators.ExplicitEuler
GeometricIntegrators.Integrators.ImplicitEuler
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
Pages   = ["integrators/rk/abstract_integrator_rk.jl",
           "integrators/rk/coefficients.jl",
           "integrators/rk/tableaus.jl",
           "integrators/rk/integrators_erk.jl",
           "integrators/rk/integrators_irk.jl",
           "integrators/rk/integrators_irk_implicit.jl",
           "integrators/rk/integrators_dirk.jl",
           "integrators/rk/integrators_eprk.jl",
           "integrators/rk/integrators_iprk.jl",
           "integrators/rk/integrators_iprk_implicit.jl",
          ]
```


## Variational Integrators

```@autodocs
Modules = [GeometricIntegrators.Integrators]
Pages   = ["integrators/vi/integrators_vprk.jl"]
```

## Degenerate Variational Integrators

```@autodocs
Modules = [GeometricIntegrators.Integrators]
Pages   = ["integrators/dvi/integrators_dvi_a.jl",
           "integrators/dvi/integrators_dvi_b.jl",
           "integrators/dvi/integrators_cmdvi.jl",
           "integrators/dvi/integrators_ctdvi.jl",
           "integrators/dvi/integrators_dvrk.jl"]
```

## Galerkin Variational Integrators

```@autodocs
Modules = [GeometricIntegrators.Integrators]
Pages   = ["integrators/cgvi/integrators_cgvi.jl",
           "integrators/dgvi/integrators_dgvi.jl"]
```
