
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

```@autodocs
Modules = [GeometricIntegrators.Integrators]
Pages   = [
            "integrators/euler/explicit_euler.jl",
            "integrators/euler/implicit_euler.jl",
        ]
```


## Runge-Kutta Methods

```@autodocs
Modules = [GeometricIntegrators.Integrators]
Pages   = ["integrators/rk/abstract.jl",
           "integrators/rk/common.jl",
           "integrators/rk/tableaus.jl",
           "integrators/rk/updates.jl",
           "integrators/rk/integrators_erk.jl",
           "integrators/rk/integrators_irk.jl",
           "integrators/rk/integrators_irk_implicit.jl",
           "integrators/rk/integrators_dirk.jl",
           "integrators/rk/integrators_eprk.jl",
           "integrators/rk/integrators_iprk.jl",
           "integrators/rk/integrators_iprk_implicit.jl",
           "integrators/rk/methods.jl",
          ]
```


## Variational Integrators

```@autodocs
Modules = [GeometricIntegrators.Integrators]
Pages   = [
            "integrators/vi/integrators_vprk.jl",
            "integrators/vi/vprk_methods.jl",
        ]
```

## Degenerate Variational Integrators

```@autodocs
Modules = [GeometricIntegrators.Integrators]
Pages   = [
            "integrators/dvi/integrators_dvi_a.jl",
            "integrators/dvi/integrators_dvi_b.jl",
            "integrators/dvi/integrators_cmdvi.jl",
            "integrators/dvi/integrators_ctdvi.jl",
            "integrators/dvi/integrators_dvrk.jl",
            "integrators/dvi/methods.jl",
        ]
```

## Galerkin Variational Integrators

```@autodocs
Modules = [GeometricIntegrators.Integrators]
Pages   = ["integrators/cgvi/integrators_cgvi.jl",
           "integrators/dgvi/integrators_dgvi.jl"]
```


## Splitting Methods

```@autodocs
Modules = [GeometricIntegrators.Integrators]
Pages   = [
            "integrators/splitting/exact_solution.jl",
            "integrators/splitting/composition_integrators.jl",
            "integrators/splitting/composition_methods.jl",
            "integrators/splitting/splitting_coefficients.jl",
            "integrators/splitting/splitting_integrator.jl",
            "integrators/splitting/splitting_methods.jl",
        ]
```
