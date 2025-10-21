# Integrators


## Geometric Integrator

```@autodocs
Modules = [GeometricIntegrators.Integrators]
Pages   = [
            "integrators/integrator_cache.jl",
            "integrators/integrator.jl",
        ]
```

## Initial Guesses

```@autodocs
Modules = [GeometricIntegrators.Integrators]
Pages   = [
            "initial_guess/initial_guess.jl",
            "initial_guess/hermite.jl",
            "initial_guess/midpoint.jl",
          ]
```

## Runge-Kutta Integrators

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
            "integrators/vi/deleqs.jl",
            "integrators/vi/deleqs_methods.jl",
            "integrators/vi/vi_methods.jl",
            "integrators/vi/position_momentum_midpoint.jl",
            "integrators/vi/position_momentum_trapezoidal.jl",
            "integrators/vi/vprk_integrator.jl",
            "integrators/vi/vprk_methods.jl",
        ]
```

## Degenerate Variational Integrators

```@autodocs
Modules = [GeometricIntegrators.Integrators]
Pages   = [
            "integrators/dvi/dvi_euler.jl",
            "integrators/dvi/dvi_midpoint.jl",
            "integrators/dvi/dvi_trapezoidal.jl",
            "integrators/dvi/dvrk.jl",
        ]
```

## Galerkin Variational Integrators

```@autodocs
Modules = [GeometricIntegrators.Integrators]
Pages   = ["integrators/cgvi/integrators_cgvi.jl",
           "integrators/dgvi/integrators_dgvi.jl"]
```

## Hamilton-Pontryagin Integrators

```@autodocs
Modules = [GeometricIntegrators.Integrators]
Pages   = [
            "integrators/hpi/hpi_methods.jl",
            "integrators/hpi/hpi_midpoint.jl",
            "integrators/hpi/hpi_trapezoidal.jl",
        ]
```

## Splitting and Composition Methods

```@autodocs
Modules = [GeometricIntegrators.Integrators]
Pages   = [
            "integrators/splitting/exact_solution.jl",
            "integrators/splitting/splitting_coefficients.jl",
            "integrators/splitting/splitting_methods.jl",
            "integrators/splitting/splitting_integrator.jl",
            "integrators/splitting/composition_integrator.jl",
            "integrators/splitting/composition_methods.jl",
        ]
```
