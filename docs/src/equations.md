# Equation and Problem Types

In *GeometricIntegrators.jl* we support three basic types of equations:

* ordinary differential equations (ODEs),
* differential algebraic equations (DAEs),
* stochastic differential equations (SDEs).

For each type, there are several subtypes:

* standard equations ([`ODEProblem`](@ref), [`DAEProblem`](@ref), [`SDEProblem`](@ref)),
* implicit equations ([`IODEProblem`](@ref), [`IDAEProblem`](@ref)),
* partitioned equations ([`PODEProblem`](@ref), [`PDAEProblem`](@ref), [`PSDEProblem`](@ref)),
* Hamiltonian equations ([`HODEProblem`](@ref), [`HDAEProblem`](@ref)),
* Lagrangian equations ([`LODEProblem`](@ref), [`LDAEProblem`](@ref)),
* split equations ([`SODEProblem`](@ref), [`SPDAEProblem`](@ref), [`SPSDEProblem`](@ref)).

Each equation holds a number of functions determining the vector field, constraints, and possibly additional information like parameters, periodicity, invariants and the Hamiltonian or Lagrangian.

To each equation type there exists a corresponding problem type, which holds the equation, initial conditions, parameters, a time span to integrate over, as well as a time step (which is typically fixed in GeometricIntegrators). In addition, these problem types provide some convenience constructors to generate the equation and the problem at once.

All of the equation and problem types are defined in the [GeometricEquations.jl](https://github.com/JuliaGNI/GeometricEquations.jl) package.
The [GeometricProblems.jl](https://github.com/JuliaGNI/GeometricProblems.jl) package implements a number of example problems for testing and benchmarking.

## Keyword Arguments

All equation and problem types take the following keyword arguments:

* `invariants = NullInvariants()`
* `parameters = NullParameters()`
* `periodicity = NullPeriodicity()`

If not set to their corresponding Null types, the user needs to pass a `NamedTuple` whose values are

* functions for invariants,
* arbitrary data structures for parameters, 
* the same data structure as the solution for periodicity.

The latter should be zero everywhere, except for those components, that are periodic, i.e.,
whose value are supposed to stay within a range `(0, max)`. Support for ranges starting
with other values than zero is currently missing but can be added if demand arises.


## Ordinary Differential Equations (ODEs)

```@eval
using GeometricEquations, Markdown
Markdown.parse(GeometricEquations.ode_equations)
```

```@eval
using GeometricEquations, Markdown
Markdown.parse(GeometricEquations.ode_example)
```


### Partitioned Ordinary Differential Equations (PODEs)

```@eval
using GeometricEquations, Markdown
Markdown.parse(GeometricEquations.pode_equations)
```


### Hamiltonian Ordinary Differential Equations (HODEs)

```@eval
using GeometricEquations, Markdown
Markdown.parse(GeometricEquations.hode_equations)
```


### Implicit Ordinary Differential Equations (IODEs)

```@eval
using GeometricEquations, Markdown
Markdown.parse(GeometricEquations.iode_equations)
```


### Lagrangian Ordinary Differential Equations (LODEs)

```@eval
using GeometricEquations, Markdown
Markdown.parse(GeometricEquations.lode_equations)
```


## Differential Algebraic Equation (DAE)

```@eval
using GeometricEquations, Markdown
Markdown.parse(GeometricEquations.dae_equations)
```

```@eval
using GeometricEquations, Markdown
Markdown.parse(GeometricEquations.dae_example)
```


### Partitioned Differential Algebraic Equation (PDAE)

```@eval
using GeometricEquations, Markdown
Markdown.parse(GeometricEquations.pdae_equations)
```


### Hamiltonian Differential Algebraic Equation (HDAE)

```@eval
using GeometricEquations, Markdown
Markdown.parse(GeometricEquations.hdae_equations)
```


### Implicit Differential Algebraic Equation (IDAE)

```@eval
using GeometricEquations, Markdown
Markdown.parse(GeometricEquations.idae_equations)
```


### Lagrangian Differential Algebraic Equation (LDAE)

```@eval
using GeometricEquations, Markdown
Markdown.parse(GeometricEquations.ldae_equations)
```


## Stochastic Differential Equations (SDEs)

```@eval
using GeometricEquations, Markdown
Markdown.parse(GeometricEquations.sde_equations)
```

```@eval
using GeometricEquations, Markdown
Markdown.parse(GeometricEquations.sde_examples)
```

### Partitioned Stochastic Differential Equation (PSDE)

```@eval
using GeometricEquations, Markdown
Markdown.parse(GeometricEquations.psde_equations)
```

```@eval
using GeometricEquations, Markdown
Markdown.parse(GeometricEquations.psde_examples)
```

### Split Partitioned Stochastic Differential Equation (SPSDE)

```@eval
using GeometricEquations, Markdown
Markdown.parse(GeometricEquations.spsde_equations)
```

```@eval
using GeometricEquations, Markdown
Markdown.parse(GeometricEquations.spsde_examples)
```
