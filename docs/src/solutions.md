## Solutions

```@setup 1
using GeometricIntegrators
using GeometricProblems.HarmonicOscillator
prob = harmonic_oscillator_ode()
```

In what we have seen so far, the solution was always automatically created by
the `integrate()` function. While this is often convenient, it is sometimes not
performant, e.g., when carrying out long-time simulations with intermediate
saving of the solution.
In such cases, it is better to preallocate a solution object by
```@example 1
sol = Solution(prob)
```
where the first argument is an equation, the second argument is the time step
and the third argument is the number of time steps that will be computed in one
integration step.
The call to the integrator is then made via
```@example 1
int = Integrator(prob, Gauss(1))
integrate!(sol, int)
```
If several integration cycles shall be performed, the `reset!()` function can be
used to copy the solution of the last time step to the initial conditions of the
solution,
```julia
for i in 1:10
    # integrate!(sol, int)
    #
    # save or process solution
    #
    # reset!(sol)
end
```
All solutions have a `t` field holding the series of time steps that has been
computed in addition to several data fields, for example `q` for an ODE solution,
`q` and `p` for a PODE solution, `q`and `λ` for a DAE solution, and `q`, `p` and
`λ` for a PDAE solution.
