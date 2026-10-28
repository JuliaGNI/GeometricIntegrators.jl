# Integrators


This section describes how to implement, extend, and test numerical integrators in GeometricIntegrators.jl. It is intended for developers who want to contribute a new integration method or modify an existing implementation.

For the mathematical formulation and user-facing interface of existing methods, see the [Integrators page: RK](../modules/integrators.md), is where the link leads, but other examples exist. This guide focuses on the software architecture and development workflow.

---
## Architecture
### 1.1 Architectural Overview

`GeometricIntegrators.jl` enforces a strict separation between problem statements, numerical tableaux/methods, and the concrete execution algorithms that advance the system state forward in time.

```
┌───────────────────────────┐     ┌───────────────────────────┐
│     GeometricProblem      │     │      GeometricMethod      │
│  (Equations + Init State) │     │  (Tableau / Formulation)  │
└─────────────┬─────────────┘     └─────────────┬─────────────┘
              │                                 │
              └────────────────┬────────────────┘
                               ▼
                ┌─────────────────────────────┐
                │     AbstractIntegrator      │
                └──────────────┬──────────────┘
                               │
       ┌───────────────────────┼───────────────────────┐
       ▼                       ▼                       ▼
ExplicitIntegrator      ImplicitIntegrator      PartitionedIntegrator
 (e.g., Explicit RK)     (e.g., FIRK, DIRK)     (e.g., Symplectic Euler)
```

The concrete `Integrator` instance binds a system of differential equations with a computational method and preallocated memory buffers, producing a callable object that executes time steps without heap allocation.

---
### 1.2 The ingredients
 
An integration is assembled from these pieces:
 
| Ingredient | Role |
|---|---|
| **Problem** | The continuous system: equations, initial conditions, time span, time step, parameters (an `AbstractProblem`, with ODE problems as `AbstractProblemODE`). |
| **Method** | A description of the discrete algorithm (tableau, basis, quadrature, flags). Carries no run-time state. |
| **Cache** | Preallocated work arrays used inside a step. |
| **Solver** | Nonlinear solver for implicit methods (default `Newton()` for implicit Runge–Kutta), from SimpleSolvers. |
| **Initial guess** | Predictor that seeds the nonlinear solver (default `HermiteExtrapolation()` for implicit Runge–Kutta). |
| **Projection** | Optional correction applied after a step, see the *Projections* developer page. |
| **Solution** | Storage for the trajectory, filled by `integrate!`. |
 
These are bound together into a `GeometricIntegrator`:
 
````julia
int = GeometricIntegrator(problem, method)   # allocates cache, solver, initial guess
sol = integrate(problem, method)             # does the above and runs the whole integration
````
 
The integrator exposes its parts through accessors used throughout the source: `problem(int)`, `equations(int)`, `cache(int)`, `cache(int, ST)`, `solver(int)`, `solverstate(int)`, `iguess(int)`, `timestep(int)`.
 
**Design rule:** a method describes, a cache stores, an integrator orchestrates. A `Method` is never mutated during a run.
---

### 1.3 Abstraction & Type Hierarchy

All integration algorithms in `GeometricIntegrators.jl` subtype `AbstractIntegrator`. Subtyping dictates how stage systems are solved, how state vectors are partitioned, and how memory caches are organized.

```
AbstractIntegrator
├── ExplicitIntegrator
│   ├── ExplicitRKIntegrator
│   └── SplittingIntegrator
├── ImplicitIntegrator
│   ├── ImplicitRKIntegrator
│   └── VariationalIntegrator
├── PartitionedIntegrator
│   ├── PRKIntegrator
│   └── SymplecticEulerIntegrator
└── DifferentialAlgebraicIntegrator
    ├── DAEIntegrator
    └── SPARKIntegrator
```

### Abstraction Categories

1. **`AbstractIntegrator`**: Root type for all time-stepping implementations.
2. **`ExplicitIntegrator`**: Algorithms where the solution at $t_{n+1}$ is computed via direct function evaluations of the vector fields without solving linear or nonlinear algebraic systems.
3. **`ImplicitIntegrator`**: Algorithms requiring nonlinear solver loops (e.g., Newton-Raphson or Quasi-Newton) or linear system solves at each time step.
4. **`PartitionedIntegrator`**: Methods tailored for split phase spaces (e.g., canonical Hamiltonian systems split into position $q$ and momentum $p$), such as Partitioned Runge-Kutta (PRK) schemes.
5. **`DifferentialAlgebraicIntegrator`**: Specialized frameworks (including Symmetric Additive Runge-Kutta / SPARK schemes) designed to enforce algebraic constraints in DAEs or projected manifold systems.

---



### 1.4 Caches, `ST`, and `DT`
 
Caches are not plain structs passed around by the user. The integrator holds a `CacheDict` that stores one cache per element type, and code retrieves the right one:
 
````julia
local C = cache(int, ST)     # in components!/residual!, ST = element type of the solver vector
````
 
The reason is automatic differentiation. When the nonlinear solver computes a Jacobian, the unknowns carry dual numbers, so stage arrays must have the *solver's* element type `ST`, not only the problem's `DT`. Rules that follow from this:
 
- Methods look up their cache type through `CacheType(ST, problem, method)`, which must be inferable **from the method alone**, because it is evaluated on every `cache(int, ST)` call. A `CacheType` that reads a value off the problem does not constant-fold and makes the Newton hot path box on every stage access. Pass dimensions to the cache *constructor*, not as type parameters.
- Quantities that depend on solver unknowns live in the `ST` cache. Quantities that do not (for example work done after the solve) can use the `DT` cache.
- Some methods keep values such as the projection parameter $\lambda$ in the cache because $\lambda$ is itself a nonlinear unknown; those arrays must be `ST`-typed.
### 1.5 Algorithm families and where they live
 
| Family | Documentation page | Typical types | Notes |
|---|---|---|---|
| Runge–Kutta (explicit, diagonally implicit, fully implicit) | Runge-Kutta | `ERK`, `DIRK`, `IRK`, [`RK`](../integrators/rk.md) | Driven by a tableau from RungeKutta.jl. |
| Partitioned Runge–Kutta | Runge-Kutta | `EPRK`, `IPRK`, `PRK` | Pair of tableaus for `q` and `p`. |
| Splitting and composition | Splitting | splitting and composition integrators | Reuses sub-integrators; exact flows are supported. |
| Variational partitioned RK | VPRK | [`VPRK`](../integrators/vprk.md), `VPRKGauss`, `VPRKLobattoIIIA`, ... | Symplectic conditions on the tableau pair. |
| Variational (discrete Euler–Lagrange) | Variational | `DiscreteEulerLagrange`, `PMVImidpoint`, `PMVItrapezoidal` | Position–momentum form from a discrete Lagrangian. |
| Degenerate variational | DVI | [`DVIA`](../integrators/dvi.md), `DVIB`, `CMDVI`, `CTDVI`, `DVRK` | Lagrangians linear in the velocities. |
| Galerkin variational | CGVI, DGVI | [`CGVI`](../integrators/cgvi.md), `CGVINodal`, DGVI variants | Need a basis and a quadrature rule. |
| Hamilton–Pontryagin | HPI, HPG | [`Hamilton–Pontryagin integrators`](../integrators/hpi.md) | |
| Additive / SPARK | SPARK | SPARK and SLRK methods | Index-2 DAEs and constrained systems. |
| Projected / energy-preserving RK | Runge-Kutta | `PGLRK`, `FLRK`, `VPRKpTableau` | Some are experimental; read their docstrings for status. |

Source files are grouped by family under `src/integrators/` (for example `rk/`, `vi/`, `dvi/`, `cgvi/`, `dgvi/`), one file per integrator family variant, with shared helpers such as `rk/common.jl` and `solver_defaults.jl`. 
 
### 1.6 The surrounding ecosystem
 
| Piece | Source | What it provides |
|---|---|---|
| Equation and problem types | GeometricEquations, GeometricProblems | `ODE`, `PODE`, `HODE`, `IODE`, `LODE`, `DAE`, `SDE`, ... and test problems. |
| Base types | GeometricBase | Shared abstractions and state variables. |
| Runge–Kutta tableaus | RungeKutta.jl | Tableaus and algebraic predicates (explicit, symplectic, ...). |
| Nonlinear solvers | SimpleSolvers.jl | `Newton`, `QuasiNewton`, `DogLeg`, `Picard`, linear solvers, line search. |
| Methods and Solutions | this package | Merged into GeometricIntegrators in v0.13. |
| Projections | this package | Module of GeometricIntegrators (revised in v0.13). |
 
A new integrator should use only the documented public interface of these modules. If an integrator needs a new capability from one of them (a new tableau predicate, a new equation evaluation), add it to that module rather than working around it locally.
 
---

## 2. Choosing an implementation strategy
 
Decide what you are adding before writing code. The answer fixes the cache layout, whether a solver is needed, and which tests apply.
 
### 2.1 Decision questions
 
1. **Problem class.** ODE, partitioned, Hamiltonian, implicit/Lagrangian, DAE, stochastic. Target the most specific class the method actually exploits.
2. **Structure to preserve.** Symplecticity, symmetry, energy, momentum, constraints. If the answer is only "order", reuse an existing Runge–Kutta integrator.
3. **Explicit or implicit.** Explicit methods are cheap and need no solver, but cannot be symplectic for general non-separable Hamiltonians. Implicit methods cost a nonlinear solve per step and give the structural properties.
4. **History.** Variational and multistep methods need data from previous steps, which affects the first step.
5. **Is it a new tableau or a new algorithm?**
### 2.2 Decision table
 
| Situation | Strategy |
|---|---|
| New Runge–Kutta or partitioned RK coefficients | Add a **tableau** and expose a constructor. Reuse `ERK`/`IRK`/`PRK`. |
| New splitting or composition scheme | Add **coefficients** and reuse the splitting integrator. |
| Separable Hamiltonian, cheap and symplectic | **Splitting/composition** or an explicit partitioned RK. |
| Non-separable Hamiltonian, symplectic | **Implicit** symplectic method (Gauss, VPRK, Galerkin). |
| Lagrangian system, momentum preservation | **Variational integrator** from a discrete Lagrangian. |
| Degenerate Lagrangian, Dirac constraints | **DVI/DVRK** or **projected VPRK**. |
| Invariants must hold exactly | **Projection** attached to a base method, or a method with exact invariants. |
| Holonomic or Dirac constraints, index-2 DAEs | **SPARK** or constrained variational integrator. |
| New algorithm with its own step equations | **New integrator type** following §3. |
 
### 2.3 Reuse before writing
 
Prefer, in order:
 
1. a new tableau for an existing integrator,
2. a projection or wrapper around an existing integrator,
3. a new integrator that reuses the existing cache helpers, solver, and initial-guess infrastructure,
4. a fully new integrator, only when the stage equations genuinely differ.
### 2.4 Nonlinear solver and initial guess
 
| Choice | Options | Guidance |
|---|---|---|
| Nonlinear solver | `Newton`, `QuasiNewton`, `DogLeg`, `Picard` | `Newton()` is the default for implicit Runge–Kutta. `Picard` needs no linear solver and is suited to small steps where the iteration contracts. |
| Linear solver | pure-Julia `LU()` by default | The package sets the linear solver of `Newton`, `QuasiNewton`, and `DogLeg` explicitly. A LAPACK LU orders operations differently, and for a numerically singular stage system can return a silently wrong result instead of raising `SingularException`. A caller who names `linear_solver_method` overrides this default. |
| Jacobian | automatic differentiation through dual numbers | This is why caches must support the solver element type `ST`. Hand-coded Jacobians are an optimization when profiling justifies them. |
| Initial guess | `HermiteExtrapolation()` by default for implicit RK | Higher-order extrapolation reduces iterations but must not change the converged answer. |
 
Set defaults per method with `default_solver(method)` and `default_iguess(method)`.
 
### 2.5 Cost note for the docstring
 
State, for the new method: number of vector-field evaluations per step, size of the nonlinear system (`solversize`), and whether a Jacobian is needed. Reviewers use this to compare against existing methods.
 
---
 
## 3. Implementing new integrators
 
### 3.1 The interface
 
These are the functions that appear in the Runge–Kutta integrators. A new integrator implements the ones that apply to it.
 
| Function | Needed for | Purpose |
|---|---|---|
| concrete `Method` type, abstract method type with traits | all | Algorithm description. |
| `isexplicit(method)`, `isimplicit(method)` | all | Traits, defined for both the instance and `Type{<:...}`. |
| `initmethod(method, problem)` | all | Adapt the method (for example its tableau element type) to the problem's number types. |
| `Cache{ST}(problem, method; kwargs...)` | all | Construct the cache. |
| `CacheType(ST, problem, method)` | all | Return the cache type, as a function of the method alone (§1.4). |
| `integrate_step!(sol, history, params, int)` | all | Advance one step in place. |
| `update!(sol, params, ...)` | all | Final state update (explicit methods call it from `integrate_step!`). |
| `components!(x, sol, params, int)` | all | Compute stage values and vector fields. For explicit methods `x` is `nothing`. |
| `internal_variables(method, problem)` | optional | Named tuple of arrays stored in each solution step for output and restart. |
| `copy_internal_variables!(solstep, cache)` | optional | Copy cache stages into the solution step. |
| `default_solver(method)`, `default_iguess(method)` | implicit | Defaults for the nonlinear solver and predictor. |
| `solversize(method, problem)` | implicit | Length of the nonlinear unknown vector. |
| `nlsolution(cache)` | implicit | The cache's solver vector. |
| `initial_guess!(sol, history, params, int)` | implicit | Fill the solver vector before the solve. |
| `residual!(b, sol, params, int)` and `residual!(b, x, sol, params, int)` | implicit | Nonlinear system; the second form calls `components!` first. |
| `check_solver_status(status, int)` | implicit | Act on the solver outcome. |
| `Base.show(io, int)` | recommended | Human-readable summary. |
 
Arguments of the step functions:
 
- `sol`: the current solution step. It carries `sol.t` and `sol.q` (and `sol.p` for partitioned problems), and is updated in place.
- `history`: previous solution steps, for example `history[1]`.
- `params`: the equation parameters, forwarded to the equation functions, as in `equations(int).v(V, t, q, params)`. These are **not** integrator parameters.
- `int`: the `GeometricIntegrator`.
### 3.2 Rules every implementation follows
 
- **In place.** `integrate_step!` writes the result into `sol`.
- **No allocations in the step.** All temporaries live in the cache. Use the helpers `create_internal_stage_vector(DT, D, S)` and `create_internal_stage_matrix(DT, D, M, S)` to allocate stage storage. `[VERIFY: argument list of create_internal_stage_matrix]`
- **Type stable.** Never hardcode `Float64`; use `DT`, `TT`, `ST`.
- **No global state.** Anything mutable belongs to the cache or solver state.
- **Use the accessors.** `timestep(int)`, `tableau(int)`, `eachstage(int)`, `equations(int)`, never fields directly.
- **Stage times.** Time at stage `i` is computed relative to `sol.t`, as `sol.t + timestep(int) * (c[i] - 1)`, because `sol.t` is already the end of the step when the step function runs. `[VERIFY: confirm convention for non-RK integrators]`
### 3.3 Annotated example: the explicit Runge–Kutta integrator
 
This is the actual structure of `integrators_erk.jl`, reduced to the parts a new explicit method must provide.
 
````julia
# 1. Abstract method type carrying traits
abstract type ERKMethod <: RKMethod end
isexplicit(method::Union{ERKMethod, Type{<:ERKMethod}}) = true
isimplicit(method::Union{ERKMethod, Type{<:ERKMethod}}) = false
 
# 2. Concrete method: immutable, holds the tableau, checks it at construction
struct ERK{TT <: Tableau} <: ERKMethod
    tableau::TT
    function ERK(tableau::TT) where {TT <: Tableau}
        @assert RungeKutta.isexplicit(tableau)
        new{TT}(tableau)
    end
end
 
# 3. Adapt the method to the problem's number types
initmethod(method::ERKMethod, ::GeometricProblem{ST, DT, TT}) where {ST, DT, TT} =
    ERK(method, TT)
 
# 4. Cache: every temporary the step uses
struct ERKCache{DT, S} <: ODEIntegratorCache{DT}
    Q::Vector{Vector{DT}}      # internal stages of the solution
    V::Vector{Vector{DT}}      # internal stages of the vector field
    function ERKCache{DT, S}(ics) where {DT, S}
        D = length(vec(ics.q))
        new(create_internal_stage_vector(DT, D, S),
            create_internal_stage_vector(DT, D, S))
    end
end
 
function Cache{ST}(problem::AbstractProblemODE, method::ERK; kwargs...) where {ST}
    S = nstages(tableau(method))
    ERKCache{ST, S}(initial_conditions(problem); kwargs...)
end
 
@inline CacheType(ST, ::AbstractProblemODE, method::ERK) =
    ERKCache{ST, nstages(tableau(method))}
 
# 5. Stages: compute Q and V without allocating
function components!(_, sol, params, int::GeometricIntegrator{<:ERK, <:AbstractProblemODE})
    local Q = cache(int).Q
    local V = cache(int).V
    for i in eachstage(int)
        tᵢ = sol.t + timestep(int) * (tableau(int).c[i] - 1)
        for k in eachindex(Q[i], V[i])
            yᵢ = 0
            for j in 1:(i - 1)
                yᵢ += tableau(int).a[i, j] * V[j][k]
            end
            Q[i][k] = sol.q[k] + timestep(int) * yᵢ
        end
        equations(int).v(V[i], tᵢ, Q[i], params)
    end
end
 
# 6. Update and step
function update!(sol, params, _, int::GeometricIntegrator{<:ERK})
    components!(nothing, sol, params, int)
    update!(sol.q, cache(int).V, tableau(int), timestep(int))
end
 
function integrate_step!(sol, history, params,
                         int::GeometricIntegrator{<:ERK, <:AbstractProblemODE})
    update!(sol, params, nothing, int)
end
````
 
What to notice:
 
- The method stores only the tableau. Everything that changes lives in the cache.
- The inner loops index preallocated arrays and call the equation function `v` in place. Nothing in `components!` allocates.
- `integrate_step!` is three lines because the work is in `components!` and `update!`, which are reused by other code paths.
- The same file defines `internal_variables` and `copy_internal_variables!` so stage values can be stored in the solution if requested.
### 3.4 Annotated example: what an implicit integrator adds
 
`integrators_irk.jl` adds a solver, an initial guess, a residual, and a solver status check.
 
````julia
abstract type IRKMethod <: RKMethod end
isexplicit(method::Union{IRKMethod, Type{<:IRKMethod}}) = false
isimplicit(method::Union{IRKMethod, Type{<:IRKMethod}}) = true
 
default_solver(::IRKMethod) = Newton()
default_iguess(::IRKMethod) = HermiteExtrapolation()
 
# Unknowns: the stacked stage increments Y_1, ..., Y_s, each of length D
solversize(method::IRKMethod, problem::AbstractProblemODE) =
    length(vec(initial_conditions(problem).q)) * nstages(method)
 
nlsolution(cache::IRKCache) = cache.x      # the cache owns the solver vector
 
# Seed the solver from the initial-guess strategy
function initial_guess!(sol, history, params,
                        int::GeometricIntegrator{<:IRK, <:AbstractProblemODE})
    for i in eachstage(int)
        soltmp = (t = history[1].t + timestep(int) * tableau(int).c[i],
                  q = cache(int).Q[i], q̇ = cache(int).V[i])
        solutionstep!(soltmp, history, problem(int), iguess(int))
    end
    # assemble x[i-th block] = h * Σ_j a[i,j] * V[j]
    ...
end
 
# Stage values from the solver vector x, using the ST cache for AD
function components!(x::AbstractVector{ST}, sol, params,
                     int::GeometricIntegrator{<:IRK, <:AbstractProblemODE}) where {ST}
    local C = cache(int, ST)
    local D = length(C.V[1])
    for i in eachindex(C.Q, C.Y)
        for k in eachindex(C.Q[i], C.Y[i])
            C.Y[i][k] = x[D * (i - 1) + k]
            C.Q[i][k] = sol.q[k] + C.Y[i][k]
        end
    end
    for i in eachindex(C.Q, C.V)
        tᵢ = sol.t + timestep(int) * (tableau(int).c[i] - 1)
        equations(int).v(C.V[i], tᵢ, C.Q[i], params)
    end
end
 
# Residual: components! first, then the algebraic stage equations
function residual!(b::AbstractVector{ST}, x::AbstractVector{ST}, sol, params,
                   int::GeometricIntegrator{<:IRK, <:AbstractProblemODE}) where {ST}
    components!(x, sol, params, int)
    residual!(b, sol, params, int)
end
 
# Step: solve, check, update
function integrate_step!(sol, history, params,
                         int::GeometricIntegrator{<:IRK, <:AbstractProblemODE})
    solverstatus = solve_with_status!(nlsolution(int), solver(int), solverstate(int),
                                      (sol, params, int))
    check_solver_status(solverstatus, int)
    update!(sol, params, nlsolution(int), int)
end
````
 
The mathematics behind this code is the stage system
 
````
Q_i = q_n + Y_i,        Y_i = h Σ_j a_ij V_j,        V_j = v(t_j, Q_j),        i = 1, …, s
````
 
with unknown vector `x = (Y_1, …, Y_s)` of length `D·s`, residual `b_i = Y_i − h Σ_j a_ij V_j`, and update `q_{n+1} = q_n + h Σ_i b_i V_i`. Note that the solver unknowns are the stage *increments* `Y`, not the stage velocities.
 
For implicit and Lagrangian problems (`IODE`, `LODE`) the stage equations are replaced by
 
````
P_i = ϑ(Q_i, V_i),     Q_i = q_n + h Σ_j a_ij V_j,
F_i = f(Q_i, V_i),     P_i = p_n + h Σ_j ā_ij F_j,
````
 
and the tableau pair `(a, ā)`, `(b, b̄)` must satisfy `b_i ā_ij + b̄_j a_ji = b_i b̄_j` and `b̄_i = b_i` for the method to be symplectic.
 
### 3.5 Step-by-step procedure for a new integrator
 
1. **Choose the strategy** (§2). If a new tableau suffices, stop here and add it.
2. **Define the method types.** Abstract type with traits, concrete type with parameters, constructor that validates its input.
3. **Define `initmethod`** so the method adopts the problem's number types.
4. **Define the cache** type, its constructor, `Cache{ST}(problem, method)`, and `CacheType`. Keep `CacheType` a function of the method alone.
5. **Implement `components!`.** Stage values and vector fields, no allocation.
6. **For implicit methods:** `solversize`, `nlsolution`, `default_solver`, `default_iguess`, `initial_guess!`, `residual!`, and the solver call with status check in `integrate_step!`.
7. **Implement `update!` and `integrate_step!`.**
8. **Optionally** `internal_variables`, `copy_internal_variables!`, and `Base.show`.
9. **Export** the method and add the file to the module's include list.
10. **Write tests and documentation** (§5, §6).
### 3.6 Notes for specific families
 
- **Runge–Kutta.** Tableau storage and validation come from RungeKutta.jl. Constructors assert the tableau class (an explicit integrator asserts an explicit tableau).
- **Partitioned and variational.** The method holds a pair of tableaus. State the symplecticity conditions in the docstring and verify them in tests. Initial momentum must be consistent where the method requires it, for example `p₀ = ϑ(q₀)` for degenerate variational Runge–Kutta methods.
- **Galerkin variational.** The method stores a basis and a quadrature rule and precomputes the mass matrix, derivative matrix, and boundary reconstruction vectors. The dimension and solver size are constructor *arguments*, not type parameters, to keep `CacheType` constant-foldable. Constructors reject bases that violate requirements (a nodal CGVI rejects interior-only nodes).
- **Splitting and composition.** Reuse existing sub-integrators; the new code manages ordering and step-size coefficients. An exact flow can be wrapped so it satisfies the ordinary integrator interface.
- **SPARK.** Works with a collection of tableaus rather than one, so the nonlinear system is larger and the cache layout reflects the additive structure.
- **Projected methods.** Attach a projection or perturb the tableau. Expose any fallback counters in the cache so users can see when the method degraded to its base method.
### 3.7 Pitfalls
 
- Allocating inside the step through closures, slices (use `@views`), or broadcasting over temporaries.
- Hardcoding `Float64`, which breaks high-precision and AD workflows.
- Storing mutable state in the method struct.
- Reading `C = cache(int)` where `cache(int, ST)` is required, so dual numbers meet `Float64` arrays and the Jacobian fails.
- Making `CacheType` depend on the problem, which defeats constant folding and causes boxing in the Newton loop.
- Warm-starting a scalar sub-solve from the previous trial, which can make the residual depend on the sequence of trials rather than on the unknown alone.
- Calling `check_solver_status` on probes that are not time steps, which rejects steps that were never meant to converge.
- Mixing up stage values and the updated state when reconstructing the solution.
---
 
## 4. Integrator lifecycle and solution interface
 
### 4.1 Entry points
 
````julia
sol = integrate(problem, method)                  # build everything and run
int = GeometricIntegrator(problem, method; solver = ..., iguess = ...)   # [VERIFY: keyword names]
integrate!(sol, int)                              # fill an existing solution
````
 
For many standard methods the integrator is selected implicitly from an equation and a tableau. The high-level functions are thin wrappers around the per-step function `integrate_step!`.
 
### 4.2 Construction
 
When a `GeometricIntegrator` is built, the package:
 
1. adapts the method to the problem with `initmethod`,
2. builds the cache dictionary through `Cache{ST}` and `CacheType`,
3. creates a solver (`default_solver`) sized by `solversize` and its state,
4. creates the initial-guess object (`default_iguess`),
5. attaches a projection when the problem or method requires one.
All of these can be overridden by the caller, which keeps maximum control with the user while providing sensible defaults.
 
### 4.3 The step
 
Each call to `integrate_step!(sol, history, params, int)`:
 
1. reads the current state from `sol` (and previous steps from `history`),
2. for implicit methods, computes an initial guess with `initial_guess!` and solves the nonlinear system with `solve_with_status!`, then checks the outcome with `check_solver_status`,
3. computes stage values and vector fields in `components!`,
4. writes the new state into `sol` in `update!`.
````
integrate!  ──►  loop over steps  ──►  integrate_step!(sol, history, params, int)
                                          │
              ┌───────────────────────────┴───────────────────────────┐
              │ implicit only                                         │ all methods
              ▼                                                       ▼
   initial_guess! → solve_with_status! → check_solver_status    components! → update!
````
 
The driver then handles storing the solution at the requested interval and applying a projection if one is attached. 
 
### 4.4 The solution interface
 
- **`sol`** is the mutable solution step. It holds `t` and `q` (plus `p` and multipliers where the problem has them), and the step function updates it in place. At the time the step function runs, `sol.t` is the end of the step and `sol.q` is the state to be advanced, which is why stage times are written `sol.t + h (c_i − 1)`.
- **`history`** gives access to previous solution steps (`history[1]` is the previous one) for initial guesses and for methods that need more than the last state.
- **`params`** are the equation parameters, passed through to the equation functions.
- **Internal variables.** If a method defines `internal_variables(method, problem)`, the returned arrays (typically stage values `Q`, `V`, `Y`) are stored in each solution step and can be used for post-processing, checkpointing, or restart. `copy_internal_variables!(solstep, cache)` copies the cache contents into the step, and only fills entries that exist.
- **Ensembles** of initial conditions are integrated through the ensemble support. 
### 4.5 Failure handling
 
A step can fail through solver non-convergence, a singular stage system, or a violated consistency condition. The convention is:
 
- the solver reports a status, and `check_solver_status(status, int)` decides whether to throw or continue,
- consistency requirements are checked at construction or initialization and warn or error clearly (for example a constructor can warn when tableau conditions fail, with a keyword such as `check_conditions = false` to silence it),
- a method that falls back to a simpler scheme records it in a counter rather than hiding it.
### 4.6 Adaptive time stepping
 
There is no general adaptivity infrastructure. An integrator that wants its own error control must still deliver output on the predefined equidistant time grid. See the *Adaptive Time Stepping* developer page.
 
---
 
## 5. Testing and mathematical verification
 
A new integrator is accepted when it is shown to be **correct** (right order), **structure preserving** (the properties it claims), and **well behaved** (no allocations, type stable, solver converges). Tests live under `test/`, with method-level checks (for example coefficient identities) and longer convergence studies kept separate. `
 
### 5.1 Construction and types
 
- Builds for each supported problem type and for several number types where meaningful.
- `@inferred` passes on construction and on a step.
- Traits are correct: `isexplicit`, `isimplicit`, and the structure predicates.
### 5.2 Order of convergence
 
**Static check.** For Runge–Kutta methods, check the order conditions on the tableau, for example
 
- order 1: `Σ b_i = 1`
- order 2: additionally `Σ b_i c_i = 1/2`
- order 3: additionally `Σ b_i c_i² = 1/3` and `Σ b_i a_ij c_j = 1/6`
and compare with the claimed order.
 
**Dynamic check.** Integrate a problem with a known solution for a fixed interval with step sizes `h, h/2, h/4, …`, measure the global error `e_k` at the final time, and compute the observed order
 
````
p ≈ log(e_k / e_{k+1}) / log(h_k / h_{k+1})
````
 
````julia
@testset "convergence order" begin
    errors = [norm(final_state(integrate(with_timestep(prob, h), method)) - exact(T), Inf)
              for h in hs]
    orders = log2.(errors[1:end-1] ./ errors[2:end])
    @test isapprox(mean(orders), target_order; atol = 0.2)
end
````
 
Use step sizes in the asymptotic regime and `BigFloat` for high orders so round-off does not hide the slope. Without a closed form, use a much finer reference solution or the method of manufactured solutions. `with_timestep` and `final_state` are placeholders for the helpers in the test utilities. 
 
### 5.3 Structure preservation
 
Test only what the method claims, and name the tests accordingly.
 
**Symplecticity.**
- *Tableau condition.* Runge–Kutta: `b_i a_ij + b_j a_ji = b_i b_j`. Partitioned: `b_i ā_ij + b̄_j a_ji = b_i b̄_j` and `b̄_i = b_i`.
- *Step-map condition.* With the Jacobian `Dφ` of the one-step map (automatic differentiation or finite differences) and the canonical matrix `J`, check `Dφᵀ J Dφ = J` to solver tolerance. For finite differences, show the defect shrinks as the difference step is refined.
- If symplecticity is *not* established, let the trait return `missing` rather than `true`, and say so in the docstring.
**Energy behavior.** Over a long run (for example `T = 10⁵` on a harmonic oscillator, pendulum, or Kepler problem), a symplectic method shows a bounded oscillation of `H` of size `O(hᵖ)` and no secular drift. Test that `max |H(t) − H(0)|` does not grow with the number of steps.
 
**Reversibility.** For a symmetric method, `Φ_h⁻¹ = Φ_{−h}`. Integrate forward `n` steps with `h` and back `n` steps with `−h`, and compare with the initial data to solver tolerance.
 
**Invariants and constraints.** Quadratic invariants (angular momentum) are preserved to round-off by Gauss methods. Constraints must hold to solver tolerance at every step.
 
**Exact conservation versus symplecticity.** For a general Hamiltonian system a method cannot be both exactly energy conserving and symplectic unless it reproduces the exact flow (Ge and Marsden, 1988). Energy-projected variants therefore set the symplecticity trait to `missing`.
 
### 5.4 Solver behavior (implicit methods)
 
- Converges in a bounded number of iterations on a benchmark at a reasonable step size.
- Fails with a clear status, not a silently wrong answer, when the step is too large or the stage system is singular.
- The converged solution does not depend on the initial-guess strategy beyond the tolerance.
### 5.5 Allocations and performance
 
````julia
integrate_step!(sol, history, params, int)               # warm up
@test (@allocated integrate_step!(sol, history, params, int)) == 0
````
 
Also compare cost and accuracy against an existing method on one reference problem (work–precision data, not just timings).
 
### 5.6 Reference problems
 
| Problem | Used for |
|---|---|
| Harmonic oscillator | Order, symplecticity, energy (closed form exists) |
| Pendulum | Nonlinear Hamiltonian, energy drift, reversibility |
| Kepler | Long-time behavior, angular momentum and energy |
| Lotka–Volterra | Noncanonical structure, invariants |
| Stiff or Van der Pol (large μ) | Stability and solver robustness |
| Cartesian pendulum | Constraints and DAE handling |
| Degenerate Lagrangian test cases | DVI, DVRK, projected VPRK |
| Linear SDE with known moments | Strong and weak order for stochastic methods |
 
### 5.7 Interpreting failures
 
| Symptom | Likely cause |
|---|---|
| Observed order one lower than expected | Inconsistent initial data or history on the first step, or a low-order initial guess feeding a reconstruction. |
| Energy drifts in a "symplectic" method | Loose solver tolerance, or the symplecticity condition is not satisfied. |
| Noisy convergence slope in `Float64` | Round-off in the asymptotic regime; use `BigFloat`. |
| Allocations in the step | Closures, array slicing, non-concrete cache fields. |
| Jacobian errors only under AD | Cache arrays typed `DT` where `ST` is needed. |
| Bisection or inner solve never converges | Warm-starting an inner solve from the previous trial. |
 
---
 
## 6. Documentation and contribution checklist
 
### 6.1 What the documentation must contain
 
Use the existing docstrings (for example for `PGLRK` and `VPRKpTableau`) as models. A complete docstring has:
 
- the equations the method solves, written out,
- the **conditions** the method requires (tableau class, basis requirements, initial-data consistency),
- the structural properties it has, and an honest **status** where something is experimental or unproven,
- the **reference**,
- the **constructors**, with keyword arguments,
- the **fields** of the cache.
Template:
 
````julia
@doc raw"""
One-line description.
 
Longer description with the stage equations
 
```math
\begin{aligned}
...
\end{aligned}
```
 
**Requirements and caveats**
- ...
 
**Status**
Order: ... Symplectic: ... (or `missing` and why).
 
**Reference**
Author. Title. Journal, Volume, Pages, Year. doi: ...
 
**Constructors**
```julia
MyMethod(tableau; kwargs...)
```
"""
````
 
Also add:
 
- an entry in the family page under *Integrators* and in the module reference,
- a bibliography entry,
- a release-note entry,
- for methods new in the literature, a short derivation.
### 6.2 Contributor checklist
 
**Design**
- [ ] Problem class and structure to preserve are stated
- [ ] Reuse of an existing integrator, tableau, or projection was considered (§2.3)
**Implementation**
- [ ] Abstract and concrete method types, traits, and validating constructor
- [ ] `initmethod`, `Cache{ST}`, `CacheType` (a function of the method alone)
- [ ] `components!`, `update!`, `integrate_step!`, allocation-free
- [ ] Implicit: `solversize`, `nlsolution`, `default_solver`, `default_iguess`, `initial_guess!`, `residual!`, solver status check
- [ ] Caches are `ST`-typed where solver unknowns enter
- [ ] Works for number types other than `Float64`
- [ ] Optional: `internal_variables`, `copy_internal_variables!`, `Base.show`
- [ ] Exported and included in the module
**Tests**
- [ ] Construction, traits, and type inference
- [ ] Order verified (tableau conditions and a convergence study)
- [ ] Each claimed structural property tested (§5.3)
- [ ] Solver behavior tested for implicit methods
- [ ] Zero-allocation step test
- [ ] Benchmark against an existing method on one reference problem
**Documentation and process**
- [ ] Docstring complete (§6.1)
- [ ] Family page, module reference, bibliography, release notes updated
- [ ] CI passes on all supported Julia versions
- [ ] PR description states the math, the claimed properties, and the evidence
### 6.3 Questions a reviewer should be able to answer from the PR
 
1. What does this method guarantee, and where is that tested?
2. Could an existing method or tableau have done the job?
3. What is the cost per step relative to comparable methods?
4. What happens when the solver fails?


