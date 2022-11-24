```@meta
CurrentModule = GeometricIntegrators
```

# Overview of Available Methods



## Integrators for ODEs

The main method types for ODEs currently implemented are Runge-Kutta methods and splitting methods.

### Runge-Kutta Methods

Any Runge-Kutta method can be selected by the [`RK`](@ref) method
```julia
rk = RK(<tableau>)
```
where `<tableau>` is any tableau from [RungeKutta.Tableaus](@ref). For most tableaus there also exist explicit shortcuts to select the method. These are listed in the following.


#### Explicit Runge-Kutta Methods

| Function                        | Order | Method                           |
|:--------------------------------|:------|:---------------------------------|
| [`ExplicitEuler`](@ref)         | 1     | Explicit / Forward Euler         |
| [`ExplicitMidpoint`](@ref)      | 2     | Explicit Midpoint                |
| [`Heun2`](@ref)                 | 2     | Heun's Method of order two       |
| [`Heun3`](@ref)                 | 3     | Heun's Method of order three     |
| [`Ralston2`](@ref)              | 2     | Ralston's Method of order two    |
| [`Ralston3`](@ref)              | 3     | Ralston's Method of order three  |
| [`Runge`](@ref)                 | 2     | Runge's Method                   |
| [`Kutta`](@ref)                 | 3     | Kutta's Method                   |
| [`RK416`](@ref)                 | 4     | Explicit 4th order Runge-Kutta (1/6 rule) |
| [`RK438`](@ref)                 | 4     | Explicit 4th order Runge-Kutta (3/8 rule) |


#### Diagonally Implicit Runge-Kutta Methods

| Function                        | Order | Method                           |
|:--------------------------------|:------|:---------------------------------|
| [`CrankNicolson`](@ref)         | 3     | Crank-Nicholson Method           |
| [`KraaijevangerSpijker`](@ref)  | 3     | Kraaijevanger & Spijker's Method |
| [`QinZhang`](@ref)              | 3     | Qin & Zhang's Method             |
| [`Crouzeix`](@ref)              | 3     | Crouzeix's Method                |


#### Fully Implicit Runge-Kutta Methods

| Function                        | Order | Method                      |
|:--------------------------------|:------|:----------------------------|
| [`ImplicitEuler`](@ref)         | 1     | Implicit / Backward Euler   |
| [`ImplicitMidpoint`](@ref)      | 2     | Implicit Midpoint           |
| [`SRK3`](@ref)                  | 4     | Symmetric Runge-Kutta s=3   |


#### Gauß, Radau and Lobatto Methods

| Function                        | Order | Method                      |
|:--------------------------------|:------|:----------------------------|
| [`Gauss`](@ref)                 | 2s    | Gauss-Legendre              |
| [`RadauIA`](@ref)               | 2s-1  | Radau-IA                    |
| [`RadauIB`](@ref)               | 2s-1  | Radau-IB                    |
| [`RadauIIA`](@ref)              | 2s-1  | Radau-IIA                   |
| [`RadauIIB`](@ref)              | 2s-1  | Radau-IIB                   |
| [`LobattoIIIA`](@ref)           | 2s-2  | Lobatto-IIIA                |
| [`LobattoIIIB`](@ref)           | 2s-2  | Lobatto-IIIB                |
| [`LobattoIIIC`](@ref)           | 2s-2  | Lobatto-IIIC                |
| [`LobattoIIIC̄`](@ref)           | 2s-2  | Lobatto-IIIC̄                |
| [`LobattoIIID`](@ref)           | 2s-2  | Lobatto-IIID                |
| [`LobattoIIIE`](@ref)           | 2s-2  | Lobatto-IIIE                |
| [`LobattoIIIF`](@ref)           | 2s    | Lobatto-IIIF                |
| [`LobattoIIIG`](@ref)           | 2s    | Lobatto-IIIG                |

All of these tableaus are generated on the fly and take the number of stages `s` as parameter.


#### Splitting Methods

| Function                               | Order | Method                       |
|:---------------------------------------|:------|:-----------------------------|
| [`LieA`](@ref)                  | 1     | Lie-Trotter Splitting A      |
| [`LieB`](@ref)                  | 1     | Lie-Trotter Splitting B      |
| [`Strang`](@ref)                | 2     | Strang / Marchuk Splitting   |
| [`Marchuk`](@ref)               | 2     | Strang / Marchuk Splitting   |
| [`StrangA`](@ref)               | 2     | Strang / Marchuk Splitting A |
| [`StrangB`](@ref)               | 2     | Strang / Marchuk Splitting B |
| [`McLachlan2`](@ref)            | 2     | McLachlan's 2nd order symmetric, minimum error composition method |
| [`McLachlan4`](@ref)            | 2     | McLachlan's 4th order symmetric, minimum error composition method |
| [`TripleJump`](@ref)            | 4     | 4th order "Triple Jump" composition method                        |
| [`SuzukiFractal`](@ref)         | 4     | Suzuki's 4th order "fractal" composition method                   |


## Integrators for partitioned ODEs

#### Partitioned Runge-Kutta Methods

| Function                        | Order | Method                      |
|:--------------------------------|:------|:----------------------------|
| [`LobattoIIIAIIIB`](@ref)       | 2s-2  | Lobatto-IIIA-IIIB           |
| [`LobattoIIIBIIIA`](@ref)       | 2s-2  | Lobatto-IIIB-IIIA           |
| [`LobattoIIICIIIC̄`](@ref)       | 2s-2  | Lobatto-IIIC-IIIC̄           |
| [`LobattoIIIC̄IIIC`](@ref)       | 2s-2  | Lobatto-IIIC̄-IIIC           |


## Integrators for implicit ODEs

All implicit Runge-Kutta and partitioned Runge-Kutta methods can also be applied to implicit ODEs.

| Function                        | Order | Method                      |
|:--------------------------------|:------|:----------------------------|
| [`ImplicitEuler`](@ref)         | 1     | Implicit / Backward Euler   |
| [`ImplicitMidpoint`](@ref)      | 2     | Implicit Midpoint           |
| [`SRK3`](@ref)                  | 4     | Symmetric Runge-Kutta s=3   |
|:--------------------------------|:------|:----------------------------|
| [`Gauss`](@ref)                 | 2s    | Gauss-Legendre              |
| [`RadauIA`](@ref)               | 2s-1  | Radau-IA                    |
| [`RadauIB`](@ref)               | 2s-1  | Radau-IB                    |
| [`RadauIIA`](@ref)              | 2s-1  | Radau-IIA                   |
| [`RadauIIB`](@ref)              | 2s-1  | Radau-IIB                   |
| [`LobattoIIIA`](@ref)           | 2s-2  | Lobatto-IIIA                |
| [`LobattoIIIB`](@ref)           | 2s-2  | Lobatto-IIIB                |
| [`LobattoIIIC`](@ref)           | 2s-2  | Lobatto-IIIC                |
| [`LobattoIIIC̄`](@ref)           | 2s-2  | Lobatto-IIIC̄                |
| [`LobattoIIID`](@ref)           | 2s-2  | Lobatto-IIID                |
| [`LobattoIIIE`](@ref)           | 2s-2  | Lobatto-IIIE                |
| [`LobattoIIIF`](@ref)           | 2s    | Lobatto-IIIF                |
| [`LobattoIIIG`](@ref)           | 2s    | Lobatto-IIIG                |
|:--------------------------------|:------|:----------------------------|
| [`LobattoIIIAIIIB`](@ref)       | 2s-2  | Lobatto-IIIA-IIIB           |
| [`LobattoIIIBIIIA`](@ref)       | 2s-2  | Lobatto-IIIB-IIIA           |
| [`LobattoIIICIIIC̄`](@ref)       | 2s-2  | Lobatto-IIIC-IIIC̄           |
| [`LobattoIIIC̄IIIC`](@ref)       | 2s-2  | Lobatto-IIIC̄-IIIC           |


## Integrators for Lagrangian ODEs

| Function                        | Method                                                                    |
|:--------------------------------|:--------------------------------------------------------------------------|
| [`VPRK`](@ref)                  | Variational Partitioned Runge-Kutta integrator                            |
| [`DegenerateVPRK`](@ref)        | Variational Partitioned Runge-Kutta integrator for degenerate Lagrangians |
| [`ProjectedVPRK`](@ref)         | Projected Variational Partitioned Runge-Kutta integrator                  |
|:--------------------------------|:--------------------------------------------------------------------------|
| [VPRKGauss](@ref)               | VPRK integrator with [TableauGauss](@ref)                                 |
| [VPRKRadauIIA](@ref)            | VPRK integrator with [TableauRadauIIA](@ref)                              |
| [VPRKRadauIIB](@ref)            | VPRK integrator with [TableauRadauIIB](@ref)                              |
| [VPRKLobattoIIIA](@ref)         | VPRK integrator with [TableauLobattoIIIA](@ref)                           |
| [VPRKLobattoIIIB](@ref)         | VPRK integrator with [TableauLobattoIIIB](@ref)                           |
| [VPRKLobattoIIIC](@ref)         | VPRK integrator with [TableauLobattoIIIC](@ref)                           |
| [VPRKLobattoIIIC̄](@ref)         | VPRK integrator with [TableauLobattoIIIC̄](@ref)                           |
| [VPRKLobattoIIID](@ref)         | VPRK integrator with [TableauLobattoIIID](@ref)                           |
| [VPRKLobattoIIIE](@ref)         | VPRK integrator with [TableauLobattoIIIE](@ref)                           |
| [VPRKLobattoIIIF](@ref)         | VPRK integrator with [TableauLobattoIIIF](@ref)                           |
| [VPRKLobattoIIIG](@ref)         | VPRK integrator with [TableauLobattoIIIG](@ref)                           |
| [VPRKLobattoIIIAIIIB](@ref)     | VPRK integrator with [TableauLobattoIIIAIIIB](@ref)                       |
| [VPRKLobattoIIIBIIIA](@ref)     | VPRK integrator with [TableauLobattoIIIBIIIA](@ref)                       |
| [VPRKLobattoIIIAIIIĀ](@ref)     | VPRK integrator with [TableauLobattoIIIAIIIĀ](@ref)                       |
| [VPRKLobattoIIIBIIIB̄](@ref)     | VPRK integrator with [TableauLobattoIIIBIIIB̄](@ref)                       |
| [VPRKLobattoIIICIIIC̄](@ref)     | VPRK integrator with [TableauLobattoIIICIIIC̄](@ref)                       |
| [VPRKLobattoIIIC̄IIIC](@ref)     | VPRK integrator with [TableauLobattoIIIC̄IIIC](@ref)                       |
| [VPRKLobattoIIIDIIID̄](@ref)     | VPRK integrator with [TableauLobattoIIIDIIID̄](@ref)                       |
| [VPRKLobattoIIIEIIIĒ](@ref)     | VPRK integrator with [TableauLobattoIIIEIIIĒ](@ref)                       |
| [VPRKLobattoIIIFIIIF̄](@ref)     | VPRK integrator with [TableauLobattoIIIFIIIF̄](@ref)                       |
| [VPRKLobattoIIIF̄IIIF](@ref)     | VPRK integrator with [TableauLobattoIIIF̄IIIF](@ref)                       |
| [VPRKLobattoIIIGIIIḠ](@ref)     | VPRK integrator with [TableauLobattoIIIGIIIḠ](@ref)                       |
|:--------------------------------|:--------------------------------------------------------------------------|
| [VPRKpInternal](@ref)           | VPRK integrator with projection on internal stages                        |
| [VPRKpLegendre](@ref)           | VPRK integrator with Legendre projection                                  |
| [VPRKpMidpoint](@ref)           | VPRK integrator with Midpoint projection                                  |
| [VPRKpSecondary](@ref)          | VPRK integrator with projection on secondary constraint                   |
| [VPRKpStandard](@ref)           | VPRK integrator with standard projection                                  |
| [VPRKpSymmetric](@ref)          | VPRK integrator with symmetric projection                                 |
| [VPRKpSymplectic](@ref)         | VPRK integrator with symplectic projection                                |
| [VPRKpVariational](@ref)        | VPRK integrator with variational projection                               |
| [VPRKpVariationalP](@ref)       | VPRK integrator with variational projection on P                          |
| [VPRKpVariationalQ](@ref)       | VPRK integrator with variational projection on Q                          |


## Integrators for DAEs

