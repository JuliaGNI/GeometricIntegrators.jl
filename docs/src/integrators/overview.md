```@meta
CurrentModule = GeometricIntegrators
```

# Overview of Available Integrators



## Integrators for ODEs




#### Explicit Runge-Kutta Methods

| Function                               | Order | Method                           |
|:---------------------------------------|:------|:---------------------------------|
| [`TableauExplicitEuler`](@ref)         | 1     | Explicit / Forward Euler         |
| [`TableauExplicitMidpoint`](@ref)      | 2     | Explicit Midpoint                |
| [`TableauHeun2`](@ref)                 | 2     | Heun's Method of order two       |
| [`TableauHeun3`](@ref)                 | 3     | Heun's Method of order three     |
| [`TableauRalston2`](@ref)              | 2     | Ralston's Method of order two    |
| [`TableauRalston3`](@ref)              | 3     | Ralston's Method of order three  |
| [`TableauRunge`](@ref)                 | 2     | Runge's Method                   |
| [`TableauKutta`](@ref)                 | 3     | Kutta's Method                   |
| [`TableauRK416`](@ref)                 | 4     | Explicit 4th order Runge-Kutta (1/6 rule) |
| [`TableauRK438`](@ref)                 | 4     | Explicit 4th order Runge-Kutta (3/8 rule) |


#### Diagonally Implicit Runge-Kutta Methods

| Function                               | Order | Method                           |
|:---------------------------------------|:------|:---------------------------------|
| [`TableauCrankNicolson`](@ref)         | 3     | Crank-Nicholson Method           |
| [`TableauKraaijevangerSpijker`](@ref)  | 3     | Kraaijevanger & Spijker's Method |
| [`TableauQinZhang`](@ref)              | 3     | Qin & Zhang's Method             |
| [`TableauCrouzeix`](@ref)              | 3     | Crouzeix's Method                |


#### Fully Implicit Runge-Kutta Methods

| Function                               | Order | Method                      |
|:---------------------------------------|:------|:----------------------------|
| [`TableauImplicitEuler`](@ref)         | 1     | Implicit / Backward Euler   |
| [`TableauImplicitMidpoint`](@ref)      | 2     | Implicit Midpoint           |
| [`TableauSRK3`](@ref)                  | 4     | Symmetric Runge-Kutta s=3   |


#### Gauß, Radau and Lobatto Methods

| Function                               | Order | Method                      |
|:---------------------------------------|:------|:----------------------------|
| [`TableauGauss`](@ref)                 | 2s    | Gauss-Legendre              |
| [`TableauRadauIA`](@ref)               | 2s-1  | Radau-IA                    |
| [`TableauRadauIB`](@ref)               | 2s-1  | Radau-IB                    |
| [`TableauRadauIIA`](@ref)              | 2s-1  | Radau-IIA                   |
| [`TableauRadauIIB`](@ref)              | 2s-1  | Radau-IIB                   |
| [`TableauLobattoIIIA`](@ref)           | 2s-2  | Lobatto-IIIA                |
| [`TableauLobattoIIIB`](@ref)           | 2s-2  | Lobatto-IIIB                |
| [`TableauLobattoIIIC`](@ref)           | 2s-2  | Lobatto-IIIC                |
| [`TableauLobattoIIIC̄`](@ref)           | 2s-2  | Lobatto-IIIC̄                |
| [`TableauLobattoIIID`](@ref)           | 2s-2  | Lobatto-IIID                |
| [`TableauLobattoIIIE`](@ref)           | 2s-2  | Lobatto-IIIE                |
| [`TableauLobattoIIIF`](@ref)           | 2s    | Lobatto-IIIF                |
| [`TableauLobattoIIIG`](@ref)           | 2s    | Lobatto-IIIG                |

All of these tableaus are generated on the fly and take the number of stages `s` as parameter.


#### Splitting Methods

| Function                               | Order | Method                       |
|:---------------------------------------|:------|:-----------------------------|
| [`TableauLieA`](@ref)                  | 1     | Lie-Trotter Splitting A      |
| [`TableauLieB`](@ref)                  | 1     | Lie-Trotter Splitting B      |
| [`TableauStrang`](@ref)                | 2     | Strang / Marchuk Splitting   |
| [`TableauMarchuk`](@ref)               | 2     | Strang / Marchuk Splitting   |
| [`TableauStrangA`](@ref)               | 2     | Strang / Marchuk Splitting A |
| [`TableauStrangB`](@ref)               | 2     | Strang / Marchuk Splitting B |
| [`TableauMcLachlan2`](@ref)            | 2     | McLachlan's 2nd order symmetric, minimum error composition method |
| [`TableauMcLachlan4`](@ref)            | 2     | McLachlan's 4th order symmetric, minimum error composition method |
| [`TableauTripleJump`](@ref)            | 4     | 4th order "Triple Jump" composition method                        |
| [`TableauSuzukiFractal`](@ref)         | 4     | Suzuki's 4th order "fractal" composition method                   |


## Integrators for partitioned ODEs

#### Partitioned Runge-Kutta Methods

| Function                               | Order | Method                      |
|:---------------------------------------|:------|:----------------------------|
| [`TableauLobattoIIIAIIIB`](@ref)       | 2s-2  | Lobatto-IIIA-IIIB           |
| [`TableauLobattoIIIBIIIA`](@ref)       | 2s-2  | Lobatto-IIIB-IIIA           |
| [`TableauLobattoIIICIIIC̄`](@ref)       | 2s-2  | Lobatto-IIIC-IIIC̄           |
| [`TableauLobattoIIIC̄IIIC`](@ref)       | 2s-2  | Lobatto-IIIC̄-IIIC           |


## Integrators for implicit ODEs




## Integrators for DAEs

