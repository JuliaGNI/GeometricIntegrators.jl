
ode = ODE(1, x -> x, [1.])
@test typeof(Solution(ode, 0.1, 10)) <: SolutionODE

pode = PODE(1, (x, y) -> x, (x, y) -> 2y, [1.], [1.])
@test typeof(Solution(pode, 0.1, 10)) <: SolutionPODE

dae = DAE(2, 1, x -> x, (x, λ) -> [λ -λ], x -> x[2]-x[1], [1.; 1.], [0.])
@test typeof(Solution(dae, 0.1, 10)) <: SolutionDAE

pdae = PDAE(1, 1, (x, y) -> x, (x, y) -> 2y, (x, y, λ) -> λ, (x, y, λ) -> -λ, (x, y) -> x[1]-y[1], [1.], [1.], [0.])
@test typeof(Solution(pdae, 0.1, 10)) <: SolutionPDAE
