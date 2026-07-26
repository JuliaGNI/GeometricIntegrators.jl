# Symplecticity of the Degenerate Variational Runge-Kutta (DVRK) method.
#
#     julia --project=scripts scripts/dvrk_symplecticity.jl
#
# DVRK preserves the noncanonical symplectic two-form ω = dϑ of a degenerate
# Lagrangian L = ϑ(q)⋅q̇ - H(q), provided that
#
#   (a) the coefficient matrix A is invertible,
#   (b) the tableau satisfies b_i a_ij + b_j a_ji = b_i b_j,
#   (c) the momentum is initialised with p₀ = ϑ(q₀),
#   (d) ϑ is given in a gauge with ϑ_μ = 0 for μ > d/2,
#   (e) the d/2 × d/2 matrix ∂ϑ_μ/∂q^ν, μ ≤ d/2 < ν, is invertible.
#
# (d) and (e) are independent: for LotkaVolterra2d, which violates (d), the matrix
# of (e) is ∂ϑ₁/∂q² = 1 + 1/q₁q₂ ≠ 0 and hence invertible.
#
# This script measures how far the one-step map q₀ ↦ q₁ is from preserving ω, as
#
#     ‖Jᵀ Ω(q₁) J - Ω(q₀)‖_∞ ,    Ω_μν = ∂ϑ_ν/∂q^μ - ∂ϑ_μ/∂q^ν ,
#
# with J the Jacobian of the map, obtained by central differences. When the
# hypotheses hold the defect sits at the finite-difference floor and does not
# grow with the step size; when (b) or (d) is violated it does.

using GeometricIntegrators
using Printf

import GeometricProblems.LotkaVolterra2d as LV                  # ϑ₂ = q₁ ≠ 0, violates (d)
import GeometricProblems.LotkaVolterra2dSingular as LVSingular  # ϑ₂ = 0, satisfies (d)

const Q₀ = [2.0, 1.5]
const PARAMS = (a₁ = 1.0, a₂ = 1.0, b₁ = -1.0, b₂ = -2.0)
const STEPS = (0.1, 0.5)

# Ω in two dimensions is determined by its single independent entry.
function Ω(m, q)
    o = m.dϑ₂dx₁(0.0, q) - m.dϑ₁dx₂(0.0, q)
    [0.0 -o; o 0.0]
end

function onestep(m, method, q, h)
    prob = m.lodeproblem(collect(q); timespan = (0.0, h), timestep = h, parameters = PARAMS)
    collect(integrate(prob, method).q[end])
end

function symplecticity_defect(m, method, q₀, h; δ = 1e-6)
    f(q) = onestep(m, method, q, h)
    J = zeros(length(q₀), length(q₀))
    for k in eachindex(q₀)
        qp = copy(q₀); qp[k] += δ
        qm = copy(q₀); qm[k] -= δ
        J[:, k] = (f(qp) .- f(qm)) ./ 2δ
    end
    maximum(abs, J' * Ω(m, f(q₀)) * J .- Ω(m, q₀))
end

const CASES = [
    ("LV singular  (in class)",        LVSingular, () -> DVRK(Gauss(1))),
    ("LV singular  (in class)",        LVSingular, () -> DVRK(Gauss(2))),
    ("LV singular  (in class)",        LVSingular, () -> DVRK(Gauss(3))),
    ("LV standard  (violates (d))",    LV,         () -> DVRK(Gauss(2))),
    ("LV singular  (violates (b))",    LVSingular, () -> DVRK(RadauIIA(2); check_conditions = false)),
]

function main()
    println("\nSymplecticity defect ‖JᵀΩ(q₁)J - Ω(q₀)‖ of the one-step DVRK map")
    @printf("q₀ = %s,  central differences with δ = 1e-6,  ‖Ω(q₀)‖ = %.3g\n\n",
            Q₀, maximum(abs, Ω(LVSingular, Q₀)))
    @printf("%-30s %-14s %12s %12s\n", "case", "method", "h = $(STEPS[1])", "h = $(STEPS[2])")
    println("-"^72)
    for (label, m, mk) in CASES
        method = mk()
        name = "$(tableau(method).name)($(tableau(method).s))"
        defects = [symplecticity_defect(m, method, Q₀, h) for h in STEPS]
        @printf("%-30s %-14s %12.2e %12.2e\n", label, name, defects...)
    end
    println("\nA symplectic map gives a defect at the finite-difference floor (~1e-11)")
    println("that does not grow with h. Growth with h indicates a genuine defect.\n")
end

main()
