# Convergence order of the Degenerate Variational Runge-Kutta (DVRK) method,
# inside and outside its hypothesis class.
#
#     julia --project=scripts scripts/dvrk_convergence.jl
#
# DVRK is built for degenerate Lagrangians L = ϑ(q)⋅q̇ - H(q) in which d/2 of the
# components of the symplectic potential ϑ vanish identically. The potential is
# only defined up to an exact one-form, so whether a given problem satisfies that
# hypothesis is a property of the *gauge*, not of the dynamics.
#
# This script measures the empirical order of DVRK(Gauss(s)) for two physical
# systems, each in a gauge that satisfies the hypothesis and in gauges that do
# not. Inside the class the full order 2s is attained; outside it the order drops
# to s. Reproduces the table in docs/src/integrators/dvi.md.

using GeometricIntegrators
using Printf

import GeometricProblems.LotkaVolterra2d as LV                    # ϑ = (q₂ + log q₂/q₁, q₁)
import GeometricProblems.LotkaVolterra2dSingular as LVSingular    # ϑ = (log q₂/q₁, 0)
import GeometricProblems.LotkaVolterra2dSymmetric as LVSymmetric  # ϑ = (log q₂/2q₁, -log q₁/2q₂)
import GeometricProblems.MasslessChargedParticle as MCP           # A₁, A₂ both ≠ 0
import GeometricProblems.MasslessChargedParticleSingular as MCPSingular  # A₂ = 0

const T = 1.0
const LV_Q₀ = [2.0, 1.0]
const LV_PARAMS = (a₁ = 1.0, a₂ = 1.0, b₁ = -1.0, b₂ = -2.0)

steps(n₀, k) = T ./ (n₀ .* 2 .^ (0:k))

"""
Least-squares slope of log(err) against log(Δt), over the points that are above
the roundoff plateau and still strictly decreasing.
"""
function loglog_slope(dts, errs; plateau = 1e-13)
    mask = [errs[i] > plateau && (i == 1 || errs[i] < errs[i - 1]) for i in eachindex(errs)]
    count(mask) < 2 && return NaN
    x, y = log.(dts[mask]), log.(errs[mask])
    x̄, ȳ = sum(x) / length(x), sum(y) / length(y)
    sum((x .- x̄) .* (y .- ȳ)) / sum((x .- x̄) .^ 2)
end

function convergence_order(build, reference, method, dts)
    errs = [relative_maximum_error(integrate(build(Δt), method).q, reference(Δt).q)
            for Δt in dts]
    (order = loglog_slope(collect(dts), errs), errs = errs)
end

# Lotka-Volterra: all gauges share the same dynamics, hence the same reference.
function lv_reference(Δt)
    integrate(
        LV.odeproblem(LV_Q₀; timespan = (0.0, T), timestep = Δt,
            parameters = LV_PARAMS), Gauss(8))
end
function lv_builder(m)
    Δt -> m.lodeproblem(LV_Q₀; timespan = (0.0, T), timestep = Δt,
        parameters = LV_PARAMS)
end

function mcp_reference(Δt)
    integrate(MCP.odeproblem(; timespan = (0.0, T), timestep = Δt), Gauss(8))
end
mcp_builder(m) = Δt -> m.lodeproblem(; timespan = (0.0, T), timestep = Δt)

# Starting resolutions per stage number. Coarser for the higher-stage methods,
# which otherwise hit machine precision before the fit has enough points above the
# roundoff plateau. The in-class cases are 2s-accurate and so saturate several
# refinements earlier than the out-of-class ones, which need a finer start to be
# in the asymptotic regime at s = 3; hence one setting per class rather than one
# global setting.
const N₀_IN = Dict(1 => 20, 2 => 10, 3 => 2)
const N₀_OUT = Dict(1 => 20, 2 => 10, 3 => 5)

const CASES = [
    ("Lotka-Volterra", "ϑ₂ = 0  (in class)", lv_builder(LVSingular), lv_reference, N₀_IN),
    ("Lotka-Volterra", "ϑ₂ = q₁ ≠ 0", lv_builder(LV), lv_reference, N₀_OUT),
    ("Lotka-Volterra", "both components ≠ 0",
        lv_builder(LVSymmetric), lv_reference, N₀_OUT),
    ("Charged particle", "A₂ = 0  (in class)",
        mcp_builder(MCPSingular), mcp_reference, N₀_IN),
    ("Charged particle", "both components ≠ 0", mcp_builder(MCP), mcp_reference, N₀_OUT)
]

function main()
    println("\nDVRK convergence order, T = $T, reference Gauss(8) on the ODE")
    println("Expected: 2s inside the hypothesis class, s outside it.\n")
    @printf("%-18s %-24s %8s %8s %8s\n", "problem", "gauge", "s = 1", "s = 2", "s = 3")
    println("-"^72)
    for (problem, gauge, build, reference, n₀) in CASES
        orders = [convergence_order(build, reference, DVRK(Gauss(s)), steps(n₀[s], 4)).order
                  for s in 1:3]
        @printf("%-18s %-24s %8.2f %8.2f %8.2f\n", problem, gauge, orders...)
    end
    println()
end

main()
