using GeometricIntegrators
using GeometricIntegrators.SPARK
using GeometricProblems.LotkaVolterra2d
using RungeKutta
using Test

include("verification_utilities.jl")

# Dynamic (multi-Δt) convergence verification for the SPARK family, mirroring the
# other test/verification/*_convergence_tests.jl suites. SPARK methods integrate the
# degenerate Lotka–Volterra system in its various DAE formulations; the reference is
# a high-order Gauß collocation of the ODE form at the matching Δt. Expected orders
# are the tableau `o` fields (Gauß 2s, Lobatto-pair 2s−2). Confirmed order
# deficiencies / divergent / singular methods are recorded as @test_broken at the
# documented order with a root-cause tag — see the "SPARK submodule" pass in
# docs/src/audit.md for the A/B/C classification.

const q₀ = [1.0, 1.0]
const params = (a₁ = 1.0, a₂ = 1.0, b₁ = -1.0, b₂ = -2.0)
const T = 1.0

steps(n0, k) = T ./ (n0 .* 2 .^ (0:k))
emq(sol, ref) = relative_maximum_error(sol.q, ref.q)

# Reference: Gauß(8) collocation of the ODE form at the same Δt.
function href(Δt)
    integrate(odeproblem(q₀; timespan = (0.0, T), timestep = Δt, parameters = params), Gauss(8))
end
ref_ode(prob) = href(timestep(prob))

build_idae(Δt) = idaeproblem(q₀; timespan = (0.0, T), timestep = Δt, parameters = params)
build_pdae(Δt) = pdaeproblem(q₀; timespan = (0.0, T), timestep = Δt, parameters = params)
build_ldae(Δt) = ldaeproblem(q₀; timespan = (0.0, T), timestep = Δt, parameters = params)
build_hdae(Δt) = hdaeproblem(q₀; timespan = (0.0, T), timestep = Δt, parameters = params)
function build_ldae_slrk(Δt)
    ldaeproblem_slrk(q₀; timespan = (0.0, T), timestep = Δt, parameters = params)
end

# @test_broken helper for methods that do not converge, diverge, or raise
# SingularException. The measurement runs at `verbosity = 0` to suppress the (correct)
# solver-divergence warnings, and any thrown exception is recorded as broken (that is
# exactly the deficiency under test).
#
# `NonlinearSolverException` no longer reaches this `catch`: since GeometricIntegratorsBase
# 0.6 the time-stepping loop catches it itself, warns naming the timestep, and returns the
# trajectory truncated there — so the deficiency now arrives as a nonsensical convergence
# order rather than as a throw, and is recorded broken either way. That warning is
# GeometricIntegratorsBase's own and `verbosity` does not gate it; six fire here, all from
# the diverging methods below. See the 0.11.0 update in docs/src/audit.md.
function broken_order(builder, method, tsteps, expected; label)
    @testset "$(label): order ≈ $(expected) (broken)" begin
        ok = try
            r = estimate_convergence_order(builder, method, tsteps;
                reference = ref_ode, errormetric = emq, plateau = 5e-13,
                integrate_options = (verbosity = 0, warn_iterations = 0))
            isapprox(r.order, expected; atol = 0.4)
        catch
            false
        end
        @test_broken ok
    end
end

function conv(builder, method, tsteps, expected; label)
    test_convergence_order(builder, method, tsteps; reference = ref_ode,
        errormetric = emq, expected = expected, atol = 0.4, plateau = 5e-13, label = label)
end

@testset "SPARK convergence" begin
    @testset "SLRK (LDAE)" begin
        conv(build_ldae_slrk, SLRKLobattoIIIAB(2),
            steps(20, 3), 2; label = "SLRKLobattoIIIAB(2)")
        conv(build_ldae_slrk, SLRKLobattoIIIBA(2),
            steps(20, 3), 2; label = "SLRKLobattoIIIBA(2)")
        conv(build_ldae_slrk, SLRKLobattoIIICC̄(2),
            steps(20, 3), 2; label = "SLRKLobattoIIICC̄(2)")
        conv(build_ldae_slrk, SLRKLobattoIIIC̄C(2),
            steps(20, 3), 2; label = "SLRKLobattoIIIC̄C(2)")
        conv(build_ldae_slrk, SLRKLobattoIIID(2),
            steps(20, 3), 2; label = "SLRKLobattoIIID(2)")
        conv(build_ldae_slrk, SLRKLobattoIIIE(2),
            steps(20, 3), 2; label = "SLRKLobattoIIIE(2)")
        conv(build_ldae_slrk, SLRKLobattoIIIAB(3),
            steps(20, 3), 4; label = "SLRKLobattoIIIAB(3)")
        conv(build_ldae_slrk, SLRKLobattoIIIBA(3),
            steps(20, 3), 4; label = "SLRKLobattoIIIBA(3)")
        conv(build_ldae_slrk, SLRKLobattoIIICC̄(3),
            steps(20, 3), 4; label = "SLRKLobattoIIICC̄(3)")
        conv(build_ldae_slrk, SLRKLobattoIIIC̄C(3),
            steps(20, 3), 4; label = "SLRKLobattoIIIC̄C(3)")
        conv(build_ldae_slrk, SLRKLobattoIIID(3),
            steps(20, 3), 4; label = "SLRKLobattoIIID(3)")
        conv(build_ldae_slrk, SLRKLobattoIIIE(3),
            steps(20, 3), 4; label = "SLRKLobattoIIIE(3)")
    end

    @testset "SPARK (IDAE)" begin
        conv(build_idae, SPARKGLRK(1), steps(20, 3), 2; label = "SPARKGLRK(1)")
        conv(build_idae, SPARKGLRK(2), steps(20, 3), 4; label = "SPARKGLRK(2)")
        conv(build_idae, SPARKGLVPRK(1), steps(20, 3), 2; label = "SPARKGLVPRK(1)")
        conv(build_idae, SPARKLobABC(2), steps(20, 3), 2; label = "SPARKLobABC(2)")
        conv(build_idae, SPARKLobABC(3), steps(20, 3), 4; label = "SPARKLobABC(3)")
        conv(build_idae, SPARKLobABD(2), steps(20, 3), 2; label = "SPARKLobABD(2)")
        conv(build_idae, SPARKLobABD(3), steps(20, 3), 4; label = "SPARKLobABD(3)")
    end

    @testset "VPARK / VSPARK primary (IDAE)" begin
        conv(build_idae, TableauGausspSymplectic(2), steps(20, 3),
            4; label = "TableauGausspSymplectic(2)")
        conv(build_idae, TableauLobattoIIIAIIIBpSymplectic(3),
            steps(20, 3), 4; label = "LobattoIIIAIIIBpSymplectic(3)")
        conv(build_idae, TableauLobattoIIIBIIIApSymplectic(3),
            steps(20, 3), 4; label = "LobattoIIIBIIIApSymplectic(3)")
        conv(build_idae, TableauVSPARKGLRKpMidpoint(2),
            steps(20, 3), 4; label = "VSPARKGLRKpMidpoint(2)")
        conv(build_idae, TableauVSPARKGLRKpSymplectic(2),
            steps(20, 3), 4; label = "VSPARKGLRKpSymplectic(2)")
        conv(build_idae, TableauVSPARKGLRKpSymmetric(2),
            steps(20, 3), 4; label = "VSPARKGLRKpSymmetric(2)")

        # Both of these were casualties of the R∞ = -1^(s+1) precedence bug, which
        # pinned R∞ to -1 at every s instead of (-1)^(s+1). At s = 3 the correct value
        # is +1, and restoring it lifts AIIIB from order 3.09 to 3.999 and BIIIA from
        # a failed measurement to 4.000. They are the only two Lobatto-pair projections
        # whose order the sign actually moved, so they are the regression guard for it.
        conv(build_idae, TableauVSPARKLobattoIIIAIIIBpSymmetric(3),
            steps(20, 3), 4; label = "VSPARKLobattoIIIAIIIBpSymmetric(3)")
        conv(build_idae, TableauVSPARKLobattoIIIBIIIApSymmetric(3),
            steps(20, 3), 4; label = "VSPARKLobattoIIIBIIIApSymmetric(3)")
    end

    @testset "VSPARK general (IDAE)" begin
        conv(build_idae, VSPARK(SPARKLobattoIIIAIIIB(3)), steps(20, 3),
            4; label = "VSPARK(SPARKLobattoIIIAIIIB(3))")
    end

    @testset "VSPARK secondary (LDAE)" begin
        conv(build_ldae, TableauVSPARKLobattoIIIAB(2),
            steps(20, 3), 2; label = "VSPARKLobattoIIIAB(2)")
        conv(build_ldae, TableauVSPARKLobattoIIIAB(3),
            steps(20, 3), 4; label = "VSPARKLobattoIIIAB(3)")
        conv(build_ldae, TableauVSPARKGLRKLobattoIIIAB(1),
            steps(20, 3), 2; label = "VSPARKGLRKLobattoIIIAB(1)")
        conv(build_ldae, TableauVSPARKGLRKLobattoIIIAB(2),
            steps(20, 3), 4; label = "VSPARKGLRKLobattoIIIAB(2)")
    end

    @testset "HPARK / HSPARK (PDAE)" begin
        conv(build_pdae, TableauHPARKGLRK(1), steps(20, 3), 2; label = "HPARKGLRK(1)")
        conv(build_pdae, HSPARK(SPARKGLRK(1)), steps(20, 3),
            2; label = "HSPARK(SPARKGLRK(1))")
        conv(build_pdae, HSPARK(SPARKGLRK(2)), steps(20, 3),
            4; label = "HSPARK(SPARKGLRK(2))")
        conv(build_pdae, HSPARK(SPARKLobABC(3)), steps(20, 3),
            4; label = "HSPARK(SPARKLobABC(3))")
        conv(build_pdae, HSPARK(SPARKLobABD(3)), steps(20, 3),
            4; label = "HSPARK(SPARKLobABD(3))")
        conv(build_pdae, TableauHSPARKLobattoIIIAIIIBpSymmetric(2),
            steps(20, 3), 2; label = "HSPARKLobattoIIIAIIIBpSymmetric(2)")
    end

    # ---------------------------------------------------------------------------
    # Known deficiencies (see docs/src/audit.md, SPARK submodule pass). These
    # are inherent properties of the methods (category A), not implementation bugs:
    #   * R(∞) ≠ 1 at even σ̃ reduces the GL-VPRK / GL-HPARK order 2s → 2;
    #   * symplectic Lobatto pairs enforcing the constraint at the solution reduce
    #     order or diverge on the degenerate Lagrangian (papers, Sec. 3.4);
    #   * the s = 2 pairs give a degenerate (singular) stage system.
    # HSPARK-secondary Lobatto (below) remains @test_broken for an implementation
    # reason (category B): a residual singularity in the ω secondary-constraint block.
    # ---------------------------------------------------------------------------
    @testset "Known deficiencies (broken)" begin
        # R(∞) = -1 order reduction 2s → 2
        broken_order(build_idae, SPARKGLVPRK(2), steps(20, 3), 4; label = "SPARKGLVPRK(2)")
        broken_order(
            build_pdae, TableauHPARKGLRK(2), steps(20, 3), 4; label = "HPARKGLRK(2)")

        # symplectic Lobatto-pair order reduction / divergence on the degenerate system
        broken_order(build_idae, SPARKLobattoIIIAIIIB(3), steps(20, 3),
            4; label = "SPARKLobattoIIIAIIIB(3)")
        broken_order(build_idae, SPARKGLRKLobattoIIIAIIIB(2), steps(20, 3),
            4; label = "SPARKGLRKLobattoIIIAIIIB(2)")
        broken_order(build_idae, SPARKLobattoIIIBIIIA(2), steps(10, 2),
            2; label = "SPARKLobattoIIIBIIIA(2)")
        broken_order(build_pdae, TableauHPARKLobattoIIIAIIIB(2),
            steps(10, 2), 2; label = "HPARKLobattoIIIAIIIB(2)")

        # singular stage system at s = 2
        broken_order(build_idae, VSPARK(SPARKLobABC(2)), steps(10, 1),
            2; label = "VSPARK(SPARKLobABC(2))")
        broken_order(build_idae, VSPARK(SPARKLobattoIIIAIIIB(2)), steps(10, 1),
            2; label = "VSPARK(SPARKLobattoIIIAIIIB(2))")

        # HSPARK-secondary: residual singularity in the ω secondary-constraint block
        broken_order(build_hdae, TableauHSPARKLobattoIIIAB(2), steps(10, 1),
            2; label = "TableauHSPARKLobattoIIIAB(2)")
        broken_order(build_hdae, TableauHSPARKGLRKLobattoIIIAB(2), steps(10, 1),
            4; label = "TableauHSPARKGLRKLobattoIIIAB(2)")
    end
end
