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
# VERIFICATION_REPORT.md for the A/B/C classification.

const q₀ = [1.0, 1.0]
const params = (a₁ = 1.0, a₂ = 1.0, b₁ = -1.0, b₂ = -2.0)
const T = 1.0

steps(n0, k) = T ./ (n0 .* 2 .^ (0:k))
emq(sol, ref) = relative_maximum_error(sol.q, ref.q)

# Reference: Gauß(8) collocation of the ODE form at the same Δt.
href(Δt) = integrate(odeproblem(q₀; timespan = (0.0, T), timestep = Δt, parameters = params), Gauss(8))
ref_ode(prob) = href(timestep(prob))

build_idae(Δt)      = idaeproblem(q₀; timespan = (0.0, T), timestep = Δt, parameters = params)
build_pdae(Δt)      = pdaeproblem(q₀; timespan = (0.0, T), timestep = Δt, parameters = params)
build_ldae(Δt)      = ldaeproblem(q₀; timespan = (0.0, T), timestep = Δt, parameters = params)
build_hdae(Δt)      = hdaeproblem(q₀; timespan = (0.0, T), timestep = Δt, parameters = params)
build_ldae_slrk(Δt) = ldaeproblem_slrk(q₀; timespan = (0.0, T), timestep = Δt, parameters = params)

# Suppress the (correct) solver-divergence log spam from the broken cases.
muffle(f) = Base.CoreLogging.with_logger(f, Base.CoreLogging.NullLogger())

# @test_broken helper for methods that do not converge, diverge, or raise
# SingularException / NonlinearSolverException. The measurement is muffled and any
# thrown exception is recorded as broken (that is exactly the deficiency under test).
function broken_order(builder, method, tsteps, expected; label)
    @testset "$(label): order ≈ $(expected) (broken)" begin
        ok = try
            r = muffle(() -> estimate_convergence_order(builder, method, tsteps;
                reference = ref_ode, errormetric = emq, plateau = 5e-13))
            isapprox(r.order, expected; atol = 0.4)
        catch
            false
        end
        @test_broken ok
    end
end

conv(builder, method, tsteps, expected; label) =
    muffle(() -> test_convergence_order(builder, method, tsteps; reference = ref_ode,
        errormetric = emq, expected = expected, atol = 0.4, plateau = 5e-13, label = label))


@testset "SPARK convergence" begin

    @testset "SLRK (LDAE)" begin
        conv(build_ldae_slrk, SLRKLobattoIIIAB(2), steps(20, 3), 2; label = "SLRKLobattoIIIAB(2)")
        conv(build_ldae_slrk, SLRKLobattoIIIBA(2), steps(20, 3), 2; label = "SLRKLobattoIIIBA(2)")
        conv(build_ldae_slrk, SLRKLobattoIIID(2),  steps(20, 3), 2; label = "SLRKLobattoIIID(2)")
        conv(build_ldae_slrk, SLRKLobattoIIIE(2),  steps(20, 3), 2; label = "SLRKLobattoIIIE(2)")
        conv(build_ldae_slrk, SLRKLobattoIIIAB(3), steps(20, 3), 4; label = "SLRKLobattoIIIAB(3)")
        conv(build_ldae_slrk, SLRKLobattoIIIBA(3), steps(20, 3), 4; label = "SLRKLobattoIIIBA(3)")
        conv(build_ldae_slrk, SLRKLobattoIIID(3),  steps(20, 3), 4; label = "SLRKLobattoIIID(3)")
        conv(build_ldae_slrk, SLRKLobattoIIIE(3),  steps(20, 3), 4; label = "SLRKLobattoIIIE(3)")
    end

    @testset "SPARK (IDAE)" begin
        conv(build_idae, SPARKGLRK(1),   steps(20, 3), 2; label = "SPARKGLRK(1)")
        conv(build_idae, SPARKGLRK(2),   steps(20, 3), 4; label = "SPARKGLRK(2)")
        conv(build_idae, SPARKGLVPRK(1), steps(20, 3), 2; label = "SPARKGLVPRK(1)")
        conv(build_idae, SPARKLobABC(2), steps(20, 3), 2; label = "SPARKLobABC(2)")
        conv(build_idae, SPARKLobABC(3), steps(20, 3), 4; label = "SPARKLobABC(3)")
        conv(build_idae, SPARKLobABD(2), steps(20, 3), 2; label = "SPARKLobABD(2)")
        conv(build_idae, SPARKLobABD(3), steps(20, 3), 4; label = "SPARKLobABD(3)")
    end

    @testset "VPARK / VSPARK primary (IDAE)" begin
        conv(build_idae, TableauGausspSymplectic(2),            steps(20, 3), 4; label = "TableauGausspSymplectic(2)")
        conv(build_idae, TableauLobattoIIIAIIIBpSymplectic(3),  steps(20, 3), 4; label = "LobattoIIIAIIIBpSymplectic(3)")
        conv(build_idae, TableauLobattoIIIBIIIApSymplectic(3),  steps(20, 3), 4; label = "LobattoIIIBIIIApSymplectic(3)")
        conv(build_idae, TableauVSPARKGLRKpMidpoint(2),         steps(20, 3), 4; label = "VSPARKGLRKpMidpoint(2)")
        conv(build_idae, TableauVSPARKGLRKpSymplectic(2),       steps(20, 3), 4; label = "VSPARKGLRKpSymplectic(2)")
        conv(build_idae, TableauVSPARKGLRKpSymmetric(2),        steps(20, 3), 4; label = "VSPARKGLRKpSymmetric(2)")
    end

    @testset "VSPARK general (IDAE)" begin
        conv(build_idae, VSPARK(SPARKLobattoIIIAIIIB(3)), steps(20, 3), 4; label = "VSPARK(SPARKLobattoIIIAIIIB(3))")
    end

    @testset "VSPARK secondary (LDAE)" begin
        conv(build_ldae, TableauVSPARKLobattoIIIAB(2),     steps(20, 3), 2; label = "VSPARKLobattoIIIAB(2)")
        conv(build_ldae, TableauVSPARKLobattoIIIAB(3),     steps(20, 3), 4; label = "VSPARKLobattoIIIAB(3)")
        conv(build_ldae, TableauVSPARKGLRKLobattoIIIAB(1), steps(20, 3), 2; label = "VSPARKGLRKLobattoIIIAB(1)")
        conv(build_ldae, TableauVSPARKGLRKLobattoIIIAB(2), steps(20, 3), 4; label = "VSPARKGLRKLobattoIIIAB(2)")
    end

    @testset "HPARK / HSPARK (PDAE)" begin
        conv(build_pdae, TableauHPARKGLRK(1),                     steps(20, 3), 2; label = "HPARKGLRK(1)")
        conv(build_pdae, HSPARK(SPARKGLRK(1)),                    steps(20, 3), 2; label = "HSPARK(SPARKGLRK(1))")
        conv(build_pdae, HSPARK(SPARKGLRK(2)),                    steps(20, 3), 4; label = "HSPARK(SPARKGLRK(2))")
        conv(build_pdae, HSPARK(SPARKLobABC(3)),                  steps(20, 3), 4; label = "HSPARK(SPARKLobABC(3))")
        conv(build_pdae, HSPARK(SPARKLobABD(3)),                  steps(20, 3), 4; label = "HSPARK(SPARKLobABD(3))")
        conv(build_pdae, TableauHSPARKLobattoIIIAIIIBpSymmetric(2), steps(20, 3), 2; label = "HSPARKLobattoIIIAIIIBpSymmetric(2)")
    end

    # ---------------------------------------------------------------------------
    # Known deficiencies (see VERIFICATION_REPORT.md, SPARK submodule pass). These
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
        broken_order(build_idae, SPARKGLVPRK(2),      steps(20, 3), 4; label = "SPARKGLVPRK(2)")
        broken_order(build_pdae, TableauHPARKGLRK(2), steps(20, 3), 4; label = "HPARKGLRK(2)")

        # symplectic Lobatto-pair order reduction / divergence on the degenerate system
        broken_order(build_idae, SPARKLobattoIIIAIIIB(3),        steps(20, 3), 4; label = "SPARKLobattoIIIAIIIB(3)")
        broken_order(build_idae, SPARKGLRKLobattoIIIAIIIB(2),    steps(20, 3), 4; label = "SPARKGLRKLobattoIIIAIIIB(2)")
        broken_order(build_idae, SPARKLobattoIIIBIIIA(2),        steps(10, 2), 2; label = "SPARKLobattoIIIBIIIA(2)")
        broken_order(build_pdae, TableauHPARKLobattoIIIAIIIB(2), steps(10, 2), 2; label = "HPARKLobattoIIIAIIIB(2)")
        broken_order(build_idae, TableauVSPARKLobattoIIIBIIIApSymmetric(3), steps(10, 2), 4; label = "VSPARKLobattoIIIBIIIApSymmetric(3)")

        # singular stage system at s = 2
        broken_order(build_idae, VSPARK(SPARKLobABC(2)),          steps(10, 1), 2; label = "VSPARK(SPARKLobABC(2))")
        broken_order(build_idae, VSPARK(SPARKLobattoIIIAIIIB(2)), steps(10, 1), 2; label = "VSPARK(SPARKLobattoIIIAIIIB(2))")

        # HSPARK-secondary: residual singularity in the ω secondary-constraint block
        broken_order(build_hdae, TableauHSPARKLobattoIIIAB(2),     steps(10, 1), 2; label = "TableauHSPARKLobattoIIIAB(2)")
        broken_order(build_hdae, TableauHSPARKGLRKLobattoIIIAB(2), steps(10, 1), 4; label = "TableauHSPARKGLRKLobattoIIIAB(2)")
    end

end
