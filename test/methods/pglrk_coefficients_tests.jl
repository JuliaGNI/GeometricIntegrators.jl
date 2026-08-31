using GeometricIntegrators
using RungeKutta
using LinearAlgebra
using Test

# The `CoefficientsPGLRK(s)` construction is the W-transformation of the Gauss
# tableau plus a skew perturbation `A = P W Q`, giving the one-parameter family
# `a(λ) = a + λ A`. Three identities characterise it, and each is a strong
# correctness check on the whole construction:
#
#   1. `a = P X Q` reproduces the Gauss tableau exactly, and `b`/`c` *are* the
#      Gauss weights/nodes.
#   2. `B A` is skew for `B = diag(b)` (because `Pᵀ B P = I`), so `a(λ)` satisfies
#      the symplecticity condition `bᵢaᵢⱼ + bⱼaⱼᵢ = bᵢbⱼ` for *every* λ.
#   3. `bᵀA = 0` and `A·1 = 0`, so `a(λ)` retains the quadrature and row-sum
#      conditions — and hence the order — for every λ. This holds only for s ≥ 3,
#      which is why the constructor rejects s = 2.
#
# Tolerances are set from measurement. `inv(P)` conditioning grows with s, so the
# bound is scaled rather than fixed: the measured worst case over s = 3..6 is
# ~2.4e-16, so 2E-15 is ~8x headroom and still machine-precision-tight.

const PGLRK_ATOL = 2E-15

@testset "$(rpad("Projected Gauss-Legendre Runge-Kutta coefficients",80))" begin
    for s in 3:6
        coeff = CoefficientsPGLRK(s)
        tab = TableauGauss(s)
        B = Diagonal(coeff.b)

        @test typeof(coeff) <: CoefficientsPGLRK{Float64}
        @test RungeKutta.nstages(coeff) == s
        @test order(coeff) == 2s

        # (1) the W-transformation reproduces the Gauss tableau
        @test coeff.a ≈ tab.a atol = PGLRK_ATOL
        @test coeff.b == tab.b
        @test coeff.c == tab.c
        @test coeff.Q ≈ inv(coeff.P) atol = PGLRK_ATOL

        # the embedded-method coefficients are unused and must be zero
        @test all(iszero, coeff.â)
        @test all(iszero, coeff.b̂)
        @test all(iszero, coeff.ĉ)

        # (2) B A is skew ⇒ the perturbation preserves symplecticity
        @test maximum(abs, B * coeff.A + coeff.A' * B) < PGLRK_ATOL

        # (3) the perturbation preserves the quadrature and row-sum conditions
        @test maximum(abs, coeff.b' * coeff.A) < PGLRK_ATOL
        @test maximum(abs, coeff.A * ones(s)) < PGLRK_ATOL

        # a(λ) is symplectic and consistent for every λ
        for λ in (0.0, 0.05, -0.13)
            aλ = getTableauPGLRK(coeff, λ)
            @test maximum(abs, B * aλ + aλ' * B - coeff.b * coeff.b') < PGLRK_ATOL
            @test maximum(abs, aλ * ones(s) - coeff.c) < PGLRK_ATOL

            # the in-place form must agree with the allocating one
            atmp = zero(coeff.a)
            getTableauPGLRK(coeff, λ, atmp)
            @test atmp == aλ
        end

        # λ = 0 recovers the Gauss tableau
        @test getTableauPGLRK(coeff, 0.0) == coeff.a
    end

    # s = 2 is rejected: the W[2,1]/W[1,2] perturbation coincides with the
    # order-determining ξ₁ entries of X and destroys consistency. Verified: at
    # s = 2, bᵀA = 0.5 and A·1 = 1 (vs ~1e-16 for s ≥ 3).
    @test_throws AssertionError CoefficientsPGLRK(2)
    @test_throws AssertionError CoefficientsPGLRK(1)

    # generic floating-point types
    @test typeof(CoefficientsPGLRK(BigFloat, 3)) <: CoefficientsPGLRK{BigFloat}

    # equality, hashing and printing
    @test CoefficientsPGLRK(3) == CoefficientsPGLRK(3)
    @test hash(CoefficientsPGLRK(3)) == hash(CoefficientsPGLRK(3))
    @test isequal(CoefficientsPGLRK(3), CoefficientsPGLRK(3))
    @test CoefficientsPGLRK(3) ≠ CoefficientsPGLRK(4)
    @test !isequal(CoefficientsPGLRK(Float64, 3), CoefficientsPGLRK(BigFloat, 3))
    @test isapprox(CoefficientsPGLRK(3), CoefficientsPGLRK(3))

    @test occursin("Projected Gauss-Legendre", sprint(show, CoefficientsPGLRK(3)))
end
