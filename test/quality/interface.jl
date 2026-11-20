using GeometricIntegrators
using GeometricIntegrators.SPARK
using GeometricBase: GeometricBase, AbstractMethod
using CompactBasisFunctions
using QuadratureRules
using InteractiveUtils: subtypes
using Test

# `GeometricBase.interface_violations` lists the properties of a method that throw, answer the
# wrong type, or are shadowed by a function of this package that does not extend GeometricBase's.
# Every method type of this package with a constructor of no argument is checked on that instance,
# and the method types without one on the instances below, built as the other tests build them.

function owned_method_types(T = AbstractMethod, types = Set{Any}())
    for S in subtypes(T)
        if !isabstracttype(S) && Base.moduleroot(parentmodule(S)) === GeometricIntegrators
            push!(types, S)
        end
        owned_method_types(S, types)
    end
    return types
end

const OWNED_METHOD_TYPES = sort!(collect(owned_method_types()); by = string)
const NULLARY_METHOD_TYPES = filter(T -> hasmethod(T, Tuple{}), OWNED_METHOD_TYPES)

QGau4 = GaussLegendreQuadrature(4)
BGau4 = Lagrange(QuadratureRules.nodes(QGau4))
QLob4 = LobattoLegendreQuadrature(4)
BLob4 = Lagrange(QuadratureRules.nodes(QLob4))
dgjump = Discontinuity(PathIntegralLinear(), LobattoLegendreQuadrature(2))

ϕ(v, q̄, q, a, Δt) = v .= (q .- q̄) ./ Δt
Dϕ(d, q̄, q, a, Δt) = d

const METHOD_INSTANCES = [
    # Runge-Kutta and partitioned Runge-Kutta methods
    "Gauss(2)" => Gauss(2),
    "RadauIA(2)" => RadauIA(2),
    "RadauIB(2)" => RadauIB(2),
    "RadauIIA(2)" => RadauIIA(2),
    "RadauIIB(2)" => RadauIIB(2),
    "LobattoIIIA(2)" => LobattoIIIA(2),
    "LobattoIIIC(2)" => LobattoIIIC(2),
    "LobattoIIIAIIIB(2)" => LobattoIIIAIIIB(2),
    "LobattoIII(2)" => LobattoIII(2),
    "LobattoIIIAIIIĀ(2)" => LobattoIIIAIIIĀ(2),
    "LobattoIIIB(2)" => LobattoIIIB(2),
    "LobattoIIIBIIIA(2)" => LobattoIIIBIIIA(2),
    "LobattoIIIBIIIB̄(2)" => LobattoIIIBIIIB̄(2),
    "LobattoIIICIIIC̄(2)" => LobattoIIICIIIC̄(2),
    "LobattoIIIC̄IIIC(2)" => LobattoIIIC̄IIIC(2),
    "LobattoIIID(2)" => LobattoIIID(2),
    "LobattoIIIDIIID̄(2)" => LobattoIIIDIIID̄(2),
    "LobattoIIIE(2)" => LobattoIIIE(2),
    "LobattoIIIEIIIĒ(2)" => LobattoIIIEIIIĒ(2),
    "LobattoIIIF(2)" => LobattoIIIF(2),
    "LobattoIIIFIIIF̄(2)" => LobattoIIIFIIIF̄(2),
    "LobattoIIIF̄(2)" => LobattoIIIF̄(2),
    "LobattoIIIF̄IIIF(2)" => LobattoIIIF̄IIIF(2),
    "LobattoIIIG(2)" => LobattoIIIG(2),
    "LobattoIIIGIIIḠ(2)" => LobattoIIIGIIIḠ(2),
    "PartitionedGauss(2)" => PartitionedGauss(2),
    "ERK(RK4())" => ERK(RK4()),
    "DIRK(Crouzeix())" => DIRK(Crouzeix()),
    "IRK(Gauss(2); implicit_update = true)" => IRK(Gauss(2); implicit_update = true),
    "EPRK(LobattoIIIAIIIB(2))" => EPRK(LobattoIIIAIIIB(2)),
    "IPRK(SymplecticEulerARK())" => IPRK(SymplecticEulerARK()),
    "IPRK(PartitionedGauss(2))" => IPRK(PartitionedGauss(2)),
    "FLRK(Gauss(2))" => FLRK(Gauss(2)),
    "PGLRK(3)" => PGLRK(3),
    "DVRK(Gauss(2))" => DVRK(Gauss(2)),
    # variational partitioned Runge-Kutta methods
    "VPRK(Gauss(2))" => VPRK(Gauss(2)),
    "VPRKGauss(2)" => VPRKGauss(2),
    "VPRKLobattoIIIAIIIB(2)" => VPRKLobattoIIIAIIIB(2),
    "VPRKRadauIIA(2)" => VPRKRadauIIA(2),
    "VPRKpTableau(4)" => VPRKpTableau(4),
    "VPRKLobattoIII(2)" => VPRKLobattoIII(2),
    "VPRKLobattoIIIA(2)" => VPRKLobattoIIIA(2),
    "VPRKLobattoIIIAIIIĀ(2)" => VPRKLobattoIIIAIIIĀ(2),
    "VPRKLobattoIIIB(2)" => VPRKLobattoIIIB(2),
    "VPRKLobattoIIIBIIIA(2)" => VPRKLobattoIIIBIIIA(2),
    "VPRKLobattoIIIBIIIB̄(2)" => VPRKLobattoIIIBIIIB̄(2),
    "VPRKLobattoIIIC(2)" => VPRKLobattoIIIC(2),
    "VPRKLobattoIIICIIIC̄(2)" => VPRKLobattoIIICIIIC̄(2),
    "VPRKLobattoIIIC̄IIIC(2)" => VPRKLobattoIIIC̄IIIC(2),
    "VPRKLobattoIIID(2)" => VPRKLobattoIIID(2),
    "VPRKLobattoIIIDIIID̄(2)" => VPRKLobattoIIIDIIID̄(2),
    "VPRKLobattoIIIE(2)" => VPRKLobattoIIIE(2),
    "VPRKLobattoIIIEIIIĒ(2)" => VPRKLobattoIIIEIIIĒ(2),
    "VPRKLobattoIIIF(2)" => VPRKLobattoIIIF(2),
    "VPRKLobattoIIIFIIIF̄(2)" => VPRKLobattoIIIFIIIF̄(2),
    "VPRKLobattoIIIF̄(2)" => VPRKLobattoIIIF̄(2),
    "VPRKLobattoIIIF̄IIIF(2)" => VPRKLobattoIIIF̄IIIF(2),
    "VPRKLobattoIIIG(2)" => VPRKLobattoIIIG(2),
    "VPRKLobattoIIIGIIIḠ(2)" => VPRKLobattoIIIGIIIḠ(2),
    "VPRKRadauIIB(2)" => VPRKRadauIIB(2),
    # Galerkin variational integrators
    "CGVI(BGau4, QGau4)" => CGVI(BGau4, QGau4),
    "CGVINodal(BLob4, QLob4)" => CGVINodal(BLob4, QLob4),
    "DGVI(BGau4, QGau4)" => DGVI(BGau4, QGau4),
    "DGVIP0(BGau4, QGau4)" => DGVIP0(BGau4, QGau4),
    "DGVIP1(BGau4, QGau4)" => DGVIP1(BGau4, QGau4),
    "DGVIEXP(BGau4, QGau4)" => DGVIEXP(BGau4, QGau4),
    "DGVIPI(BGau4, QGau4, dgjump)" => DGVIPI(BGau4, QGau4, dgjump),
    # Hamilton-Pontryagin integrators
    "HPImidpoint(ϕ, Dϕ, Dϕ, Dϕ, Float64[])" => HPImidpoint(ϕ, Dϕ, Dϕ, Dϕ, Float64[]),
    "HPItrapezoidal(ϕ, Dϕ, Dϕ, Dϕ, Float64[])" => HPItrapezoidal(ϕ, Dϕ, Dϕ, Dϕ, Float64[]),
    # splitting and composition methods
    "Composition(Strang())" => Composition(Strang()),
    "Splitting([1, 2], [0.5, 0.5])" => Splitting([1, 2], [0.5, 0.5]),
    # projection methods
    "StandardProjection([0, 1], [0, 1])" => StandardProjection([0, 1], [0, 1]),
    "PostProjection(Gauss(1))" => PostProjection(Gauss(1)),
    "MidpointProjection(Gauss(1))" => MidpointProjection(Gauss(1)),
    "SymmetricProjection(Gauss(1))" => SymmetricProjection(Gauss(1)),
    # SPARK methods
    "SPARKGLRK(1)" => SPARKGLRK(1),
    "TableauHPARKGLRK(1)" => TableauHPARKGLRK(1),
    "TableauGausspSymplectic(1)" => TableauGausspSymplectic(1),
    "SLRKLobattoIIIAB(2)" => SLRKLobattoIIIAB(2),
    "HSPARK(SPARKGLRK(1))" => HSPARK(SPARKGLRK(1)),
    "VSPARK(SPARKLobattoIIIAIIIB(2))" => VSPARK(SPARKLobattoIIIAIIIB(2)),
    "TableauHSPARKGLRKpSymmetric(1)" => TableauHSPARKGLRKpSymmetric(1),
    "TableauHSPARKGLRKLobattoIIIAB(1)" => TableauHSPARKGLRKLobattoIIIAB(1),
    "TableauVSPARKGLRKpMidpoint(1)" => TableauVSPARKGLRKpMidpoint(1),
    "TableauVSPARKLobattoIIIAB(2)" => TableauVSPARKLobattoIIIAB(2)
]

@testset "Interface conformance" begin
    @test !isempty(NULLARY_METHOD_TYPES)

    # Each nullary method is built once, so the coverage check below counts the instances checked here.
    nullary_methods = [T() for T in NULLARY_METHOD_TYPES]

    @testset "$T()" for (T, method) in zip(NULLARY_METHOD_TYPES, nullary_methods)
        @test isempty(GeometricBase.interface_violations(method))
    end

    @testset "$label" for (label, method) in METHOD_INSTANCES
        @test isempty(GeometricBase.interface_violations(method))
    end

    # A new method type without a constructor of no argument and without an instance above fails here.
    @testset "every owned method type is covered" begin
        checked = vcat(nullary_methods, last.(METHOD_INSTANCES))
        @testset "$T" for T in OWNED_METHOD_TYPES
            @test any(method -> method isa T, checked)
        end
    end
end
