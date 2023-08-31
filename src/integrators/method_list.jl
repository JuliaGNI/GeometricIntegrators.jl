
meta_methods = (
    ERK,
    IRK,
    DIRK,
    EPRK,
    IPRK,
    # FLRK,
    VPRK,
    DVRK,
    ProjectedMethod,
    HPImidpoint,
    HPItrapezoidal,
)

explicit_rungekutta_methods = (
    ForwardEuler,
    ExplicitEulerRK,
    ExplicitMidpoint,
    Heun2,
    Heun3,
    Kutta3,
    Ralston2,
    Ralston3,
    RK21,
    RK22,
    RK31,
    RK32,
    RK4,
    RK41,
    RK42,
    RK416,
    RK438,
    RK5,
    Runge2,
    SSPRK2,
    SSPRK3,
)

diagonally_implicit_rungekutta_methods = (
    CrankNicolson,
    Crouzeix,
    KraaijevangerSpijker,
    QinZhang,
)

fully_implicit_rungekutta_methods = (
    BackwardEuler,
    ImplicitEulerRK,
    ImplicitMidpoint,
    SRK3,
)

implicit_rungekutta_methods = (
    diagonally_implicit_rungekutta_methods...,
    fully_implicit_rungekutta_methods...,
)

runge_kutta_methods = (
    explicit_rungekutta_methods...,
    implicit_rungekutta_methods...,
)

runge_kutta_families = (
    Gauss,
    LobattoIII,
    LobattoIIIA,
    LobattoIIIB,
    LobattoIIIC,
    LobattoIIID,
    LobattoIIIE,
    LobattoIIIF,
    LobattoIIIF̄,
    LobattoIIIG,
    RadauIA,
    RadauIB,
    RadauIIA,
    RadauIIB,
)

partitioned_runge_kutta_methods = (
    SymplecticEulerA,
    SymplecticEulerB,
)

partitioned_runge_kutta_families = (
    PartitionedGauss,
    LobattoIIIAIIIB,
    LobattoIIIBIIIA,
    LobattoIIIAIIIĀ,
    LobattoIIIBIIIB̄,
    LobattoIIICIIIC̄,
    LobattoIIIC̄IIIC,
    LobattoIIIDIIID̄,
    LobattoIIIEIIIĒ,
    LobattoIIIFIIIF̄,
    LobattoIIIF̄IIIF,
    LobattoIIIGIIIḠ,
)

variational_partitioned_runge_kutta_families = (
    VPSRK3,
    VPRKGauss,
    VPRKRadauIIA,
    VPRKRadauIIB,
    VPRKLobattoIII,
    VPRKLobattoIIIA,
    VPRKLobattoIIIB,
    VPRKLobattoIIIC,
    VPRKLobattoIIID,
    VPRKLobattoIIIE,
    VPRKLobattoIIIF,
    VPRKLobattoIIIF̄,
    VPRKLobattoIIIG,
    VPRKLobattoIIIAIIIB,
    VPRKLobattoIIIBIIIA,
    VPRKLobattoIIIAIIIĀ,
    VPRKLobattoIIIBIIIB̄,
    VPRKLobattoIIICIIIC̄,
    VPRKLobattoIIIC̄IIIC,
    VPRKLobattoIIIDIIID̄,
    VPRKLobattoIIIEIIIĒ,
    VPRKLobattoIIIFIIIF̄,
    VPRKLobattoIIIF̄IIIF,
    VPRKLobattoIIIGIIIḠ,
)

variational_integrators = (
    PMVImidpoint,
    PMVItrapezoidal,    
)

degenerate_variational_integrators = (
    DVIA,
    DVIB,
    CMDVI,
    CTDVI,
)

splitting_methods = (
    LieA,
    LieB,
    Strang,
    Marchuk,
    StrangA,
    StrangB,
    McLachlan2,
    McLachlan4,
    TripleJump,
    SuzukiFractal,
)

method_groups = (
    meta_methods,
    runge_kutta_methods,
    runge_kutta_families,
    partitioned_runge_kutta_methods,
    partitioned_runge_kutta_families,
    variational_partitioned_runge_kutta_families,
    variational_integrators,
    degenerate_variational_integrators,
    splitting_methods,
)

@inline tuplejoin(x) = x
@inline tuplejoin(x, y) = (x..., y...)
@inline tuplejoin(x, y, z...) = tuplejoin(tuplejoin(x, y), z...)

method_list = tuplejoin(method_groups...)

for m in nameof.(method_list)
    @eval export $m
end

# projection methods

export ProjectionMethod
export ProjectedMethod

# export InternalStageProjection
# export LegendreProjection
export MidpointProjection
export PostProjection
# export SecondaryProjection
export StandardProjection
export SymmetricProjection
export SymplecticProjection
export VariationalProjection
export VariationalProjectionOnP
export VariationalProjectionOnQ

# projected VPRK methods
# export VPRKpInternal
# export VPRKpLegendre
# export VPRKpMidpoint
# export VPRKpSecondary
# export VPRKpStandard
# export VPRKpSymmetric
# export VPRKpSymplectic
# export VPRKpVariational
# export VPRKpVariationalP
# export VPRKpVariationalQ


_display_property(p::Bool) = p ? "✓" : "✗"
_display_property(p::Missing) = p

_true(m) = true

function _row(m, refs, selector)
    row = [
        refs ? "[`$(m)`](@ref)" : string(m),
        string(order(m)),
        _display_property(isexplicit(m)),
        _display_property(issymmetric(m)),
        _display_property(issymplectic(m)),
    ]

    if selector == _true
        row = [row...,
            _display_property(isodemethod(m)),
            _display_property(ispodemethod(m)),
            _display_property(ishodemethod(m)),
            _display_property(isiodemethod(m)),
            _display_property(islodemethod(m)),
            _display_property(issodemethod(m)),
            _display_property(isdaemethod(m)),
            _display_property(ispdaemethod(m)),
            _display_property(ishdaemethod(m)),
            _display_property(isidaemethod(m)),
            _display_property(isldaemethod(m)),
        ]
    end
    
    return row
end

struct MethodList{MD}
    header
    data

    function MethodList(list::Tuple = method_list; markdown::Bool = false, selector = _true, refs = false)
        header = [
            "Method",
            "Order",
            "Explicit",
            "Symmetric",
            "Symplectic",
        ]

        if selector == _true
            header = [header...,
            "ODE",
            "PODE",
            "HODE",
            "IODE",
            "LODE",
            "SODE",
            "DAE",
            "PDAE",
            "HDAE",
            "IDAE",
            "LDAE",
        ]
        end
    
        data = []
        
        for m in list
            if selector(m)
                data = [data..., _row(m, refs, selector)]
            end
        end

        new{markdown}(header, permutedims(hcat(data...)))
    end
end

Base.length(ml::MethodList) = size(ml.data, 1)

"""
```julia
Base.show(io::IO, ml::MethodList)
Base.show(io::IO, ::MIME"text/markdown", ml::MethodList)
```

Pretty-print MethodList.
"""
function Base.show(io::IO, ml::MethodList)
    pretty_table(io, ml.data; header = ml.header, limit_printing = false)
end

function Base.show(io::IO, ml::MethodList{true})
    pretty_table(io, ml.data; header = ml.header, limit_printing = false, tf = tf_markdown)
end

function Base.show(io::IO, ::MIME"text/markdown", ml::MethodList)
    pretty_table(io, ml.data; header = ml.header, limit_printing = false, tf = tf_markdown)
end
