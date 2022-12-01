
methods = (
    RK,
    PRK,
    # explicit Runge-Kutta methods
    ForwardEuler,
    ExplicitEuler,
    ExplicitMidpoint,
    Heun2,
    Heun3,
    Kutta3,
    Ralston2,
    Ralston3,
    RK4,
    RK416,
    RK438,
    Runge2,
    SSPRK2,
    SSPRK3,
    # diagonally implicit Runge-Kutta methods
    CrankNicolson,
    Crouzeix,
    KraaijevangerSpijker,
    QinZhang,
    # fully implicit Runge-Kutta methods
    BackwardEuler,
    ImplicitEuler,
    ImplicitMidpoint,
    SRK3,
    # Runge-Kutta methods with variable number of stages
    Gauss,
    LobattoIIIA,
    LobattoIIIB,
    LobattoIIIC,
    LobattoIIIC̄,
    LobattoIIID,
    LobattoIIIE,
    LobattoIIIF,
    LobattoIIIF̄,
    LobattoIIIG,
    RadauIA,
    RadauIB,
    RadauIIA,
    RadauIIB,
    # partitioned Runge-Kutta methods with variable number of stages
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
    # variational partitioned Runge-Kutta methods
    VPRK,
    VPSRK3,
    VPRKGauss,
    VPRKRadauIIA,
    VPRKRadauIIB,
    VPRKLobattoIIIA,
    VPRKLobattoIIIB,
    VPRKLobattoIIIC,
    VPRKLobattoIIIC̄,
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
    # projected VPRK methods
    ProjectedVPRK,
    # VPRKpInternal,
    # VPRKpLegendre,
    # VPRKpMidpoint,
    # VPRKpSecondary,
    # VPRKpStandard,
    # VPRKpSymmetric,
    # VPRKpSymplectic,
    # VPRKpVariational,
    # VPRKpVariationalP,
    # VPRKpVariationalQ,
    # degenerate VPRK methods
    DegenerateVPRK,
    # degenerate VI methods
    DVIA,
    DVIB,
    CMDVI,
    CTDVI,
    # splitting methods
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

for m in nameof.(methods)
    @eval export $m
end


_display_property(p::Bool) = p ? "✓" : "✗"
_display_property(p::Missing) = p


struct MethodList{MD}
    header
    data

    function MethodList(list::Tuple = methods; markdown::Bool = false)
        header = [
            "Method",
            "Order",
            "Explicit",
            "Symmetric",
            "Symplectic",
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
    
        data = [[
            string(m),
            string(order(m)),
            _display_property(isexplicit(m)),
            _display_property(issymmetric(m)),
            _display_property(issymplectic(m)),
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
        ] for m in list]

        new{markdown}(header, permutedims(hcat(data...)))
    end
end


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
