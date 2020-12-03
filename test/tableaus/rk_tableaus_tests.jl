
@testset "$(rpad("Runge-Kutta Tableaus",80))" begin

    # explicit midpoint
    a = Rational{Int64}[[0    0]
                        [1//2 0]]
    b = Rational{Int64}[0, 1   ]
    c = Rational{Int64}[0, 1//2]
    o = 2

    @test TableauERK(:explicit_midpoint, o, a, b, c) == TableauExplicitMidpoint()

    tab_explicit_midpoint1 = @test_nowarn TableauERK(:explicit_midpoint, o, a, b, c)
    tab_explicit_midpoint2 = @test_logs (:warn, r"Initializing TableauDIRK with explicit tableau explicit_midpoint.*") TableauDIRK(:explicit_midpoint, o, a, b, c)
    tab_explicit_midpoint3 = @test_logs (:warn, r"Initializing TableauFIRK with explicit tableau explicit_midpoint.*") TableauFIRK(:explicit_midpoint, o, a, b, c)

    @test tab_explicit_midpoint1 == tab_explicit_midpoint2
    @test tab_explicit_midpoint2 == tab_explicit_midpoint3
    @test tab_explicit_midpoint3 == tab_explicit_midpoint1

    @test !isequal(tab_explicit_midpoint1, tab_explicit_midpoint2)
    @test !isequal(tab_explicit_midpoint2, tab_explicit_midpoint3)
    @test !isequal(tab_explicit_midpoint3, tab_explicit_midpoint1)

    tab_explicit_midpoint1 = TableauERK(:explicit_midpoint, o, a, b, c)
    tab_explicit_midpoint2 = TableauERK(:explicit_midpoint, o, a, b, c)
    @test tab_explicit_midpoint1 == tab_explicit_midpoint2
    @test isequal(tab_explicit_midpoint1, tab_explicit_midpoint2)

    tmp = mktempdir()
    tab_explicit_midpoint1 = TableauERK(:explicit_midpoint, o, a, b, c)
    writeTableauToFile(tmp, tab_explicit_midpoint1)
    tab_explicit_midpoint2 = readTableauERKFromFile(tmp, "explicit_midpoint")
    rm(tmp, recursive=true)
    @test tab_explicit_midpoint1 == tab_explicit_midpoint2
    @test !isequal(tab_explicit_midpoint1, tab_explicit_midpoint2) # fails because data arrays in tab_explicit_midpoint1 are of type Rational{Int64} while data arrays in tab_explicit_midpoint2 are of type Float


    # Heun
    a = [[0.0 0.0]
         [1.0 0.0]]
    b = [0.5, 0.5]
    c = [0.0, 1.0]
    o = 2

    @test_nowarn tab_explicit_heun = TableauERK(:heun, o, a, b, c)
    @test_logs (:warn, r"Initializing TableauDIRK with explicit tableau heun.*") TableauDIRK(:heun, o, a, b, c)
    @test_logs (:warn, r"Initializing TableauFIRK with explicit tableau heun.*") TableauFIRK(:heun, o, a, b, c)

    @test TableauERK(:heun, o, a, b, c) == TableauHeun2()

    tmp = mktempdir()
    tab_explicit_heun1 = TableauERK(:heun, o, a, b, c)
    writeTableauToFile(tmp, tab_explicit_heun1)
    tab_explicit_heun2 = readTableauERKFromFile(tmp, "heun")
    rm(tmp, recursive=true)
    @test tab_explicit_heun1 == tab_explicit_heun2
    @test isequal(tab_explicit_heun1, tab_explicit_heun2)


    # Crouzeix
    fac = 0.5/√3
    a = [[ 0.5+fac 0.0    ]
         [-2.0*fac 0.5+fac]]
    b = [0.5,     0.5    ]
    c = [0.5+fac, 0.5-fac]
    o = 3

    @test TableauDIRK(:crouzeix, o, a, b, c) == TableauCrouzeix()

    @test_throws AssertionError tab_explicit_crouzeix = TableauERK(:crouzeix, o, a, b, c)
    @test_nowarn tab_explicit_crouzeix = TableauDIRK(:crouzeix, o, a, b, c)
    @test_logs (:warn, r"Initializing TableauFIRK with diagonally implicit tableau crouzeix.*") tab_explicit_crouzeix = TableauFIRK(:crouzeix, o, a, b, c)


    # implicit midpoint
    a = ones(Rational{Int64}, 1, 1) // 2
    b = Rational{Int64}[1   ]
    c = Rational{Int64}[1//2]
    o = 2

    @test TableauFIRK(:implicit_midpoint, o, a, b, c) == TableauImplicitMidpoint()
    @test TableauFIRK(:implicit_midpoint, o, a, b, c) == TableauGLRK(1)

    @test_throws AssertionError tab_implicit_midpoint = TableauERK(:implicit_midpoint, o, a, b, c)
    @test_throws AssertionError tab_implicit_midpoint = TableauDIRK(:implicit_midpoint, o, a, b, c)
    tab_implicit_midpoint = TableauFIRK(:implicit_midpoint, o, a, b, c)


    # Gauss-Legendre
    a = Array{Float64}(@dec128 [
         [0.25       0.25-√3/6]
         [0.25+√3/6  0.25     ]
        ])
    b = Array{Float64}(@dec128 [0.5,      0.5     ])
    c = Array{Float64}(@dec128 [0.5-√3/6, 0.5+√3/6])
    o = 4

    @test TableauFIRK(:glrk2, o, a, b, c) == TableauGLRK(2)

    @test_throws AssertionError tab_explicit_glrk2 = TableauERK(:glrk2, o, a, b, c)
    @test_throws AssertionError tab_explicit_glrk2 = TableauDIRK(:glrk2, o, a, b, c)
    tab_explicit_glrk2 = TableauFIRK(:glrk2, o, a, b, c)


    # instantiate all explicit tableaus
    @test typeof(TableauExplicitEuler())        <: TableauERK
    @test typeof(TableauForwardEuler())         <: TableauERK
    @test typeof(TableauExplicitMidpoint())     <: TableauERK
    @test typeof(TableauHeun2())                <: TableauERK
    @test typeof(TableauHeun3())                <: TableauERK
    @test typeof(TableauRalston2())             <: TableauERK
    @test typeof(TableauRalston3())             <: TableauERK
    @test typeof(TableauRunge2())               <: TableauERK
    @test typeof(TableauRunge())                <: TableauERK
    @test typeof(TableauKutta3())               <: TableauERK
    @test typeof(TableauKutta())                <: TableauERK
    @test typeof(TableauERK4())                 <: TableauERK
    @test typeof(TableauERK416())               <: TableauERK
    @test typeof(TableauERK438())               <: TableauERK

    # instantiate all diagonally implicit tableaus
    @test typeof(TableauCrankNicolson())        <: TableauDIRK
    @test typeof(TableauKraaijevangerSpijker()) <: TableauDIRK
    @test typeof(TableauQinZhang())             <: TableauDIRK
    @test typeof(TableauCrouzeix())             <: TableauDIRK

    # instantiate all fully implicit tableaus
    @test typeof(TableauImplicitEuler())        <: TableauFIRK
    @test typeof(TableauBackwardEuler())        <: TableauFIRK
    @test typeof(TableauImplicitMidpoint())     <: TableauFIRK
    @test typeof(TableauSRK3())                 <: TableauFIRK

    @test typeof(TableauGLRK(1))                <: TableauFIRK
    @test typeof(TableauGLRK(2))                <: TableauFIRK

    @test typeof(TableauRadauIA(2))             <: TableauFIRK
    @test typeof(TableauRadauIA(3))             <: TableauFIRK
    @test typeof(TableauRadauIA(4))             <: TableauFIRK

    @test typeof(TableauRadauIIA(2))            <: TableauFIRK
    @test typeof(TableauRadauIIA(3))            <: TableauFIRK
    @test typeof(TableauRadauIIA(4))            <: TableauFIRK

    @test typeof( @test_logs (:warn, r"Initializing TableauFIRK with diagonally implicit tableau LobattoIIIA\(2\).*") TableauLobattoIIIA(2) ) <: TableauFIRK
    @test typeof(TableauLobattoIIIA(3)) <: TableauFIRK
    @test typeof(TableauLobattoIIIA(4)) <: TableauFIRK
    @test typeof( @test_logs (:warn, r"Initializing TableauFIRK with explicit tableau LobattoIIIB\(2\).*") TableauLobattoIIIB(2) ) <: TableauFIRK
    @test typeof(TableauLobattoIIIB(3)) <: TableauFIRK
    @test typeof(TableauLobattoIIIB(4)) <: TableauFIRK
    @test typeof(TableauLobattoIIIC(2)) <: TableauFIRK
    @test typeof(TableauLobattoIIIC(3)) <: TableauFIRK
    @test typeof(TableauLobattoIIIC(4)) <: TableauFIRK
    @test typeof( @test_logs (:warn, r"Initializing TableauFIRK with explicit tableau LobattoIIIC̄\(2\).*") TableauLobattoIIIC̄(2) ) <: TableauFIRK
    @test typeof( @test_logs (:warn, r"Initializing TableauFIRK with diagonally implicit tableau LobattoIIIC̄\(3\).*") TableauLobattoIIIC̄(3) ) <: TableauFIRK
    @test typeof(TableauLobattoIIIC̄(4)) <: TableauFIRK
    @test typeof(TableauLobattoIIID(2)) <: TableauFIRK
    @test typeof(TableauLobattoIIID(3)) <: TableauFIRK
    @test typeof(TableauLobattoIIID(4)) <: TableauFIRK
    @test typeof( @test_logs (:warn, r"Initializing TableauFIRK with diagonally implicit tableau LobattoIIIE\(2\).*") TableauLobattoIIIE(2) ) <: TableauFIRK
    @test typeof(TableauLobattoIIIE(3)) <: TableauFIRK
    @test typeof(TableauLobattoIIIE(4)) <: TableauFIRK
    @test typeof(TableauLobattoIIIF(2)) <: TableauFIRK
    @test typeof(TableauLobattoIIIF(3)) <: TableauFIRK
    @test typeof(TableauLobattoIIIF(4)) <: TableauFIRK
    @test typeof(TableauLobattoIIIG(2)) <: TableauFIRK
    @test typeof(TableauLobattoIIIG(3)) <: TableauFIRK
    @test typeof(TableauLobattoIIIG(4)) <: TableauFIRK

    # instatiate all partitioned tableaus
    @test typeof(TableauSymplecticEulerA())     <: TableauEPRK
    @test typeof(TableauSymplecticEulerB())     <: TableauEPRK
    @test typeof(TableauLobattoIIIAIIIB2())     <: TableauEPRK
    @test typeof(TableauLobattoIIIBIIIA2())     <: TableauEPRK

    # test instatiation of partioned tableau by composition of two RK tableaus
    @test typeof(TableauEPRK(:PERK4, 4, TableauERK4().q, TableauERK4().q)) <: TableauEPRK
    @test TableauEPRK(:PERK4, 4, TableauERK4().q, TableauERK4().q) == TableauEPRK(:PERK4, 4, TableauERK4().q)


    # TODO Add tests for TableauIPRK, TableauSARK, TableauSPARK and TableauGLM.

end
