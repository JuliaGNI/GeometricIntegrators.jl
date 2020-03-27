
@testset "$(rpad("Runge-Kutta Tableaus",80))" begin

    # explicit midpoint
    a = Rational{Int64}[[0    0]
                        [1//2 0]]
    b = Rational{Int64}[0, 1   ]
    c = Rational{Int64}[0, 1//2]
    o = 2

    @test TableauERK(:explicit_midpoint, o, a, b, c) == getTableauExplicitMidpoint()

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

    @test TableauERK(:heun, o, a, b, c) == getTableauHeun()

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
    o = 2

    @test TableauDIRK(:crouzeix, o, a, b, c) == getTableauCrouzeix()

    @test_throws AssertionError tab_explicit_crouzeix = TableauERK(:crouzeix, o, a, b, c)
    @test_nowarn tab_explicit_crouzeix = TableauDIRK(:crouzeix, o, a, b, c)
    @test_logs (:warn, r"Initializing TableauFIRK with diagonally implicit tableau crouzeix.*") tab_explicit_crouzeix = TableauFIRK(:crouzeix, o, a, b, c)


    # implicit midpoint
    a = ones(Rational{Int64}, 1, 1) // 2
    b = Rational{Int64}[1   ]
    c = Rational{Int64}[1//2]
    o = 2

    @test TableauFIRK(:implicit_midpoint, o, a, b, c) == getTableauImplicitMidpoint()
    @test TableauFIRK(:implicit_midpoint, o, a, b, c) == getTableauGLRK(1)

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

    @test TableauFIRK(:glrk2, o, a, b, c) == getTableauGLRK(2)

    @test_throws AssertionError tab_explicit_glrk2 = TableauERK(:glrk2, o, a, b, c)
    @test_throws AssertionError tab_explicit_glrk2 = TableauDIRK(:glrk2, o, a, b, c)
    tab_explicit_glrk2 = TableauFIRK(:glrk2, o, a, b, c)


    # instantiate all explicit tableaus
    @test typeof(getTableauExplicitEuler()) <: TableauERK
    @test typeof(getTableauExplicitMidpoint()) <: TableauERK
    @test typeof(getTableauHeun()) <: TableauERK
    @test typeof(getTableauKutta()) <: TableauERK
    @test typeof(getTableauRunge()) <: TableauERK
    @test typeof(getTableauERK4()) <: TableauERK
    @test typeof(getTableauERK438()) <: TableauERK

    # instantiate all diagonally implicit tableaus
    @test typeof(getTableauCrouzeix()) <: TableauDIRK

    # instantiate all fully implicit tableaus
    @test typeof(getTableauImplicitEuler()) <: TableauFIRK
    @test typeof(getTableauImplicitMidpoint()) <: TableauFIRK
    @test typeof(getTableauGLRK(1)) <: TableauFIRK
    @test typeof(getTableauGLRK(2)) <: TableauFIRK
    @test typeof(getTableauSRK3()) <: TableauFIRK

    @test typeof( @test_logs (:warn, r"Initializing TableauFIRK with diagonally implicit tableau LobIIIA2.*") getTableauLobIIIA2() ) <: TableauFIRK
    @test typeof(getTableauLobIIIA3()) <: TableauFIRK
    @test typeof(getTableauLobIIIA4()) <: TableauFIRK
    @test typeof( @test_logs (:warn, r"Initializing TableauFIRK with explicit tableau LobIIIB2.*") getTableauLobIIIB2() ) <: TableauFIRK
    @test typeof(getTableauLobIIIB3()) <: TableauFIRK
    @test typeof(getTableauLobIIIB4()) <: TableauFIRK
    @test typeof(getTableauLobIIIC2()) <: TableauFIRK
    @test typeof(getTableauLobIIIC3()) <: TableauFIRK
    @test typeof(getTableauLobIIIC4()) <: TableauFIRK
    @test typeof(getTableauLobIIID2()) <: TableauFIRK
    @test typeof(getTableauLobIIID3()) <: TableauFIRK
    @test typeof(getTableauLobIIID4()) <: TableauFIRK
    @test typeof( @test_logs (:warn, r"Initializing TableauFIRK with diagonally implicit tableau LobIIIE2.*") getTableauLobIIIE2() ) <: TableauFIRK
    @test typeof(getTableauLobIIIE3()) <: TableauFIRK
    @test typeof(getTableauLobIIIE4()) <: TableauFIRK
    @test typeof(getTableauLobIIIF2()) <: TableauFIRK
    @test typeof(getTableauLobIIIF3()) <: TableauFIRK
    @test typeof(getTableauLobIIIF4()) <: TableauFIRK
    @test typeof(getTableauLobIIIG2()) <: TableauFIRK
    @test typeof(getTableauLobIIIG3()) <: TableauFIRK
    @test typeof(getTableauLobIIIG4()) <: TableauFIRK

    @test typeof(getTableauRadIIA2()) <: TableauFIRK
    @test typeof(getTableauRadIIA3()) <: TableauFIRK

    # instatiate all partitioned tableaus
    @test typeof(getTableauSymplecticEulerA()) <: TableauEPRK
    @test typeof(getTableauSymplecticEulerB()) <: TableauEPRK
    @test typeof(getTableauLobattoIIIAIIIB2()) <: TableauEPRK
    @test typeof(getTableauLobattoIIIBIIIA2()) <: TableauEPRK

    # test instatiation of partioned tableau by composition of two RK tableaus
    @test typeof(TableauEPRK(:PERK4, 4, getTableauERK4().q, getTableauERK4().q)) <: TableauEPRK
    @test TableauEPRK(:PERK4, 4, getTableauERK4().q, getTableauERK4().q) == TableauEPRK(:PERK4, 4, getTableauERK4().q)


    # TODO Add tests for TableauIPRK, TableauSARK, TableauSPARK and TableauGLM.

end
