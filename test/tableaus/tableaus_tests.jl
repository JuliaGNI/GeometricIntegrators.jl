
# explicit midpoint
a = Rational{Int64}[[0    0]
                    [1//2 0]]
b = Rational{Int64}[0, 1   ]
c = Rational{Int64}[0, 1//2]
o = 2

@test TableauERK(:explicit_midpoint, o, a, b, c) == getTableauExplicitMidpoint()

tab_explicit_midpoint1 = TableauERK(:explicit_midpoint, o, a, b, c)
tab_explicit_midpoint2 = TableauDIRK(:explicit_midpoint, o, a, b, c)
tab_explicit_midpoint3 = TableauFIRK(:explicit_midpoint, o, a, b, c)

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

tab_explicit_heun = TableauERK(:heun, o, a, b, c)
tab_explicit_heun = TableauDIRK(:heun, o, a, b, c)
tab_explicit_heun = TableauFIRK(:heun, o, a, b, c)

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
tab_explicit_crouzeix = TableauDIRK(:crouzeix, o, a, b, c)
tab_explicit_crouzeix = TableauFIRK(:crouzeix, o, a, b, c)


# implicit midpoint
a = ones(Rational{Int64}, 1, 1) // 2
b = Rational{Int64}[1   ]
c = Rational{Int64}[1//2]
o = 2

@test TableauFIRK(:implicit_midpoint, o, a, b, c) == getTableauImplicitMidpoint()
@test TableauFIRK(:implicit_midpoint, o, a, b, c) == getTableauGLRK1()

@test_throws AssertionError tab_implicit_midpoint = TableauERK(:implicit_midpoint, o, a, b, c)
@test_throws AssertionError tab_implicit_midpoint = TableauDIRK(:implicit_midpoint, o, a, b, c)
tab_implicit_midpoint = TableauFIRK(:implicit_midpoint, o, a, b, c)


# Gauss-Legendre
a = [[0.25       0.25-√3/6]
     [0.25+√3/6  0.25     ]]
b = [0.5,      0.5     ]
c = [0.5-√3/6, 0.5+√3/6]
o = 4

@test TableauFIRK(:glrk2, o, a, b, c) == getTableauGLRK2()

@test_throws AssertionError tab_explicit_glrk2 = TableauERK(:glrk2, o, a, b, c)
@test_throws AssertionError tab_explicit_glrk2 = TableauDIRK(:glrk2, o, a, b, c)
tab_explicit_glrk2 = TableauFIRK(:glrk2, o, a, b, c)


# instantiate all explicit tableaus
@test typeof(getTableauExplicitEuler()) <: TableauERK
@test typeof(getTableauExplicitMidpoint()) <: TableauERK
@test typeof(getTableauHeun()) <: TableauERK
@test typeof(getTableauKutta()) <: TableauERK
@test typeof(getTableauERK4()) <: TableauERK
@test typeof(getTableauERK438()) <: TableauERK

# instantiate all diagonally implicit tableaus
@test typeof(getTableauCrouzeix()) <: TableauDIRK

# instantiate all fully implicit tableaus
@test typeof(getTableauImplicitEuler()) <: TableauFIRK
@test typeof(getTableauImplicitMidpoint()) <: TableauFIRK
@test typeof(getTableauGLRK1()) <: TableauFIRK
@test typeof(getTableauGLRK2()) <: TableauFIRK
@test typeof(getTableauGLRK3()) <: TableauFIRK

# instatiate all partitioned tableaus
@test typeof(getTableauSymplecticEulerA()) <: TableauEPRK
@test typeof(getTableauSymplecticEulerB()) <: TableauEPRK

# test instatiation of partioned tableau by composition of two RK tableaus
@test typeof(TableauEPRK(:PERK4, 4, getTableauERK4().q, getTableauERK4().q)) <: TableauEPRK

# test computation of Gauss-Legendre Runge-Kutta tableaus
@test getTableauGLRK(1) == getTableauGLRK1()

glrk2_tab1 = getTableauGLRK(2)
glrk2_tab2 = getTableauGLRK2()
@test glrk2_tab1.q.a ≈ glrk2_tab1.q.a atol=2 * eps()
@test glrk2_tab1.q.b == glrk2_tab1.q.b
@test glrk2_tab1.q.c == glrk2_tab1.q.c

glrk3_tab1 = getTableauGLRK(3)
glrk3_tab2 = getTableauGLRK3()
@test glrk3_tab1.q.a ≈ glrk3_tab1.q.a atol=2 * eps()
@test glrk3_tab1.q.b == glrk3_tab1.q.b
@test glrk3_tab1.q.c == glrk3_tab1.q.c


# TODO Add tests for TableauIPRK, TableauSARK, TableauSPARK and TableauGLM.
