
using GeomDAE
using Base.Test

# explicit midpoint
a = Rational{Int64}[[0    0]
                    [1//2 0]]
b = Rational{Int64}[0, 1   ]
c = Rational{Int64}[0, 1//2]
o = 2

@test TableauERK(:explicit_midpoint, o, a, b, c) == getTableauExplicitMidpoint()

tab_explicit_midpoint1 = TableauERK(:explicit_midpoint, o, a, b, c)
tab_explicit_midpoint2 = TableauIRK(:explicit_midpoint, o, a, b, c)
tab_explicit_midpoint3 = TableauNLIRK(:explicit_midpoint, o, a, b, c)

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
tab_explicit_heun = TableauIRK(:heun, o, a, b, c)
tab_explicit_heun = TableauNLIRK(:heun, o, a, b, c)

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

@test TableauIRK(:crouzeix, o, a, b, c) == getTableauCrouzeix()

@test_throws AssertionError tab_explicit_crouzeix = TableauERK(:crouzeix, o, a, b, c)
tab_explicit_crouzeix = TableauIRK(:crouzeix, o, a, b, c)
tab_explicit_crouzeix = TableauNLIRK(:crouzeix, o, a, b, c)


# implicit midpoint
a = ones(Rational{Int64}, 1, 1) // 2
b = Rational{Int64}[1   ]
c = Rational{Int64}[1//2]
o = 2

@test TableauNLIRK(:implicit_midpoint, o, a, b, c) == getTableauImplicitMidpoint()
@test TableauNLIRK(:implicit_midpoint, o, a, b, c) == getTableauGLRK1()

@test_throws AssertionError tab_implicit_midpoint = TableauERK(:implicit_midpoint, o, a, b, c)
@test_throws AssertionError tab_implicit_midpoint = TableauIRK(:implicit_midpoint, o, a, b, c)
tab_implicit_midpoint = TableauNLIRK(:implicit_midpoint, o, a, b, c)


# Gauss-Legendre
a = [[0.25       0.25-√3/6]
     [0.25+√3/6  0.25     ]]
b = [0.5,      0.5     ]
c = [0.5-√3/6, 0.5+√3/6]
o = 4

@test TableauNLIRK(:glrk2, o, a, b, c) == getTableauGLRK2()

@test_throws AssertionError tab_explicit_glrk2 = TableauERK(:glrk2, o, a, b, c)
@test_throws AssertionError tab_explicit_glrk2 = TableauIRK(:glrk2, o, a, b, c)
tab_explicit_glrk2 = TableauNLIRK(:glrk2, o, a, b, c)


# instantiate all explicit tableaus
@test typeof(getTableauExplicitEuler()) <: TableauERK
@test typeof(getTableauExplicitMidpoint()) <: TableauERK
@test typeof(getTableauHeun()) <: TableauERK
@test typeof(getTableauKutta()) <: TableauERK
@test typeof(getTableauERK4()) <: TableauERK
@test typeof(getTableauERK438()) <: TableauERK

# instantiate all diagonally implicit tableaus
@test typeof(getTableauCrouzeix()) <: TableauIRK

# instantiate all fully implicit tableaus
@test typeof(getTableauImplicitEuler()) <: TableauNLIRK
@test typeof(getTableauImplicitMidpoint()) <: TableauNLIRK
@test typeof(getTableauGLRK1()) <: TableauNLIRK
@test typeof(getTableauGLRK2()) <: TableauNLIRK
@test typeof(getTableauGLRK3()) <: TableauNLIRK


# TODO Add tests for TableauPRK, TableauSPARK and TableauGLM.
