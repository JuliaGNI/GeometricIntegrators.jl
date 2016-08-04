
using GeomDAE
using Base.Test

# explicit midpoint
a = Rational{Int64}[[0    0]
                    [1//2 0]]
b = Rational{Int64}[0, 1   ]
c = Rational{Int64}[0, 1//2]

tab_explicit_midpoint = TableauERK(:explicit_midpoint, 2, a, b, c)
tab_explicit_midpoint = TableauIRK(:explicit_midpoint, 2, a, b, c)
tab_explicit_midpoint = TableauNLIRK(:explicit_midpoint, 2, a, b, c)


# Heun
a = [[0.0 0.0]
     [1.0 0.0]]
b = [0.5, 0.5]
c = [0.0, 1.0]

tab_explicit_heun = TableauERK(:heun, 2, a, b, c)
tab_explicit_heun = TableauIRK(:heun, 2, a, b, c)
tab_explicit_heun = TableauNLIRK(:heun, 2, a, b, c)


# Crouzeix
fac = 0.5/√3
a = [[ 0.5+fac 0.0    ]
     [-2.0*fac 0.5+fac]]
b = [0.5,     0.5    ]
c = [0.5+fac, 0.5-fac]

@test_throws AssertionError tab_explicit_crouzeix = TableauERK(:crouzeix, 2, a, b, c)
tab_explicit_crouzeix = TableauIRK(:crouzeix, 2, a, b, c)
tab_explicit_crouzeix = TableauNLIRK(:crouzeix, 2, a, b, c)


# implicit midpoint
a = reshape(Rational{Int64}[1//2], 1, 1) # TODO How to create a 2dim 1x1 array?
b = Rational{Int64}[1   ]
c = Rational{Int64}[1//2]

@test_throws AssertionError tab_implicit_midpoint = TableauERK(:implicit_midpoint, 2, a, b, c)
tab_implicit_midpoint = TableauIRK(:implicit_midpoint, 2, a, b, c)
tab_implicit_midpoint = TableauNLIRK(:implicit_midpoint, 2, a, b, c)


# Gauss-Legendre
fac = √3/6.
a = [[0.25     0.25-fac]
     [0.25+fac 0.25    ]]
b = [0.5,     0.5    ]
c = [0.5-fac, 0.5+fac]

@test_throws AssertionError tab_explicit_glrk2 = TableauERK(:glrk2, 4, a, b, c)
@test_throws AssertionError tab_explicit_glrk2 = TableauIRK(:glrk2, 4, a, b, c)
tab_explicit_glrk2 = TableauNLIRK(:glrk2, 4, a, b, c)
