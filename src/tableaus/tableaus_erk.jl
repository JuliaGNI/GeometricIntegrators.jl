
"Tableau for explicit Euler method"
function TableauExplicitEuler()
    a = zeros(Float64, 1, 1)
    b = [1.0]
    c = [0.0]
    o = 1

    TableauERK(:explicit_euler, o, a, b, c)
end

"Alias for [`TableauExplicitEuler`](@ref)"
TableauForwardEuler = TableauExplicitEuler

"Tableau for explicit midpoint method"
function TableauExplicitMidpoint()
    a = [[0.0  0.0]
         [0.5  0.0]]
    b =  [0.0, 1.0]
    c =  [0.0, 0.5]
    o = 2

    TableauERK(:explicit_midpoint, o, a, b, c)
end

"Tableau of Heun's two-stage, 2nd order method"
function TableauHeun2(T=Float64)
    a = [[ 0.0   0.0 ]
         [ 1.0   0.0 ]]
    b =  [ 0.5,  0.5 ]
    c =  [ 0.0,  1.0 ]
    o = 2

    TableauERK(:heun2, o, a, b, c)
end

"Tableau of Heun's three-stage, 3rd order method"
function TableauHeun3(T=Float64)
    a = [[ 0.0   0.0   0.0 ]
         [ 1/3   0.0   0.0 ]
         [ 0.0   2/3   0.0 ]]
    b =  [ 1/4,  0.0,  3/4 ]
    c =  [ 0.0,  1/3,  2/3 ]
    o = 3

    TableauERK(:heun3, o, a, b, c)
end

"Tableau of Ralston's two-stage, 2nd order method"
function TableauRalston2(T=Float64)
    a = [[ 0.0   0.0 ]
         [ 2/3   0.0 ]]
    b =  [ 1/4,  3/4 ]
    c =  [ 0.0,  2/3 ]
    o = 2

    TableauERK(:ralston2, o, a, b, c)
end

"Tableau of Ralston's three-stage, 3rd order method"
function TableauRalston3(T=Float64)
    a = [[ 0.0  0.0  0.0 ]
         [ 1/2  0.0  0.0 ]
         [ 0.0  3/4  0.0 ]]
    b =  [ 2/9, 3/9, 4/9 ]
    c =  [ 0.0, 1/2, 3/4 ]
    o = 3

    TableauERK(:ralston3, o, a, b, c)
end

"Tableau for Runge's method"
function TableauRunge()
    a = [[0.0  0.0]
         [1.0  0.0]]
    b =  [0.5, 0.5]
    c =  [0.0, 1.0]
    o = 2

    TableauERK(:runge, o, a, b, c)
end

"Alias for [`TableauRunge`](@ref)"
TableauRunge2 = TableauRunge

"Tableau for Kutta's method of order three"
function TableauKutta()
    a = [[ 0.0  0.0  0.0]
         [ 0.5  0.0  0.0]
         [-1.0  2.0  0.0]]
    b =  [ 1/6, 4/6, 1/6]
    c =  [ 0.0, 0.5, 1.0]
    o = 3

    TableauERK(:kutta, o, a, b, c)
end

"Alias for [`TableauKutta`](@ref)"
TableauKutta3 = TableauKutta

"Tableau for explicit Runge-Kutta method of order four (1/6 rule)"
function TableauERK416()
    a = [[0.0  0.0  0.0  0.0]
         [0.5  0.0  0.0  0.0]
         [0.0  0.5  0.0  0.0]
         [0.0  0.0  1.0  0.0]]
    b =  [1/6, 1/3, 1/3, 1/6]
    c =  [0.0, 0.5, 0.5, 1.0]
    o = 4

    TableauERK(:erk4, o, a, b, c)
end

"Alias for [`TableauRK416`](@ref)"
TableauERK4 = TableauERK416

"Tableau for explicit Runge-Kutta method of order four (3/8 rule)"
function TableauERK438()
    a = [[ 0.0  0.0  0.0  0.0]
         [ 1/3  0.0  0.0  0.0]
         [-1/3  1.0  0.0  0.0]
         [ 1.0 -1.0  1.0  0.0]]
    b = [1/8, 3/8, 3/8, 1/8]
    c = [0.0, 1/3, 2/3, 1.0]
    o = 4

    TableauERK(:erk438, o, a, b, c)
end
