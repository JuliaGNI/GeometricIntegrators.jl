
"Tableau for explicit Euler method"
function TableauExplicitEuler()
    a = zeros(Float64, 1, 1)
    b = [1.0]
    c = [0.0]
    o = 1

    TableauERK(:explicit_euler, o, a, b, c)
end

"Tableau for explicit midpoint method"
function TableauExplicitMidpoint()
    a = [[0.0 0.0]
         [0.5 0.0]]
    b = [0.0, 1.0]
    c = [0.0, 0.5]
    o = 2

    TableauERK(:explicit_midpoint, o, a, b, c)
end

"Tableau for Runge's method"
function TableauRunge()
    a = [[0.0 0.0]
         [1.0 0.0]]
    b = [0.5, 0.5]
    c = [0.0, 1.0]
    o = 2

    TableauERK(:runge, o, a, b, c)
end

"Tableau for Heun's method"
function TableauHeun()
    a = [[0.0 0.0]
         [1.0 0.0]]
    b = [0.5, 0.5]
    c = [0.0, 1.0]
    o = 2

    TableauERK(:heun, o, a, b, c)
end

"Tableau for Kutta's method of order three"
function TableauKutta()
    a = [[ 0.0 0.0 0.0]
         [ 0.5 0.0 0.0]
         [-1.0 2.0 0.0]]
    b = [1/6, 4/6, 1/6]
    c = [0.0, 0.5, 1.0]
    o = 3

    TableauERK(:kutta, o, a, b, c)
end

"Tableau for explicit Runge-Kutta method of order four (1/6 rule)"
function TableauERK416()
    a = [[0.0 0.0 0.0 0.0]
         [0.5 0.0 0.0 0.0]
         [0.0 0.5 0.0 0.0]
         [0.0 0.0 1.0 0.0]]
    b = [1/6, 1/3, 1/3, 1/6]
    c = [0.0, 0.5, 0.5, 1.0]
    o = 4

    TableauERK(:erk4, o, a, b, c)
end

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
