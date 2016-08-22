
"Tableau for symplectic Euler-A method"
function getTableauSymplecticEulerA()
    a_q = [[0.0 0.0]
           [1.0 0.0]]
    b_q = [1.0, 0.0]
    c_q = [0.0, 1.0]

    a_p = [[0.0 0.0]
           [0.0 1.0]]
    b_p = [0.0, 1.0]
    c_p = [0.0, 1.0]

    o = 1

    TableauPRK(:symplectic_euler_a, o, a_q, a_p, b_q, b_p, c_q, c_p)
end

"Tableau for symplectic Euler-B method"
function getTableauSymplecticEulerB()
    a_q = [[0.0 0.0]
           [0.0 1.0]]
    b_q = [1.0, 0.0]
    c_q = [0.0, 1.0]

    a_p = [[0.0 0.0]
           [1.0 0.0]]
    b_p = [0.0, 1.0]
    c_p = [0.0, 1.0]

    o = 1

    TableauPRK(:symplectic_euler_b, o, a_q, a_p, b_q, b_p, c_q, c_p)
end
