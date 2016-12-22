

function getTableauSymplecticEulerForward()
    a = [[0.0 0.0]
         [1.0 0.0]]
    b = [1.0, 0.0]
    c = [0.0, 1.0]

    o = 1

    TableauRK(:symplectic_euler_backward, o, a, b, c)
end

function getTableauSymplecticEulerBackward()
    a = [[0.0 0.0]
         [0.0 1.0]]
    b = [0.0, 1.0]
    c = [0.0, 1.0]

    o = 1

    TableauRK(:symplectic_euler_forward, o, a, b, c)
end

"Tableau for symplectic Euler-A method"
function getTableauSymplecticEulerA()
    TableauEPRK(:symplectic_euler_a, 1, getTableauSymplecticEulerForward(), getTableauSymplecticEulerBackward())
end

"Tableau for symplectic Euler-B method"
function getTableauSymplecticEulerB()
    TableauEPRK(:symplectic_euler_b, 1, getTableauSymplecticEulerBackward(), getTableauSymplecticEulerForward())
end
