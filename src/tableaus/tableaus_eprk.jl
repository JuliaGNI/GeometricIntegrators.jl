
function getCoefficientsSymplecticEulerForward()
    a = [[0.0 0.0]
         [1.0 0.0]]
    b = [1.0, 0.0]
    c = [0.0, 1.0]

    o = 1

    CoefficientsRK(:symplectic_euler_forward, o, a, b, c)
end

function getCoefficientsSymplecticEulerBackward()
    a = [[0.0 0.0]
         [0.0 1.0]]
    b = [0.0, 1.0]
    c = [0.0, 1.0]

    o = 1

    CoefficientsRK(:symplectic_euler_backward, o, a, b, c)
end

"Tableau for symplectic Euler-A method"
function getTableauSymplecticEulerA()
    TableauEPRK(:symplectic_euler_a, 1, getCoefficientsSymplecticEulerForward(), getCoefficientsSymplecticEulerBackward())
end

"Tableau for symplectic Euler-B method"
function getTableauSymplecticEulerB()
    TableauEPRK(:symplectic_euler_b, 1, getCoefficientsSymplecticEulerBackward(), getCoefficientsSymplecticEulerForward())
end


"Tableau for Gauss-Lobatto IIIAIIIB method with s=2 stages"
function getTableauLobattoIIIAIIIB2()
    TableauEPRK(:lobatto_IIIA_IIIB_2, 2, getCoefficientsLobIIIA2(), getCoefficientsLobIIIB2())
end


"Tableau for Gauss-Lobatto IIIBIIIA method with s=2 stages"
function getTableauLobattoIIIBIIIA2()
    TableauEPRK(:lobatto_IIIB_IIIA_2, 2, getCoefficientsLobIIIB2(), getCoefficientsLobIIIA2())
end

