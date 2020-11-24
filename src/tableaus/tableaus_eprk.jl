
function CoefficientsSymplecticEulerForward()
    a = [[0.0 0.0]
         [1.0 0.0]]
    b = [1.0, 0.0]
    c = [0.0, 1.0]

    o = 1

    CoefficientsRK(:symplectic_euler_forward, o, a, b, c)
end

function CoefficientsSymplecticEulerBackward()
    a = [[0.0 0.0]
         [0.0 1.0]]
    b = [0.0, 1.0]
    c = [0.0, 1.0]

    o = 1

    CoefficientsRK(:symplectic_euler_backward, o, a, b, c)
end

"Tableau for symplectic Euler-A method"
function TableauSymplecticEulerA()
    TableauEPRK(:symplectic_euler_a, 1, CoefficientsSymplecticEulerForward(), CoefficientsSymplecticEulerBackward())
end

"Tableau for symplectic Euler-B method"
function TableauSymplecticEulerB()
    TableauEPRK(:symplectic_euler_b, 1, CoefficientsSymplecticEulerBackward(), CoefficientsSymplecticEulerForward())
end


"Tableau for Gauss-Lobatto IIIAIIIB method with s=2 stages"
function TableauLobattoIIIAIIIB2()
    TableauEPRK(:lobatto_IIIA_IIIB_2, 2, CoefficientsLobattoIIIA(2), CoefficientsLobattoIIIB(2))
end


"Tableau for Gauss-Lobatto IIIBIIIA method with s=2 stages"
function TableauLobattoIIIBIIIA2()
    TableauEPRK(:lobatto_IIIB_IIIA_2, 2, CoefficientsLobattoIIIB(2), CoefficientsLobattoIIIA(2))
end

